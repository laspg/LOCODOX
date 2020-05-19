% DOXY_argo_prepare prepare the argo float data for Doxy correction.
%
% SYNTAX
% [argoWork,Work] = DOXY_argo_prepare(argo,Work,LogDir,fig_nbr,argoP, argoWorkP)
%
% DESCRIPTION
% DOXY_argo_prepare prepare the argo float data for Doxy correction, through
% different steps :
% * Getting the right PTS (from doxy sensor or not)
% * Remove negative and constant pressure
% * Compute Depth and Density
% * Apply quality flags to DOXY
% * Convert Doxy data to new unit
% 
%
% INPUT
%     argo (structure)      Argo float structure (get directly from data).
%                           Example :
%                               argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%
%                               argo.pres=
%                                       name: 'DOXY'
%                                        dim: {'N_PROF'  'N_LEVELS'}
%                                       data: [85x120 single]
%                                  long_name: 'DISSOLVED OXYGEN'
%                              standard_name: 'moles_of_oxygen_per_unit_mass_in_sea_water'
%                                 FillValue_: 99999
%                                      units: 'micromole/kg'
%                                  valid_min: 0
%                                  valid_max: 650
%                                   C_format: '%9.3f'
%                             FORTRAN_format: 'F9.3'
%                                 resolution: 0.0010
%                                       type: 5
%
%     argoP (structure)     Argo float structure (get directly from data)
%                           for primary sampling
%
%     Work (structure)      Doxy correction working structure, issued and
%                           computed from argo float data.
%                           Example:
%                             Work = 
%                                   readme: [1x1 struct]
%                                     unit: 'mumol/kg'
%                                      wmo: 1901205
%                                 DOXY_RAW: [85x120 single]
%                                  CAPTEUR: {'Optode'}
%                                  DOXY_QC: [85x120 char]
%
%     LogDir (string)       the full path of the Log directory (defined in the
%                           configuration file)
%
%     fig_nbr               Number of the figure to plot (for title and
%                           position)
%     
%     argoWorkP (structure) float structure on the argo model, with the PTS
%                           adjusted field
%                           Example :
%                             argoWorkP =
%                                 pres_adjusted: [1x1 struct]
%                                 temp_adjusted: [1x1 struct]
%                                 psal_adjusted: [1x1 struct]
%
%
% OUTPUT
%     argowork (structure)  Float working structure, issued and computed
%                           from argo float data.Modified and
%                           completed
%                           Example:
%                           argoWork = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%
%     Work (structure)      Doxy correction working structure, issued and
%                           computed from the float data. Modified and
%                           completed
%
% CALL : 
%   DOXY_get_PTS, DOXY_PLOT_QC, prsenz, sw_pden, DOXY_QC, DOXY_convert,
%   sw_satO2, 02cto02s
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 26/11/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%              v3    26/01/2017  Anne Piron, Altran Ouest
%                                argo - ameliorations 2017 phase 1
%              v3.2  13/09/2017  Emilie Brion, Altran Ouest
%                                - argoWork.doxy_adjusted is no more set to
%                                NaN when argoWork.an_dens is NaN. It's QC
%                                is consequently not set to 4 for the same
%                                reason.
%                                - Work.DEPTH computed from non NaN values
%                                of pres_adjusted and an_dens
%              v3.3  19/09/2017  Emilie Brion, Altran Ouest
%                                Add the call to DOXY_unselect    
%              v3.3  25/10/2017  Emilie Brion, Altran Ouest
%                                Take into account that Primary Sampling
%                                could be inexistent (if the float doesn't
%                                go deeper than the cutoff pressure)
%              v3.0  20/02/2019  Marine GALLIAN, Altran Ouest
%                                Delete the hook removal
%                                Change figure showing data and their QC
%                                Do not plot data when all is NaN (for BGC
%                                floats)
%                                Introduction of the variable "fig_nbr" to
%                                document which sampling is displayed

function [argoWork,Work] = DOXY_argo_prepare(argo,Work,LogDir,fig_nbr,argoP, argoWorkP)
% =========================================================================
%% Initialisation
% =========================================================================
wmo        = Work.wmo;
unit       = Work.unit;

%Figure name and position
fig_name={'DATA PRESENTATION : PRIMARY SAMPLING','DATA PRESENTATION : SECONDARY SAMPLING','DATA PRESENTATION : NEAR SURFACE SAMPLING','DATA PRESENTATION : OTHER SAMPLING'};
fig_pos={[0 0.6 0.45 0.4],[0 0.6 0.45 0.4],[0 0.25 0.45 0.4],[0 0.25 0.45 0.4]};

% Doxy raw, if all NaN do not plot data %--MG
Work.DOXY_RAW = argo.doxy.data;
if all(all(isnan(argo.doxy.data)))
    if Work.makePlot==1
        Work.makePlot=0;
        mem_plot=1;
    else
        mem_plot=0;
    end
end

% Calcul du jour de l'annee (1-365) en fonction du jour julien
Work.timar = datetime(datevec(argo.juld.data + datenum(1950,1,1,0,0,0)));
tabyear = datestr(argo.juld.data + datenum(1950,1,1),10);
Work.datat = (argo.juld.data + datenum(1950,1,1)) - datenum(str2num(tabyear),1,1)+1; %#ok<ST2NM>

% =========================================================================
%% Get the suitable PTS data (from doxy sensor, or from adjusted, or normal)
% =========================================================================
[argoWork,Work] = DOXY_get_PTS(argo,Work,LogDir);

if Work.savePlot == 1
    Work.savePlot = 0;
    mem = 1;
else
    mem = 0;
end

% Control Plot
if Work.makePlot
    hfig = figure('unit','normalized','Outerposition',fig_pos{fig_nbr},...
        'Name',sprintf('%s - %d',fig_name{fig_nbr},Work.wmo),'NumberTitle','off');
end

% =========================================================================
%% Compute raw density
% If unpumped psal_adjusted value has only NaN values (so QC = 4), get the
% psal field instead. In delayed mode, all the QC of the unpumped
% psal_ajusted takes QC = 4 and data = NaN, because of the unpumped state.
% But for DOXY correction, we does use the real time psal (and temp, pres),
% with non NaN values du to an artificial flag 4.
% =========================================================================

if (all(all(isnan(argo.psal.data))) || all(all(isnan(argo.temp.data))))
    isok = ismember(argoP.cycle_number.data,argo.cycle_number.data);    
    Work.DENS = real(sw_pden(argoP.psal.data(isok,:), ...
        argoP.temp.data(isok,:), argoP.pres.data(isok,:),0));
    argoWork.raw.an_dens.data = Work.DENS - 1000;
    argoWork.raw.density.data = Work.DENS;
else
    Work.DENS = real(sw_pden(argo.psal.data, ...
        argo.temp.data, argo.pres.data,0));
    argoWork.raw.an_dens.data = Work.DENS - 1000;
    argoWork.raw.density.data = Work.DENS;
end

% =========================================================================
%% Compute Initial Saturation and Convert Argo Doxy data in the choosen unit
%% and keep the doxy raw field
% =========================================================================
% Saturation and saturation percentage
%   Convert in ml/L to compute the % of saturation (%SAT)
unit = argo.doxy.units;
tmpDoxy = DOXY_convert(argo.doxy.data,unit,'mumol/L',argoWork.raw.an_dens.data);
if (all(all(isnan(argo.psal.data))) || all(all(isnan(argo.temp.data))))
    isok = ismember(argoP.cycle_number.data,argo.cycle_number.data);
    if Work.presEff == 1, P = argoP.pres.data(isok,:); else, P = 0; end
    argoWork.raw.psat.data = O2ctoO2s(tmpDoxy,argoP.temp.data(isok,:),argoP.psal.data(isok,:),P); % in mL/L
else
    if Work.presEff == 1, P = argo.pres.data; else, P = 0; end
    argoWork.raw.psat.data = O2ctoO2s(tmpDoxy,argo.temp.data,argo.psal.data,P); % in mL/L
end
tmpDoxy = DOXY_convert(argo.doxy.data,unit,'ml/L',argoWork.raw.an_dens.data);
argoWork.raw.sat.data = (tmpDoxy./argoWork.raw.psat.data)*100;


% =========================================================================
%% Remove negative and constant pressure
% =========================================================================
% No negative pressure
inok = argoWork.pres_adjusted.data < 0;
argoWork.pres_adjusted.data(inok) = NaN;

% No constant pressure
myDim = argoWork.pres_adjusted.dim;
okDim = find(strcmp(myDim,'N_PROF'));
nbprof = size(argoWork.pres_adjusted.data,okDim);
inok = false(size(argoWork.pres_adjusted.data));
inok(1:nbprof,1:end-1) = (diff(argoWork.pres_adjusted.data(1:nbprof,1:end),1,2) == 0);
argoWork.pres_adjusted.data(inok) = NaN;
inok(1:nbprof,2:end) = (diff(argoWork.pres_adjusted.data(1:nbprof,end:-1:1),1,2) == 0);
argoWork.pres_adjusted.data(inok) = NaN;

% =========================================================================
%% Interpolate temp and psal to secondary pres if VSS = secondary sampling    %%%%MARINE %18/06/19
% =========================================================================
if strcmp(argo.VSS,'Secondary sampling')
    if all(all(isnan(argoWork.temp_adjusted.data))) || all(all(isnan(argoWork.psal_adjusted.data)))
        i=1;
        %maximum in dbar for extrapolation of psal and temp from primary profile on levels of secondary (near surface and
        %at the bottom)
        interp_surf_max=5;
        interp_prof_max=5;
        while i <= argo.n_prof %Changement pour un while %marine 17/06/19
            isokP = argoP.cycle_number.data == argo.cycle_number.data(i);
            %marine 17/06/19
            if length(find(isokP))==1
                isokNoNan = ~isnan(argoWorkP.pres_adjusted.data(isokP,:));
                if any(find(isokNoNan))
                    X = argoWorkP.pres_adjusted.data(isokP,isokNoNan);
                    Y = argoWorkP.temp_adjusted.data(i,isokNoNan); %marine 11/06/19
                    YY = argoWorkP.psal_adjusted.data(i,isokNoNan); %marine 11/06/19
                    Xq = argoWork.pres_adjusted.data(i,:);
                    %marine 25/06/19 : suppression des doublons pour
                    %l'interpolation
                    [~,X_ind]=unique(X); X_ind=sort(X_ind);
                    X_tmp=X(X_ind); Y_tmp=Y(X_ind); YY_tmp=YY(X_ind);
                    [~,Y_ind]=unique(Y_tmp); Y_ind=sort(Y_ind);
                    [~,YY_ind]=unique(YY_tmp); YY_ind=sort(YY_ind);
                    argoWork.temp_adjusted.data(i,:) = interp1(X_tmp(Y_ind),Y_tmp(Y_ind),Xq);
                    argoWork.psal_adjusted.data(i,:) = interp1(X_tmp(YY_ind),YY_tmp(YY_ind),Xq);
                    
                    %Valeur inférieur au minimum du primary D
                    is_inf_secondary=Xq<min(X);
                    ind_diff_max=abs(Xq-min(X))<interp_surf_max;
                    [~,ind_min_P]=min(X);
                    argoWork.temp_adjusted.data(i,and(ind_diff_max,is_inf_secondary))=Y(ind_min_P); %marine 25/06/19
                    argoWork.psal_adjusted.data(i,and(ind_diff_max,is_inf_secondary))=YY(ind_min_P); %marine 25/06/19
                    
                    %Valeur supérieur au maximum du primary D
                    is_sup_secondary=Xq>max(X);
                    ind_diff_max=abs(Xq-max(X))<interp_prof_max;
                    [~,ind_max_P]=max(X);
                    argoWork.temp_adjusted.data(i,and(ind_diff_max,is_sup_secondary))=Y(ind_max_P); %marine 25/06/19
                    argoWork.psal_adjusted.data(i,and(ind_diff_max,is_sup_secondary))=YY(ind_max_P); %marine 25/06/19
                end
                i = i+1;
            elseif length(find(isokP))==2 %Cas avec des profils descendants
                ind_isokP=find(isokP);
                isokP1=isokP; isokP1(ind_isokP(2))=0; %Indice du profil montant
                isokP2=isokP; isokP2(ind_isokP(1))=0; %Indice du profil descendant (D)
                
                %%%%%%%%%%%%%%%%%%%%%% Profil montant %%%%%%%%%%%%%%%%%%%%%
                isokNoNan = ~isnan(argoWorkP.pres_adjusted.data(isokP1,:));                
                if any(find(isokNoNan))
                    X = argoWorkP.pres_adjusted.data(isokP1,isokNoNan);
                    Y = argoWorkP.temp_adjusted.data(isokP1,isokNoNan); %marine 11/06/19
                    YY = argoWorkP.psal_adjusted.data(isokP1,isokNoNan); %marine 11/06/19
                    Xq = argoWork.pres_adjusted.data(isokP1,:);
                    
                    %marine 25/06/19 : suppression des doublons pour
                    %l'interpolation
                    [~,X_ind]=unique(X); X_ind=sort(X_ind);
                    X_tmp=X(X_ind); Y_tmp=Y(X_ind); YY_tmp=YY(X_ind);
                    [~,Y_ind]=unique(Y_tmp); Y_ind=sort(Y_ind);
                    [~,YY_ind]=unique(YY_tmp); YY_ind=sort(YY_ind);
                    argoWork.temp_adjusted.data(isokP1,:) = interp1(X_tmp(Y_ind),Y_tmp(Y_ind),Xq);
                    argoWork.psal_adjusted.data(isokP1,:) = interp1(X_tmp(YY_ind),YY_tmp(YY_ind),Xq);
                    
                    %Valeur inférieur au minimum du primary D
                    is_inf_secondary=Xq<min(X);
                    ind_diff_max=abs(Xq-min(X))<interp_surf_max;
                    [~,ind_min_P]=min(X);
                    argoWork.temp_adjusted.data(isokP1,and(ind_diff_max,is_inf_secondary))=Y(ind_min_P);
                    argoWork.psal_adjusted.data(isokP1,and(ind_diff_max,is_inf_secondary))=YY(ind_min_P);
                    
                    %Valeur supérieur au maximum du primary D
                    is_sup_secondary=Xq>max(X);
                    ind_diff_max=abs(Xq-max(X))<interp_prof_max;
                    [~,ind_max_P]=max(X);
                    argoWork.temp_adjusted.data(isokP1,and(ind_diff_max,is_sup_secondary))=Y(ind_max_P);
                    argoWork.psal_adjusted.data(isokP1,and(ind_diff_max,is_sup_secondary))=YY(ind_max_P);
                end
                i = i+sum(isokP1);
                
                
                %%%%%%%%%%%%%%%%%%% Profil descendant %%%%%%%%%%%%%%%%%%%%%
                isokNoNan = ~isnan(argoWorkP.pres_adjusted.data(isokP2,:));
                if any(find(isokNoNan))
                    X = argoWorkP.pres_adjusted.data(isokP2,isokNoNan);
                    Y = argoWorkP.temp_adjusted.data(isokP2,isokNoNan); %marine 11/06/19
                    YY = argoWorkP.psal_adjusted.data(isokP2,isokNoNan); %marine 11/06/19
                    Xq = argoWork.pres_adjusted.data(isokP2,:);
                    
                    
                    
                    if length(find(isokNoNan))>1  %marine 18/06/19
                        %marine 25/06/19 : suppression des doublons pour
                        %l'interpolation
                        [~,X_ind]=unique(X); X_ind=sort(X_ind);
                        X_tmp=X(X_ind); Y_tmp=Y(X_ind); YY_tmp=YY(X_ind);
                        [~,Y_ind]=unique(Y_tmp); 
                        Y_ind=sort(Y_ind);
                        [~,YY_ind]=unique(YY_tmp); YY_ind=sort(YY_ind);

                        if length(Y_ind)>1 %01/07/2019 marine
                            argoWork.temp_adjusted.data(isokP2,:) = interp1(X_tmp(Y_ind),Y_tmp(Y_ind),Xq);
                        else
                            test_isok_Y=abs(Xq-X_tmp(Y_ind))<5;
                            argoWork.temp_adjusted.data(isokP2,test_isok_Y)=Y_tmp(Y_ind);
                        end
                        
                        if length(YY_ind)>1
                            argoWork.psal_adjusted.data(isokP2,:) = interp1(X_tmp(YY_ind),YY_tmp(YY_ind),Xq);
                        else
                            test_isok_YY=abs(Xq-X_tmp(YY_ind))<5;
                            argoWork.psal_adjusted.data(isokP2,test_isok_YY) = YY_tmp(YY_ind);
                        end
                        
                        %Valeur inférieur au minimum du primary D
                        is_inf_secondary=Xq<min(X);
                        ind_diff_max=abs(Xq-min(X))<interp_surf_max;
                        [~,ind_min_P]=min(X);
                        argoWork.temp_adjusted.data(isokP2,and(ind_diff_max,is_inf_secondary))=Y(ind_min_P);
                        argoWork.psal_adjusted.data(isokP2,and(ind_diff_max,is_inf_secondary))=YY(ind_min_P);
                        
                        %Valeur supérieur au maximum du primary D
                        is_sup_secondary=Xq>max(X);
                        ind_diff_max=abs(Xq-max(X))<interp_prof_max;
                        [~,ind_max_P]=max(X);
                        argoWork.temp_adjusted.data(isokP2,and(ind_diff_max,is_sup_secondary))=Y(ind_max_P);
                        argoWork.psal_adjusted.data(isokP2,and(ind_diff_max,is_sup_secondary))=YY(ind_max_P);
                    elseif length(X)==1
                        test_isok=abs(Xq-X)<5;
                        argoWork.temp_adjusted.data(isokP2,test_isok)=Y;
                        argoWork.psal_adjusted.data(isokP2,test_isok)=YY;
                    else %Cas ou il n'y a pas de PSAL ni de TEMP dans les profils descendant primary
                        disp('No PSAL and TEMP data in primary descending profile (--> need to be implemented in LOCODOX in DOXY_argo_prepare function)')
                        %                 %Recherche du profil montant le plus proche en temps
                        %                 is_D='D'==argoP.direction.data;
                        %                 juld_primary=argoP.juld.data;
                        %                 juld_primary(is_D)=NaN;
                        %                 [~,ind_min]=min(abs(argo.juld.data(isokP2)-juld_primary));
                        %                 isokNoNan = ~isnan(argoWorkP.pres_adjusted.data(ind_min,:));
                        %                 X = argoWorkP.pres_adjusted.data(ind_min,isokNoNan);
                        %                 Y = argoWorkP.temp_adjusted.data(ind_min,isokNoNan); %marine 11/06/19
                        %                 YY = argoWorkP.psal_adjusted.data(ind_min,isokNoNan); %marine 11/06/19
                        %                 Xq = argoWork.pres_adjusted.data(isokP2,:);
                        %                 argoWork.temp_adjusted.data(isokP2,:) = interp1(X,Y,Xq);
                        %                 argoWork.psal_adjusted.data(isokP2,:) = interp1(X,YY,Xq);
                    end
                end
                i = i+sum(isokP2);
            else
                i=i+1;
            end
        end
    end
end

% =========================================================================
%% Compute density, Compute depth (m) from pressure
% =========================================================================
%Compute density from psal, temp and pres.
Work.DENSITY = real(sw_pden(argoWork.psal_adjusted.data, ...
              argoWork.temp_adjusted.data, argoWork.pres_adjusted.data,0));
argoWork.an_dens.data = Work.DENSITY - 1000;
argoWork.density.data = Work.DENSITY;

% Control Plot : pressure VS density
if Work.makePlot
    DOXY_PLOT_QC(hfig,1,Work,argo,argoWork);
end

% Compute depth from pres
Work.DEPTH = NaN*ones(size(argoWork.pres_adjusted.data));
for i=1:size(argoWork.pres_adjusted.data,1)
    % compute DEPTH at non NaN values. if all values are NaN, Work.DEPTH
    % remain to NaN.
    isnok = isnan(argoWork.pres_adjusted.data(i,:)) | isnan(argoWork.an_dens.data(i,:));
    if sum(isnok) ~= size(argoWork.pres_adjusted.data,2)
        Work.DEPTH(i,~isnok) = prsenz(argoWork.pres_adjusted.data(i,~isnok),...
            argoWork.an_dens.data(i,~isnok), argo.latitude.data(i));
    end
end


% =========================================================================
%% Quality Flag on Doxy
% =========================================================================

% Remove data with QC flag = 4
[argoWork,Work] = DOXY_QC(wmo,argo,argoWork,Work);

% Control Plot
if Work.makePlot
    DOXY_PLOT_QC(hfig,2,Work,argo,argoWork);
end

%% Compute density
% =========================================================================
% If unpumped psal_adjusted value has only NaN values (so QC = 4), get the
% psal field instead. In delayed mode, all the QC of the unpumped
% psal_ajusted takes QC = 4 and data = NaN, because of the unpumped state.
% But for DOXY correction, we does use the real time psal (and temp, pres),
% with non NaN values du to an artificial flag 4.
if strcmp(Work.whichCorr,'INAIR') && all(all(isnan(argoWork.psal_adjusted.data)))
    Work.DENS = real(sw_pden(argo.psal.data, ...
              argo.temp.data, argo.pres.data,0));    
else
    Work.DENS = real(sw_pden(argoWork.psal_adjusted.data, ...
              argoWork.temp_adjusted.data, argoWork.pres_adjusted.data,0));
end
argoWork.an_dens.data = Work.DENS - 1000;
argoWork.density.data = Work.DENS;


% =========================================================================
%% Compute Saturation and Convert Argo Doxy data in the choosen unit
% =========================================================================
% Saturation and saturation percentage
%   Convert in ml/L to compute the % of saturation (%SAT)
unit = argoWork.doxy_adjusted.units;
tmpDoxy = DOXY_convert(argoWork.doxy_adjusted.data,unit,'mumoL/L',argoWork.an_dens.data);
%   Saturation et saturation percentage
if Work.presEff == 1
    P = argoWork.pres_adjusted.data; 
else
    P = 0; 
end
argoWork.psat.data = O2ctoO2s(tmpDoxy,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data, P); % in mL/L
tmpDoxy = DOXY_convert(argoWork.doxy_adjusted.data,unit,'mL/L',argoWork.an_dens.data);
argoWork.sat.data = (tmpDoxy./argoWork.psat.data)*100;

% Convert Doxy in the choosen unit
argoWork.doxy_adjusted.data = DOXY_convert(argoWork.doxy_adjusted.data,argo.doxy.units,unit,argoWork.an_dens.data);

% =========================================================================
%% Compute PPOX from DOXY
% =========================================================================
tmpDoxy = DOXY_convert(argoWork.doxy_adjusted.data,unit,'mumol/L',argoWork.an_dens.data);
argoWork.ppox_adjusted.data = O2ctoO2p(tmpDoxy,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);

% =========================================================================
%% Finalize Work
% =========================================================================
Work.DOXY_QC = argoWork.doxy_adjusted_qc.data; %Variable des QC
Work.PSAT = argoWork.psat.data;
Work.DOXY = argoWork.doxy_adjusted.data;
Work.PPOX = argoWork.ppox_adjusted.data;
Work.PSAL = argoWork.psal_adjusted.data;
Work.TEMP = argoWork.temp_adjusted.data;
Work.PRES = argoWork.pres_adjusted.data;

% =========================================================================
%% Control Plot
% =========================================================================
if mem == 1
    Work.savePlot = 1;
end

if Work.makePlot
    DOXY_PLOT_QC(hfig,3,Work,argo,argoWork)
end

if exist('mem_plot','var')
    Work.makePlot=mem_plot;
end

end

