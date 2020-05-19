% DOXY_drift corrects data from sensor drift
%
% SYNTAX
% [Work,DRIFT] = DOXY_drift(Work, argoWork, argo, refCyc, whichCorr, argoTrajWork)
%
% DESCRIPTION
% DOXY_drift compare the doxy argo data to the WOA data or NCEP data and estimates a
% sensor drift. If needed, a correction is apply. The drift estimation is
% made at depth below the user choice for WOA data. The WOA data and the argo data are both
% interpolate on a pressure vertical grid every 100m from the depth chosen by the user 
% (see the configuration file) to 6000m.
% For drift computes on NCEP data, in-air measurement are used.
%
% INPUT
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA)
%                           Example:
%                           Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                           where
%                           Work.DOXY_WOA =
%                           name: 'DOXY_ADJUSTED'
%                            dim: {'N_PROF'  'N_LEVELS'}
%                           data: [85x120 double]
%                      long_name: 'DISSOLVED OXYGEN adjusted and converted'
%                  standard_name: 'moles_of_oxygen_per_unit_mass_in_sea_water'
%                     FillValue_: 99999
%                          units: 'mumol/kg'
%                      valid_min: 0
%                      valid_max: 650
%                       C_format: '%9.3f'
%                 FORTRAN_format: 'F9.3'
%                     resolution: 0.0010
%                           type: 5
%
%   argowork (structure)  float working structure issued from the argo data
%                         Example:
%                         argoWork =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                                       sat: [1x1 struct]
%                                      psat: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%
%   argo (structure)      Argo float structure (get directly from data).
%                         Example :
%                         argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%
%                               argo.pres
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
%  whichCorr (string)      Indicates which kind of correction is used
%                          Example: whichCorr='REF', whichCorr='WOA',
%                          whichCorr = 'INAIR';
%
%  argoTrajWork (struct)   Float trajectory data structure
%
% OUTPUT
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA),
%                        completed with the annual drift.
%                        Example:
%                        Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                                     DRIFT: 4.8400
%
%   DRIFT (structure)    Drift estimating structure, from argo data and
%                        climatology data.
%                        Example:
%                        DRIFT =
%                                      doxy: [46x85 double]
%                                   doxyref: [46x85 double]
%                                   deltaO2: [46x85 double]
%                              deltaO2_Mean: [1x85 double]
%                            doxy_corrected: [85x120 double]
%                          fitRegression: [1x85 double]
%
%
%
% CALL :
%   DOXY_PLOT_drift, O2ptoO2c, O2ptoO2s, O2stoO2c
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 07/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%             25/03/2019        Marine GALLIAN, Altran Ouest
%                               Fill the summary file
%                               Manage onvertion for drift computed on NCEP
%                               and WOA
%                               Improve message box
%             30/03/2019        Marine GALLIAN, Altran Ouest
%                               Now drift is computed using PRES and no
%                               more using the DEPTH variable
%       v3.4 10.04.2020         T.Reynaud temporal ==> time


function [Work, DRIFT] = DOXY_drift(Work, argoWork, argo, argoTrajWork, whichCorr)
% =========================================================================
%% Initialisation
% =========================================================================
if ~isfield(Work,'refCyc')
    okdim = strcmp(argoWork.pres_adjusted.dim,'N_PROF');
    dim = size(argoWork.pres_adjusted.data);
    refCyc = 1:dim(okdim);
else
    refCyc = Work.refCyc;
end

fprintf('\t Drift correction is done with PRES\n');


%nbr_lin_reg= Minimum number of points to compute a drift. Default value :10.
nbr_lin_reg=10;

% for saving plot
if Work.savePlot == 1
    Work.savePlot = 0;
    mem = 1;
else
    mem = 0;
end

% depth at which to compute the drift and unit
if  ~strcmp(whichCorr,'INAIR')
    maxdepth=abs(nanmin(-Work.PRES(:)));
    unit = argoWork.doxy_adjusted.units;
    Zq=Work.min_drift_depth:100:6000; %28/06/19
else %drift compute with NCEP data
    Zq = 1;
    unit = argoTrajWork.ppox_doxy_adjusted.units;
end

% WOA and REF correction : find the DOXY reference interpolated on argo :
% DOXY_REF_interpP or DOXY_WOA_interP
doxyRefField = sprintf('DOXY_%s_interpP',whichCorr);

% initialize
tmpTab = NaN(length(refCyc),length(Zq));
DRIFT.doxy = tmpTab;
DRIFT.doxyref = tmpTab;
DRIFT.deltaO2 = tmpTab;
DRIFT.deltaO2_Mean = NaN(length(refCyc),1);

% =========================================================================
%% Prepare data for drift computation
% =========================================================================
% -------------------------------------------------------------------------
% For drift compute on WOA
% Interpolate data at deep depth, every 100m
% If no data under min(Zq), no interpolation (NaN)
% -------------------------------------------------------------------------
if strcmp(whichCorr,'WOA') || strcmp(whichCorr,'REF') %strcmp(whichDrift,'WOA')
    
    if length(refCyc)<nbr_lin_reg
        fprintf('\t No Drift Computation: not enougth data for this float');
        waitfor(warndlg(sprintf('\t No Drift Computation: not enougth data for this float (%d cycle) \n', length(refCyc)) ,...
            'DOXY_drift.m')); %marine 11/06/19
        Work.DODRIFT = false;
        if Work.makePlot
            % Nb of day after deployment
            daydiff = (argo.juld.data - argo.juld.data(1));
            DOXY_PLOT_drift(DRIFT, Work, daydiff,unit,1);
            if mem == 1  %marine 11/06/19
                Work.savePlot = 1;
            end
            return
        end
        
    else
        for i = refCyc
            isok = find(~isnan(Work.PRES(i,:)) & Work.PRES(i,:) >= Work.min_drift_depth);
            Zar = Work.PRES(i,isok);
            
            % At least two points to allow interpolate
            if  length(Zar) < 2
                DRIFT.doxy(i,:) = NaN;
                DRIFT.doxyref(i,:) = NaN;
            else
                Var = Work.DOXY(i,isok);
                Vref = Work.(doxyRefField)(i,isok);
                DRIFT.doxy(i,:) = interp1(Zar,Var,Zq,'linear');
                DRIFT.doxyref(i,:) = interp1(Zar,Vref,Zq,'linear');
            end
            DRIFT.deltaO2(i,:) = DRIFT.doxy(i,:) - DRIFT.doxyref(i,:);
            DRIFT.deltaO2_Mean(i,:) = nanmean(DRIFT.deltaO2(i,:));
        end
        
        if all(all(isnan(DRIFT.doxy)))
            fprintf('\t No Drift Computation: no data below %d m\n', min(Zq));
            waitfor(warndlg(sprintf('\t No Drift Computation: no data below %d m\n --> The minimum depth for drift calculation can be changed in the configuration file\n', min(Zq)),...
                'DOXY_drift.m')); %marine 11/06/19
            Work.DODRIFT = false;
            if Work.makePlot
                % Nb of day after deployment
                daydiff = (argo.juld.data - argo.juld.data(1));
                DOXY_PLOT_drift(DRIFT, Work, daydiff,unit,1);
                if mem == 1  %marine 11/06/19
                    Work.savePlot = 1;
                end
                return
            end
        end
    end
    
    % Nb of day after deployment
    daydiff = (argo.juld.data - argo.juld.data(1));
end

% -------------------------------------------------------------------------
% For INAIR correction
% -------------------------------------------------------------------------

if strcmp(whichCorr,'INAIR') %strcmp(whichDrift,'NCEP')
    DRIFT.doxy = Work.surface.argo.pO2_air;
    DRIFT.doxyref = Work.surface.ncep.pO2;
    
    DRIFT.deltaO2 = (DRIFT.doxy - DRIFT.doxyref)';
    DRIFT.deltaO2_Mean = nanmean(DRIFT.deltaO2);
    
    % Nb of day after deployment
    %   daydiff = (argoTrajWork.juld.data - min(argoTrajWork.juld.data))';
    %--EB : daymin was set to argoTrajWork.juld.data{1}(1), whereas it
    % could be null or empty
    if Work.NSIAfloat
        daymin = Work.surface.argo.juld_air(1);
        daydiff = Work.surface.argo.juld_air - daymin;
        nbr_meas=Work.nbrMeas.inAir;
    else
        daymin = argoTrajWork.juld.data(1);
        daydiff = argoTrajWork.juld.data - daymin;
    end
end

%--EB - 20180418
% Uniformize the data format (no single and double in the same field)
if iscell(DRIFT.doxy)
    driftFields = fieldnames(DRIFT);
    nbfield = numel(driftFields);
    for f = 1:nbfield
        tmpfield = driftFields{f};
        DRIFT.(tmpfield) = cellfun(@double,DRIFT.(tmpfield),'UniformOutput',false);
    end
end
%--EB - 20180418

% -------------------------------------------------------------------------
%% Control Plot
% Plot interpolated doxy from argo and doxy from WOA respect to the number
% of day after deployement, and plot the difference.
% Or plot the ppox from argo inflated and from NCEP respect to the number
% of day after deployement, and plot the difference.
% -------------------------------------------------------------------------
if Work.makePlot
    [ax,~] = DOXY_PLOT_drift(DRIFT, Work, daydiff,unit,1);
end

if iscell(daydiff)
    daydiff = cell2mat(daydiff);
end

% =========================================================================
%% Drift computation
%  Compute the linear regression off the doxy difference timeserie
%  (argo - woa, or argo - NCEP).
%  Based on Takeshita
% =========================================================================

if Work.drift_spec == 1
    v =  ver;
    if all(cellfun('isempty',strfind({v.Name},'Statistics')))
        input = inputdlg({'Statistics Toolbox unavailable : Matlab can''t apply your drift equation. The drift will be compute with a polynomial equation. The default degree is 1. Do you wan''t to change it ?'},...
            'DOXY_drift : non statistic toolbox',1,{'1'});
        Work.drift_spec = 0;
        Work.drift_fitPolynomialDegree = str2num(char(input)); %#ok<ST2NM>
    end
end
if Work.drift_spec == 0
    coeffFit = NaN(Work.drift_fitPolynomialDegree+1,length(Zq));
else
    coeffFit = NaN(numcoeffs(Work.drift_fittype),length(Zq));
end
cpt = 0;
noDriftComputation = false;

for z=1:length(Zq)
    % linear regression
    if strcmp(whichCorr,'WOA') || strcmp(whichCorr,'REF')
        
        isok = ~isnan(daydiff) & ~isnan(DRIFT.deltaO2(:,z));
        if sum(isok) >= nbr_lin_reg && any(~isnan(DRIFT.deltaO2(:,z)))
            if Work.drift_spec == 0
                fitResult = polyfit(double(daydiff(isok)),double(DRIFT.deltaO2(isok,z)),Work.drift_fitPolynomialDegree);
                coeffFit(:,z) = fitResult;
            else
                fitResult = fit(double(daydiff(isok)),double(DRIFT.deltaO2(isok,z)),Work.drift_fittype);
                coeffFit(:,z) = coeffvalues(fitResult);
            end
            
            % control plot : check the coefficient and the linear regression at
            % a local depth
            if Work.makePlot
                if Work.drift_spec == 0
                    DRIFT.poliv_c = polyval(coeffFit(:,z), double(daydiff));
                else
                    DRIFT.poliv_c = fitResult(double(daydiff));
                end
                [ax,~] = DOXY_PLOT_drift(DRIFT,Work,daydiff,unit,2,ax);
                DRIFT = rmfield(DRIFT,'poliv_c');
            end
            
        else
            cpt = cpt+1;
            if z==length(Zq) && cpt==length(Zq)
                errorMess = strcat('Less than [', num2str(nbr_lin_reg) ,'] valid data below 1500m for each profile');
                noDriftComputation = true;
            end
        end
        
    elseif strcmp(whichCorr,'INAIR')
        data = double(DRIFT.deltaO2)';
        isok = ~isnan(data) & ~isnan(daydiff);
        ind_conv=0;
        
        if isfield(Work,'fitResult_WOA')
            coeffFit(:,z) = Work.fitResult_WOA;
            ind_conv=1; %PPOX has to be converted to apply the drift
        elseif sum(isok) >= nbr_lin_reg %|| any(~cellfun(@isnan,DRIFT.deltaO2,'UniformOutput',false))
            if ~Work.drift_spec
                fitResult = polyfit(daydiff(isok),data(isok),Work.drift_fitPolynomialDegree);
                coeffFit(:,z) = fitResult;
            else
                % fit doesn't accept row vectors.
                if ~iscolumn(daydiff)
                    daydiff = daydiff';
                end
                if ~iscolumn(data)
                    data = data';
                end
                fitResult = fit(daydiff(isok),data(isok),Work.drift_fittype);
                coeffFit(:,z) = coeffvalues(fitResult);
            end
            % control plot : check the coefficient and the linear regression at
            % a local depth
            if Work.makePlot
                if Work.drift_spec == 0
                    DRIFT.poliv_c = polyval(coeffFit(:,z), double(daydiff));
                else
                    DRIFT.poliv_c = fitResult(double(daydiff));
                end
                [ax,~] = DOXY_PLOT_drift(DRIFT,Work,daydiff,unit,2,ax);
                DRIFT = rmfield(DRIFT,'poliv_c');
            end
        else
            cpt = cpt+1;
            if z==length(Zq) && cpt==length(Zq)
                errorMess = 'Less than 10 valid data at surface';
                noDriftComputation = true;
            end
        end
    end
end

% =========================================================================
%Compute drift only if it is possible
% =========================================================================
if ~noDriftComputation
    % -------------------------------------------------------------------------
    % Slope averaging
    % -------------------------------------------------------------------------
    coeffFit_mean = nanmean(coeffFit,2);
    
    % For polynomial regression for drift calculation, no offset is taking into
    % account. for example, if the regression is like ax+b, the drift is only
    % the slope : a.
    if Work.drift_spec == 0
        coeffFit_mean_slope = coeffFit_mean(1:end-1);
        % apply the slope, not the offset.
        coeffFit_mean_slope = [coeffFit_mean_slope; 0];
    else
        coeffFit_mean_slope = coeffFit_mean;
    end
    
    isok = ~isnan(daydiff);
    % control plot : check the mean coefficient and the linear regression for
    % all depth
    if Work.makePlot
        if mem == 1
            Work.savePlot = 1;
        end
        if Work.drift_spec == 0
            DRIFT.fitRegression = polyval(coeffFit_mean,double(daydiff));
        else
            listCoef = [];
            for i=1:length(coeffFit_mean)
                listCoef = [listCoef,sprintf('coeffFit_mean(%d),',i)];
            end
            listCoef(end) = '';
            eval(['c = cfit(Work.drift_fittype,' listCoef ');']);
            DRIFT.fitRegression = feval(c,double(daydiff));
        end
        [~,Work]=DOXY_PLOT_drift(DRIFT,Work,daydiff(isok),unit,3, ax);
    end
    
    %Recompute fitRegression with offset = 0 now
    if strcmp(whichCorr,'INAIR') && ind_conv~=1 %Case INAIR correction and drift computed on NCEP
        %Change of time we use profile time and no more in air or in water
        daydiff = (argo.juld.data - argo.juld.data(1));
        
        if Work.drift_spec == 0
            DRIFT.fitRegression = polyval(coeffFit_mean_slope,double(daydiff));
        else
            listCoef = [];
            for i=1:length(coeffFit_mean)
                listCoef = [listCoef,sprintf('coeffFit_mean_slope(%d),',i)];
            end
            listCoef(end) = '';
            eval(['c = cfit(Work.drift_fittype,' listCoef ');']);
            DRIFT.fitRegression = feval(c,double(daydiff));
        end
        if size(DRIFT.fitRegression,2) == size(Work.DOXY,1)
            DRIFT.fitRegression = DRIFT.fitRegression';
        end
    elseif Work.drift_spec == 0
        DRIFT.fitRegression = polyval(coeffFit_mean_slope,double(daydiff));
    else
        listCoef = [];
        for i=1:length(coeffFit_mean_slope)
            listCoef = [listCoef,sprintf('coeffFit_mean_slope(%d),',i)];
        end
        listCoef(end) = '';
        eval(['c = cfit(Work.drift_fittype,' listCoef ');']);
        DRIFT.fitRegression = feval(c,double(daydiff));
    end
    
    % Compute drift for in air and in water data
    % Change time reference for in air and in water data
    if strcmp(whichCorr,'INAIR') && ~isfield(Work,'isnotINAIRcorr')
        
        if Work.NSIAfloat
            daymin_airwater = min(Work.surface.argo.juld_air);
            daydiff_airwater = Work.surface.argo.juld_air - daymin_airwater;
        else
            daymin_airwater = min(min(argoTrajWork.juld.data));
            daydiff_airwater = argoTrajWork.juld.data - daymin_airwater;
        end
        %Recompute drift
        if Work.drift_spec == 0
            DRIFT.fitRegression_airwater = polyval(coeffFit_mean_slope,double(daydiff_airwater));
        else
            listCoef_airwater = [];
            for i=1:length(coeffFit_mean)
                listCoef_airwater = [listCoef_airwater,sprintf('coeffFit_mean_slope(%d),',i)];
            end
            listCoef_airwater(end) = '';
            eval(['c = cfit(Work.drift_fittype,' listCoef_airwater ');']);
            DRIFT.fitRegression_airwater = feval(c,double(daydiff_airwater));
        end
        if size(DRIFT.fitRegression,2) == size(Work.DOXY,1)
            DRIFT.fitRegression_airwater = DRIFT.fitRegression_airwater';
        end
        
        %Apply constant drift
        if isfield(Work,'ind_drift_stop')
            for i=1:length(daydiff_airwater)
                if daydiff_airwater(i)>Work.ind_drift_stop(2) && Work.ind_drift_stop(2)>-1
                    Work.ind_drift_stop_airwater=[i,Work.ind_drift_stop(2)];
                    break
                end
            end
            DRIFT.fitRegression_airwater(Work.ind_drift_stop_airwater(1):end)=DRIFT.fitRegression_airwater(Work.ind_drift_stop_airwater(1));
        end
        
        
    end
    
    % -------------------------------------------------------------------------
    % Annual drift
    % -------------------------------------------------------------------------
    if ~strcmp(whichCorr,'INAIR')
        Work.DRIFT = (DRIFT.fitRegression(end) - DRIFT.fitRegression(1)) / ...
            ((nanmax(argo.juld.data)-nanmin(argo.juld.data))/365);
    else
        if iscell(argoTrajWork.juld.data)
            Work.DRIFT = (DRIFT.fitRegression(end) - DRIFT.fitRegression(1)) / ...
                ((nanmax(cell2mat(argoTrajWork.juld.data))-nanmin(cell2mat(argoTrajWork.juld.data)))/365);
        else
            Work.DRIFT = (DRIFT.fitRegression(end) - DRIFT.fitRegression(1)) / ...
                ((nanmax(argoTrajWork.juld.data)-nanmin(argoTrajWork.juld.data))/365);
        end
    end
    
    
    % =========================================================================
    %% Asking user for drift
    % =========================================================================
    if isfield(Work,'ind_drift_stop')
        explanation3=sprintf('A constant drift is applied from day number %d',Work.ind_drift_stop(2));
    else
        explanation3=' ';
    end
    
    aux_year=find((abs(daydiff-365.25))==min(abs(daydiff-365.25))); %% position after a year life
    nyears_floatlife = (argo.juld.data(end)-argo.juld.data(1))/365;
    
    %Display leading coefficient %marine
    if strcmp(whichCorr,'INAIR')
        if Work.drift_fitPolynomialDegree == 1
            explanation0 = sprintf('The best fitting was: PO2_corr(t) ~= PO2(t) + %2.4f * t \n',coeffFit_mean_slope(1));
        elseif Work.drift_fitPolynomialDegree == 2
            explanation0 = sprintf('The best fitting was: PO2_corr(t) ~= PO2(t) + %2.4f * t + %2.4f * t^2\n',coeffFit_mean_slope(1),coeffFit_mean_slope(2));
        else
            explanation0 = '';
        end
    else
        if Work.drift_fitPolynomialDegree == 1
            explanation0 = sprintf('The best fitting was : DOXY_corr(t) = DOXY(t) + %2.4f * t \n',coeffFit_mean_slope(1));
        elseif Work.drift_fitPolynomialDegree == 2
            explanation0 = sprintf('The best fitting was : DOXY_corr(t) = DOXY(t) + %2.4f * t + %2.4f * t^2\n',coeffFit_mean_slope(1),coeffFit_mean_slope(2));
        elseif Work.drift_fitPolynomialDegree == 3
            explanation0 = sprintf('The best fitting was : DOXY_corr(t) = DOXY(t) + %2.4f * t + %2.4f * t^2 + %2.4f * t^3\n',coeffFit_mean_slope(1),coeffFit_mean_slope(2),coeffFit_mean_slope(3)); '';
        end
    end
    
    %Display drift after a year and at the end of life
    if nyears_floatlife>=1
        explanation1 = sprintf('The drift after a year is %4.3f %s \n',DRIFT.fitRegression(aux_year), unit);
        explanation2 = sprintf('The drift at the end of the float life ( %4.2f years) is %4.3f %s \n', nyears_floatlife, DRIFT.fitRegression(end),unit);
    else
        explanation1 = sprintf('The drift at the end of the float life ( %4.2f years) is %4.3f %s \n', nyears_floatlife, DRIFT.fitRegression(end),unit);
        explanation2 = sprintf(' \n');
    end
    
    %Message of sugestion depending of the availability of data
    deltaO2MeanNok = sum(isnan(DRIFT.deltaO2_Mean));
    if daydiff(end) < 365
        % Less than one year of data
        driftCorrSuggestion = sprintf('* Suggestion: NO DRIFT CORRECTION:');
        explanation = sprintf('* State: less than 1 year of data');
    elseif length(argo.cycle_number.data) - deltaO2MeanNok<nbr_lin_reg+1
        % less than 5 profils with depth > 1500m
        driftCorrSuggestion = sprintf('* Suggestion: NO DRIFT CORRECTION:');
        explanation =sprintf( '* State: less than 10 profiles > 1500m');
    else
        explanation = ' ';
        driftCorrSuggestion= ' ';
    end
    
    %Choice of user
    if strcmp(whichCorr,'INAIR') && strcmp(Work.whichDrift,'WOA')
        if isfield(Work,'fitResult_WOA')
            applyDriftC='YES';
        else
            applyDriftC='NO';
        end
    else
        applyDriftC = questdlg({sprintf('Float #%d',Work.wmo);...
            ' ---------------------------';
            'The statistical method results in: '; ...
            ' ---------------------------';explanation0;
            explanation1; explanation2; ...
            '---------------------------';...
            driftCorrSuggestion; explanation; explanation3;...
            ' ---------------------------';...
            'Apply time drift correction ?';''},...
            'TIME DRIFT CORRECTION','YES','NO','YES');
    end
    
else
    if mem==1
        Work.savePlot=1;
    end
    driftCorrSuggestion = '* DRIFT CORRECTION NOT POSSIBLE : no data';
    explanation = sprintf('* State: %s',errorMess);
    applyDriftC='NO';
    % to small drift
    waitfor(warndlg({driftCorrSuggestion;explanation},'NO DRIFT COMPUTATION')); %marine 12/06/19
end

% =========================================================================
%% apply drift to data
% =========================================================================
switch applyDriftC
    case 'YES'
        fprintf('\t The drift is applied \n')       
        % =========================================================================
        %% Compute drift correction over profile and trajectory
        % =========================================================================
        Work.DODRIFT = true;
        % -------------------------------------------------------------------------
        % Compute the profile correction
        % -------------------------------------------------------------------------
        
        %Variable Initialisation
        okdim = strcmp(argoWork.pres_adjusted.dim,'N_PROF');
        dim = size(argoWork.pres_adjusted.data);
        DRIFT.doxy_corrected = NaN(dim);
        if strcmp(whichCorr,'INAIR')
            DRIFT.ppox_corrected = NaN(size(argoWork.pres_adjusted.data));
        end
        
        %     if size(DRIFT.fitRegression,2) == size(Work.DOXY,1)
        %         DRIFT.fitRegression = DRIFT.fitRegression';
        %     end
    
        if isfield(Work,'isnotREFcorr') %Case NCEP and drift on WOA
            Work.DOXY_DRIFTCORR_COEF = coeffFit_mean_slope;
            
        elseif isfield(Work,'isnotINAIRcorr') % case WOA or REF and drift on NCEP %marine 20/06/19
            %convert DOXY to PPOX
            if Work.presEff == 1, P = Work.pres_convert; else, P = 0; end
            tmpDoxy = DOXY_convert(Work.DOXY_convert,Work.units_convert,'mumol/L',Work.dens_convert);
            PPOX_convert=O2ctoO2p(tmpDoxy,Work.temp_convert, Work.psal_convert,P);
            daydiff = (argo.juld.data - argo.juld.data(1));
            DRIFT.fitRegression_convert = polyval(coeffFit_mean_slope,double(daydiff));
            %Compute constant drift
            if isfield(Work,'ind_drift_stop')
                for i=1:length(daydiff)
                    if daydiff(i)>Work.ind_drift_stop(2) && Work.ind_drift_stop(2)>-1
                        Work.ind_drift_stop=[i,Work.ind_drift_stop(2)];
                        break
                    end
                end
                DRIFT.fitRegression_convert(Work.ind_drift_stop(1):end)=DRIFT.fitRegression_convert(Work.ind_drift_stop(1));
            end           
            %Apply drift
            for i = 1:dim(okdim)
                DRIFT.ppox_corrected(i,:) = PPOX_convert(i,:) - DRIFT.fitRegression_convert(i);
            end
            %reconvert PPOX to DOXY and PSAT
            DRIFT.doxy_corrected=O2ptoO2c(DRIFT.ppox_corrected,Work.temp_convert,Work.psal_convert,P);
            DRIFT.psat_corrected=O2ptoO2s(DRIFT.ppox_corrected,Work.temp_convert,Work.psal_convert,P);
            %Modifie output of the function
            Work.PPOX_DRIFTCORR_COEF = coeffFit_mean_slope;
            drift_on='NCEP';
            
        elseif ~strcmp(whichCorr,'INAIR') %case WOA or REF with drift on WOA
            if Work.presEff == 1, P = Work.pres_convert; else, P = 0; end
            %Compute constant drift
            if isfield(Work,'ind_drift_stop')
                DRIFT.fitRegression(Work.ind_drift_stop(1):end)=DRIFT.fitRegression(Work.ind_drift_stop(1));
            end
            %Apply drift
            for i = 1:dim(okdim)
                DRIFT.doxy_corrected(i,:) = Work.DOXY(i,:) - DRIFT.fitRegression(i);
            end
            tmpDoxy = DOXY_convert(DRIFT.doxy_corrected,argoWork.doxy_adjusted.units,'mumol/L',argoWork.an_dens.data);
            DRIFT.psat_corrected = O2ctoO2s(tmpDoxy,argoWork.temp_adjusted.data, argoWork.psal_adjusted.data,P);
            Work.DOXY_DRIFTCORR_COEF = coeffFit_mean_slope;
            drift_on='WOA';
        elseif (strcmp(whichCorr,'INAIR') && ind_conv==1) %case INAIR with drift on WOA
            %Compute constant drift            
            if isfield(Work,'ind_drift_stop')
                DRIFT.fitRegression(Work.ind_drift_stop(1):end)=DRIFT.fitRegression(Work.ind_drift_stop(1));
            end         
            %Apply drift            
            for i = 1:dim(okdim)
                DRIFT.doxy_corrected(i,:) = Work.DOXY(i,:) - DRIFT.fitRegression(i);
            end
            if Work.presEff == 1, P = argoWork.pres_adjusted.data; else, P = 0; end
            %convert all data
            tmpDoxy = DOXY_convert(DRIFT.doxy_corrected,argoWork.doxy_adjusted.units,'mumol/L',argoWork.an_dens.data);
            DRIFT.ppox_corrected =O2ctoO2p(tmpDoxy,argoWork.temp_adjusted.data, argoWork.psal_adjusted.data,P);
            DRIFT.psat_corrected = O2ctoO2s(tmpDoxy,argoWork.temp_adjusted.data, argoWork.psal_adjusted.data,P);
            
            %Apply drift to in air and in water data (if there is a
            %constant drift, it is already computed in DRIFT.fitRegression_airwater
            tmp_doxy_inair=O2ptoO2c(Work.surface.argo.pO2_air,Work.surface.argo.temp_air, 0,0);
            tmp_doxy_inwater=O2ptoO2c(Work.surface.argo.pO2_water,Work.surface.argo.temp_water, Work.surface.argo.sss,0);
            tmp_doxy_inwater_corr = tmp_doxy_inwater - DRIFT.fitRegression_airwater;
            tmp_doxy_inair_corr = tmp_doxy_inair - DRIFT.fitRegression_airwater;
            DRIFT.pO2_corrected.inair = O2ctoO2p(tmp_doxy_inair_corr,Work.surface.argo.temp_air, 0,0);
            DRIFT.pO2_corrected.inwater = O2ctoO2p(tmp_doxy_inwater_corr,Work.surface.argo.temp_water, Work.surface.argo.sss,0);
            
            %Modifie output of the function
            Work.PPOX_DRIFTCORR = DRIFT.ppox_corrected;
            Work.PPOX_DRIFTCORR_inair = DRIFT.pO2_corrected.inair;
            Work.PPOX_DRIFTCORR_inwater = DRIFT.pO2_corrected.inwater;
            Work.DOXY_DRIFTCORR_COEF = coeffFit_mean_slope;
            drift_on='WOA';
            
        elseif strcmp(whichCorr,'INAIR') && ind_conv~=1 %Case INAIR correction and drift computed on NCEP
            
            if isfield(Work,'ind_drift_stop')
                DRIFT.fitRegression(Work.ind_drift_stop(1):end)=DRIFT.fitRegression(Work.ind_drift_stop(1));
            end   
            
            %Apply drift to ppox data and convert to psat and doxy (near
            %surface sampling)
            for i = 1:dim(okdim)
                DRIFT.ppox_corrected(i,:) = Work.PPOX(i,:) - DRIFT.fitRegression(i);
            end
            if Work.presEff == 1, P = argoWork.pres_adjusted.data; else, P = 0; end
            
            DRIFT.doxy_corrected = O2ptoO2c(DRIFT.ppox_corrected,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);
            DRIFT.psat_corrected = O2ptoO2s(DRIFT.ppox_corrected,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);
            
            %Apply drift to in air and in water data
            %Change time reference for in air and in water data
            DRIFT.pO2_corrected.inair = Work.surface.argo.pO2_air - DRIFT.fitRegression_airwater;
            DRIFT.pO2_corrected.inwater = Work.surface.argo.pO2_water - DRIFT.fitRegression_airwater;
            
            %Modifie output of the function
            Work.PPOX_DRIFTCORR = DRIFT.ppox_corrected;
            Work.PPOX_DRIFTCORR_inair = DRIFT.pO2_corrected.inair;
            Work.PPOX_DRIFTCORR_inwater = DRIFT.pO2_corrected.inwater;
            Work.PPOX_DRIFTCORR_COEF = coeffFit_mean_slope;
            drift_on='NCEP';
        end
        
        
        if ~isfield(Work,'isnotREFcorr')
            %Modifie output of the function
            DRIFT.psat_corrected(DRIFT.psat_corrected<=0) = NaN;
            DRIFT.doxy_corrected(DRIFT.doxy_corrected<=0) = NaN;
            Work.DOXY_DRIFTCORR = DRIFT.doxy_corrected;
            Work.PSAT_DRIFTCORR = DRIFT.psat_corrected;
            
            % -------------- SUMMARY FILE -----------------
            if ~(isfield(Work,'fitResult_WOA') && strcmp(whichCorr,'INAIR')) %condition to avoid writing in double the drift info
                fprintf(Work.log_summ,'Time drift correction applied : \n');
                if Work.drift_spec == 1
                    fprintf(Work.log_summ,sprintf('Time drift computed on %s and estimated by applying a fitting proposed by the user. The fitting equation is %s\n',drift_on,formula(Work.drift_fittype)));
                    fprintf(Work.log_summ,explanation1);
                    fprintf(Work.log_summ,explanation2);
                elseif Work.drift_fitPolynomialDegree == 1
                    fprintf(Work.log_summ,sprintf('Time drift computed on %s and estimated by applying a linear regression.\n',drift_on));
                    fprintf(Work.log_summ,explanation0);
                    fprintf(Work.log_summ,explanation1);
                    fprintf(Work.log_summ,explanation2);
                elseif Work.drift_fitPolynomialDegree == 2
                    fprintf(Work.log_summ, sprintf('Time drift computed on %s and estimated by applying a polynomial regression of order two.\n',drift_on));
                    fprintf(Work.log_summ,explanation0);
                    fprintf(Work.log_summ,explanation1);
                    fprintf(Work.log_summ,explanation2);
                elseif Work.drift_fitPolynomialDegree == 3
                    fprintf(Work.log_summ, sprintf('Time drift computed on %s and estimated by applying a polynomial regression of order three.\n',drift_on));
                    fprintf(Work.log_summ,explanation0);
                    fprintf(Work.log_summ,explanation1);
                    fprintf(Work.log_summ,explanation2);
                end
                
                if isfield(Work,'ind_drift_stop')
                    fprintf(Work.log_summ,sprintf('A constant drift was applied from day number %d\n',Work.ind_drift_stop(2)));
                end
                fprintf(Work.log_summ,'\n');
                fprintf(Work.log_summ,'-----------------------------------------------------------\n');
            end
        end
        
        
        
    case 'NO'
        fprintf('\t Time drift is not applied \n') 
        Work.DODRIFT = false;
        DRIFT.psat_corrected = Work.PSAT;
        DRIFT.doxy_corrected = Work.DOXY;
        DRIFT.ppox_corrected = Work.PPOX;
        
        % -------------- SUMMARY FILE -----------------
        fprintf(Work.log_summ,'Time drift correction not applied. \n');
        fprintf(Work.log_summ,'\n');
        fprintf(Work.log_summ,'-----------------------------------------------------------\n');
        
end
