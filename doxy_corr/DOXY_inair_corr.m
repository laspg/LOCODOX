% DOXY_inair_corr do the correction of PPOX argo data using NCEP data
%
% SYNTAX
% [argo1Struct, argo2Struct, argo3Struct, argo4Struct, WOA,NCEP, Work,goProg] 
%   = DOXY_inair_corr(CONFIG, WOA, NCEP, argo1Struct, argo2Struct, ...
%   argo3Struct, argo4Struct, argo, argoWork, Work, argoTrajWork)
%
% DESCRIPTION
% DOXY_inair_corr do the DOXY argo data correction, through different steps :
% * Colocalize NCEP data over the argo time, lat, lon for InAir data
% * Interpolate WOA data over the argo positions and its vertical grid
% * Compute data drift from WOA data or NCEP data
% * Compute non linear regression between inflated data (PO2 argo in air),
%   deflated data (PO2 argo in water subsurface) and air data (PO2 from NCEP)
%   Bittig et al., 2015 :  Po2inflated = c*PO2deflated + [(1-c)/m]PO2air
% * Apply correction to both ascending and descending profile
% * Apply quality flags to DOXY
%
% INPUT
%   CONFIG (struct)       Configuration Structure
%
%   WOA (struct)          Climatology WOA structure.
%                         Example:
%                          WOA =
%                                init: [1x1 struct]
%                              interp: [1x1 struct]
%                          WOA.init =
%                                dimorder: 'C'
%                                latitude: [1x1 struct]
%                               longitude: [1x1 struct]
%                                   depth: [1x1 struct]
%                                    time: [1x1 struct]fitnlm
%                                 doxywoa: [1x1 struct]
%                                           ...
%
%   NCEP (struct)         NCEP atmospheric structure.
%                         Example:
%                          NCEP =
%                               init: [1x1 struct]
%                              coloc: [1x1 struct]
%                          NCEP.init
%                             nodata: 0
%                                slp: [1×1 struct]
%                                air: [1×1 struct]
%                               rhum: [1×1 struct]
%                                           ...
%   argo1Struct (struct)    Three data structures organised as the main profile
%   argo2Struct (struct)    for working is in argo1Struct (vertical
%   argo3Struct (struct)    sampling shceme of interest for the kind of
%   argo4Struct (struct)    correction), and the other in the second ones.
%                           The program DOXY_argo_prepare_main.m describes
%                           the choice of the main and other profiles.
%                           Example:
%                             argo1Struct = fitnlm
%                                     argo: [1x1 struct]
%                                      Dim: [1x1 struct]
%                                     Work: [1x1 struct]
%                                 argoWork: [1x1 struct]
%                                      VSS: 'Near-surface sampling'
%
%   argo (structure)      Argo float structure (get directly from data),
%                         for ascending profile and Primary profile.
%                           Example :
%                               argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%
%    Work (structure)     Doxy correction working structure, issued and
%                         computed from argo float data.
%                         Example:
%                             Work =
%                                   readme: [1x1 struct]
%                                     unit: 'mumol/kg'
%                                      wmo: 1901205
%                                 DOXY_RAW: [85x120 single]
%                                  DOXY_QC: [85x120 char]
%
%    argowork (structure)  Float working structure issued from the argo data
%                         argoWork = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   ..........
%
%    argoTrajWork (struct) float data structure (argoTrajWork) choose for the
%                        Doxy correction, read in the trajectory files.
%                        Example: argoTrajWork
%                          float_serial_no: [1x1 struct]
%                            wmo_inst_type: [1x1 struct]
%                                     juld: [1x1 struct]
%                                 latitude: [1x1 struct]
%                                longitude: [1x1 struct]
%                                   ..........
%
% OUTPUT
%    argo1Struct (structure)  argo1Struct, argo2Struct and argo3Struct completed
%    argo2Struct
%    argo3Struct
%    argo4Struct
%
%    WOA (struct)          WOA.interp completed
%
%    NCEP (struct)         NCEP.coloc completed
%
%    Work (struct)         Work completed, in particular with coefficient
%                          of DOXY linear correction.
%
%    goProg (double)       0\1. If this index is 0, LOCODOX stop the
%                          treatment in progress
%
% CALL : 
% DOXY_NCEP_colocalize, DOXY_interp_WOA_main, DOXY_interp_REF_main, 
% DOXY_get_primary_PTS_for_traj, DOXY_get_profile_temp_for_traj, 
% DOXY_drift, DOXY_corr_compute, DOXY_corr_apply_main,
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 03/05/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%              v2 09/11/2017   Emilie Brion, Altran Ouest
%                              argument WOA and NCEP (input and output) was
%                              melt and/or forgotten
%               22/03/2019     Marine GALLIAN, Altran Ouest 
%                              Add possibility to compute drift on WOA or
%                              on NCEP 
%                              Compute reference profile if it exists
%              12/03/2020      Thierry Reynaud
%                              whichDrift was missing for argo4Struct

function [argo1Struct,argo2Struct,argo3Struct,argo4Struct,WOA,NCEP, Work,goProg] = DOXY_inair_corr(CONFIG, WOA,...
    NCEP, argo1Struct, argo2Struct, argo3Struct, argo4Struct, argo, argoWork, Work, argoTrajWork)


% =========================================================================
%% Colocalize NCEP data over the argo time, lat, lon for InAir data
% =========================================================================
[NCEP] = DOXY_NCEP_colocalize(NCEP,argoTrajWork,CONFIG);

% =========================================================================
%% Interpolate WOA data over the argo positions and its vertical
% grid
% =========================================================================
% -------------------------------------------------------------------------
% Step1 : Useful to the correction application, in the folowing.
% -------------------------------------------------------------------------
% Main profile, both ascending and descending profiles
[~, argo1Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo1Struct.argo, WOA.init, ...
    argo1Struct.argoWork, argo1Struct.Work);

% Secondary profile, both ascending and descending profiles
if argo2Struct.argo.n_prof ~= 0
    if strcmp(Work.whichCorr,'INAIR') && strcmp(argo2Struct.VSS,'Primary sampling')
        argo2Struct.Work.makePlot=1;
        [~, argo2Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo2Struct.argo, WOA.init, ...
        argo2Struct.argoWork, argo2Struct.Work);
        argo2Struct.Work.makePlot=0;
    else
        [~, argo2Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo2Struct.argo, WOA.init, ...
        argo2Struct.argoWork, argo2Struct.Work);
    end
end
if argo3Struct.argo.n_prof ~= 0
    [~, argo3Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo3Struct.argo, WOA.init, ...
        argo3Struct.argoWork, argo3Struct.Work);
end
if argo4Struct.argo.n_prof ~= 0
    [~, argo4Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo4Struct.argo, WOA.init, ...
        argo4Struct.argoWork, argo4Struct.Work);
end

%--MG   Compute reference profile if exist
if isfield(CONFIG,'REF_ARGO')
    [~, argo2Struct.Work.plotREF] = DOXY_interp_REF_main(CONFIG.bddFile, CONFIG.REF_ARGO, argo2Struct.argo, ...
    argo2Struct.argoWork, argo2Struct.Work);
end


% -------------------------------------------------------------------------
% Step2 : To be able to apply correction
% -------------------------------------------------------------------------
% Main profile, ascending profiles only
Work.makePlot=0;
[WOA.interp, Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo, WOA.init, ...
    argoWork, Work);
Work.makePlot=1;

% =========================================================================
%% Get PTS in the primary profile
% PSAL
% the CTD pump is switched off in the near-surface and surface phases.
% Moreover, no PSAL could be measured in the "in-air" phase.
% Therefore, PSAL is get as the first correct value (QC ok) in the primary
% profile between Pcutoff and a maximum depth defined in the configuration
% file (ex: 20m).
%
% TEMP is get as the first correct value (QC ok) in the primary profile
% between Pcutoff and a maximum depth defined in the configuration file
% (ex: 20m).
%
% PRES : the pressure associated with TEMP
% =========================================================================
[argoTrajWork, goProg] = DOXY_get_primary_PTS_for_traj(CONFIG, Work,...
    argoTrajWork,argo1Struct,argo2Struct,argo3Struct, argo4Struct);
if goProg == 0
    return
end

% =========================================================================
%% Get temp and pres from the main profiles
% TEMP is get as the first correct value (QC ok) in the main profile
% between Pcutoff and a maximum depth defined in the configuration file
% (ex: 20m).
% PRES : the pressure associated with TEMP
% Remark : if no data could be find (all NaN), LOCODOX keeps the values
% find just above, in the primary profile.
% =========================================================================
[argoTrajWork, goProg] = DOXY_get_profile_temp_for_traj(CONFIG, Work,...
    argoTrajWork,argoWork,argo);

if goProg == 0
    return
end
% =========================================================================
%% Get doxy from main profile
% =========================================================================
argoTrajWork.profile.doxy_adjusted.data = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.doxy_adjusted.juld = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.cellShape.doxy_adjusted.data = cell(1,length(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.cellShape.doxy_adjusted.juld = cell(1,length(argoTrajWork.ppox_doxy_adjusted.data));

for icycle = 1:length(argoTrajWork.cycle_number.data)
    if ~isempty(argoWork.doxy_adjusted.data(icycle))
        % compare cycle_number
        okCyc = argo.cycle_number.data == argoTrajWork.cycle_number.data(icycle);
        % set psal_adjusted to the primary profil value
        if any(okCyc)
            idx = find(~isnan(argoWork.doxy_adjusted.data(okCyc,:)),1,'first');
            argoTrajWork.profile.doxy_adjusted.data(icycle) = argoWork.doxy_adjusted.data(okCyc,idx);
            argoTrajWork.profile.doxy_adjusted.pres(icycle) = argoWork.pres_adjusted.data(okCyc,idx);
            argoTrajWork.profile.doxy_adjusted.cycle_number(icycle) = argo.cycle_number.data(okCyc);
            argoTrajWork.profile.doxy_adjusted.juld(icycle) = argo.juld.data(okCyc);
            if ~isempty(argoTrajWork.ppox_doxy_adjusted.data{icycle})
                argoTrajWork.profile.cellShape.doxy_adjusted.data{icycle} = repmat(argoWork.doxy_adjusted.data(okCyc,idx),size(argoTrajWork.temp_adjusted.data{icycle}));
                argoTrajWork.profile.cellShape.doxy_adjusted.juld{icycle} = repmat(argo.juld.data(okCyc),size(argoTrajWork.temp_adjusted.data{icycle}));
            end
        end
    end
end

% =========================================================================
%% Compute PO2 from NCEP and Argo
% =========================================================================
% -------------------------------------------------------------------------
% PO2 from Argo : inwater and inair
% -------------------------------------------------------------------------
% 20180731 : dedicated measurement_code for inwater traj PPOX
% New argo format decision : In-Air and near-surface measurements will have
% dedicated measurement_code in the trajectory. If no Near-Surface
% measurement_code available, then take the first correct value from the
% Near-surface profile.
Work.NSIAfloat = 0;
%Save number of measurement in air and in water
Work.nbrMeas.inAir=[];
Work.nbrMeas.inWater=[];
for i = 1:length(argoTrajWork.ppox_doxy_adjusted.data)
    isInWater = ismember(argoTrajWork.measurement_code.data{i},CONFIG.inWaterMC);
    isInAir = ismember(argoTrajWork.measurement_code.data{i},CONFIG.inAirMC);
    if ~isempty(isInWater)
        Work.nbrMeas.inWater=[Work.nbrMeas.inWater length(find(isInWater==1))];
        Work.NSIAfloat = 1;
        Work.surface.argo.pO2_water{i} = argoTrajWork.ppox_doxy_adjusted.data{i}(isInWater);
        Work.surface.argo.temp_water{i} = argoTrajWork.temp_adjusted.data{i}(isInWater);
        Work.surface.argo.juld_water{i} = argoTrajWork.juld.data{i}(isInWater);
    end
    if ~isempty(isInAir)
        Work.nbrMeas.inAir=[Work.nbrMeas.inAir length(find(isInAir==1))];
        Work.surface.argo.pO2_air{i} = argoTrajWork.ppox_doxy_adjusted.data{i}(isInAir);
        Work.surface.argo.temp_air{i} = argoTrajWork.temp_adjusted.data{i}(isInAir);
        Work.surface.argo.juld_air{i} = argoTrajWork.juld.data{i}(isInAir);
    end
end

if Work.NSIAfloat == 0
    Work.surface.argo.pO2_water = {[]};
    Work.surface.argo.temp_water = argoTrajWork.profile.cellShape.temp_adjusted.data;
    Work.surface.argo.juld_water = argoTrajWork.profile.cellShape.doxy_adjusted.juld;
end

if Work.NSIAfloat
    for icycle = 1:length(argoTrajWork.cycle_number.data)
        isInWater = ismember(argoTrajWork.measurement_code.data{icycle},CONFIG.inWaterMC);
        nbr_inWater=length(find(isInWater==1));
        Work.surface.argo.sss{icycle} = repmat(argoTrajWork.profile.psal_adjusted.data(icycle),1,nbr_inWater); 
    end
end
tmpField = fieldnames(Work.surface.argo);
for ifield = 1:length(tmpField)
    Work.surface.argo.(tmpField{ifield}) = cell2mat(cellfun(@(x) double(x),Work.surface.argo.(tmpField{ifield}),'UniformOutput',false));
end

% -------------------------------------------------------------------------
% Switch from cell to double format
% -------------------------------------------------------------------------
trajFields = fieldnames(argoTrajWork);
okfield = ~ismember(trajFields,{'cycle_number','profile'});
trajFields = trajFields(okfield);
argoTrajWork2 = argoTrajWork;

for i=1:length(trajFields)
    tmpField = trajFields{i};
    argoTrajWork2.(tmpField).data = cell2mat(cellfun(@(x) double(x),argoTrajWork.(tmpField).data,'UniformOutput',false));
    argoTrajWork2.(tmpField).dim = {''};
end

varField = fieldnames(argoTrajWork2.profile.cellShape);
for i=1:length(varField)
    profileField = fieldnames(argoTrajWork2.profile.cellShape.(varField{i}));
    for j=1:length(profileField)
        argoTrajWork2.profile.cellShape.(varField{i}).(profileField{j}) = cell2mat(cellfun(@(x) double(x),argoTrajWork2.profile.cellShape.(varField{i}).(profileField{j}),'UniformOutput',false));
    end
end

argoTrajWork = argoTrajWork2;

% -------------------------------------------------------------------------
% PO2 from Argo : inwater and inair if not trueIA
% -------------------------------------------------------------------------
if Work.NSIAfloat == 0
    argoTrajWork.DENS = real(sw_pden(argoTrajWork.profile.cellShape.psal_adjusted.data, ...
        argoTrajWork.profile.cellShape.temp_adjusted.data, argoTrajWork.profile.cellShape.psal_adjusted.pres,0));
    argoTrajWork.inwater.an_dens.data = argoTrajWork.DENS - 1000;
    argoTrajWork.inwater.density.data = argoTrajWork.DENS;
    
    % Saturation and saturation percentage
    unit = argoWork.doxy_adjusted.units;
    tmpDoxy = DOXY_convert(argoTrajWork.profile.cellShape.doxy_adjusted.data,unit,'mumol/L',argoTrajWork.inwater.an_dens.data);
    if Work.presEff == 1, P = argoTrajWork.profile.cellShape.temp_adjusted.pres; else, P = 0; end
    argoTrajWork.inwater.psat.data = O2ctoO2s(tmpDoxy,argoTrajWork.profile.cellShape.temp_adjusted.data,argoTrajWork.profile.cellShape.psal_adjusted.data,P); % in mL/L
    tmpDoxy = DOXY_convert(argoTrajWork.profile.cellShape.doxy_adjusted.data,unit,'ml/L',argoTrajWork.inwater.an_dens.data);
    argoTrajWork.inwater.sat.data = (tmpDoxy./argoTrajWork.inwater.psat.data)*100;
    
    Work.surface.argo.pO2_water = O2stoO2p(argoTrajWork.inwater.psat.data,argoTrajWork.profile.cellShape.temp_adjusted.data,...
        argoTrajWork.profile.cellShape.psal_adjusted.data,P);
end

% -------------------------------------------------------------------------
% Take SSS and SST from "inwater" argo values
% -------------------------------------------------------------------------
% Work.surface.argo.sss = argoTrajWork.psal_adjusted.data;
if ~Work.NSIAfloat
    Work.surface.argo.sss = argoTrajWork.profile.cellShape.psal_adjusted.data;
end

% -------------------------------------------------------------------------
% Scale humidity :
% Scale humidity with logarithmic profile between 10 m NCEP and sea-surface
% water vapor pressure
% -------------------------------------------------------------------------

% Roughness length for humidity (m)
z0q=1e-4;

% Argo sea surface water vapor pressure (mBar)
Work.surface.argo.Pvap = 1013.25*watervapor(Work.surface.argo.temp_water,Work.surface.argo.sss);

% ncep water vapor pressure scaled to the optode heigth (0.2m), in mBar
% WARNING : you have to know the height of the optode. Here, 0.2m
NCEP.coloc.phum = watervapor(NCEP.coloc.air,0).* NCEP.coloc.rhum/100*1013.25;
NCEP.coloc.Pvap = Work.surface.argo.Pvap +(NCEP.coloc.phum - Work.surface.argo.Pvap).*log(.2/z0q)./log(10/z0q);

% -------------------------------------------------------------------------
% PO2 from NCEP
% -------------------------------------------------------------------------
Work.surface.ncep.pO2 = (NCEP.coloc.slp - NCEP.coloc.Pvap).*0.20946;
NCEP.coloc.pO2 = Work.surface.ncep.pO2;

% -------------------------------------------------------------------------
% Quality Control : No negative values, in air values in [180 - 240]
% -------------------------------------------------------------------------
allowedInAirValues = [150 240];
Work.surface.argo.pO2_air(Work.surface.argo.pO2_air < 0) = NaN;
Work.surface.argo.pO2_air(Work.surface.argo.pO2_air < allowedInAirValues(1) | Work.surface.argo.pO2_air > allowedInAirValues(2)) = NaN;

Work.surface.argo.pO2_water(Work.surface.argo.pO2_water < 0) = NaN;

Work.surface.ncep.pO2(Work.surface.ncep.pO2 < 0) = NaN;

NCEP.coloc.pO2(NCEP.coloc.pO2 < 0) = NaN;

%% =========================================================================
% Compute non linear regression between inflated data (PO2 argo in air),
% deflated data (PO2 argo in water subsurface) and air data (PO2 from NCEP)
% Bittig et al., 2015 : 
% Tackling Oxygen Optode Drift: Near-Sruface and
% In-Air Oxygen Optode Measurements on a float Provide and Accurate in
% Situ Reference
% =========================================================================

% =========================================================================
% Check data drift, computed on WOA or NCEP data
% =========================================================================
fprintf('INFO >>>>>    DRIFT : \n')
Work.whichDrift=questdlg('Do you want to compute the drift from WOA data or NCEP data ? ',sprintf('DRIFT COMPUTATION - %d',Work.wmo)','WOA','NCEP','WOA');

if strcmp(Work.whichDrift,'WOA')
    fprintf('\t The drift is computed on WOA data \n')
    vect_WOA=[1 3:(argo2Struct.Dim.n_prof.dimlength)];
    Work_WOA= argo2Struct.Work;
    Work_WOA.whichCorr='REF';
    %--MG
    Work_WOA.makePlot=1;
    Work_WOA.savePlot=Work.savePlot;
    argoWork_WOA=argo2Struct.argoWork;   
    argo_WOA=argo2Struct.argo;
    argo_WOA.juld.data = argo_WOA.juld.data(vect_WOA);
    Work_WOA.DEPTH=Work_WOA.DEPTH(vect_WOA,:);
    Work_WOA.DOXY=Work_WOA.DOXY(vect_WOA,:);
    Work_WOA.DOXY_WOA_interpP=Work_WOA.DOXY_WOA_interpP(vect_WOA,:) ;
    Work_WOA.whichDrift=Work.whichDrift;
    argoWork_WOA.pres_adjusted.data=argoWork_WOA.pres_adjusted.data(vect_WOA,:);
    argoWork_WOA.an_dens.data=  argoWork_WOA.an_dens.data(vect_WOA,:);
    argoWork_WOA.psal_adjusted.data =argoWork_WOA.psal_adjusted.data(vect_WOA,:);
    argoWork_WOA.temp_adjusted.data=argoWork_WOA.temp_adjusted.data(vect_WOA,:);
    Work_WOA.isnotREFcorr=1;
    Work_WOA.min_drift_depth=Work.min_drift_depth;
    [Work_WOA, ~] = DOXY_drift(Work_WOA, argoWork_WOA, argo_WOA,[],'WOA');  
    
    %--MG
    if isfield(Work_WOA,'DOXY_DRIFTCORR_COEF')    
        Work.fitResult_WOA=Work_WOA.DOXY_DRIFTCORR_COEF;
    end
    if isfield(Work_WOA,'ind_drift_stop')
        Work.ind_drift_stop=Work_WOA.ind_drift_stop;
    end
    if Work.savePlot == 1
        Work.savePlot = 0;
        mem = 1;
    else
        mem = 0;
    end
    Work.makePlot=0;
else 
    Work.whichDrift='NCEP';
    mem=Work.savePlot;    
    fprintf('\t The drift is computed on NCEP data \n')
end 

[Work, ~] = DOXY_drift(Work, argoWork, argo, argoTrajWork,Work.whichCorr);
Work.makePlot=1;

% =========================================================================
% Make the non-linear regression
% Po2inflated = c*PO2deflated + [(1-c)/m]PO2air
% =========================================================================
[Work, CORR, hFig, ~] = DOXY_corr_compute_main(Work.PresOrDens, Work.whichCorr, ...
    Work, argoWork, argo, argoTrajWork);

% =========================================================================
% Apply correction to both ascending and descending profile,
% and propagate the correction computed with the main profile
% (ascending) to the secondary profile.
% =========================================================================
for n=1:4
    argoStruct = eval(['argo' num2str(n) 'Struct;']);
    argoStruct.Work.whichO2quantity = Work.whichO2quantity;
    argoStruct.Work.DODRIFT = Work.DODRIFT;
    argoStruct.Work.(['FIT_slope_' Work.whichO2quantity]) = Work.(['FIT_slope_' Work.whichO2quantity]);
    argoStruct.Work.(['FIT_intercept_' Work.whichO2quantity]) = Work.(['FIT_intercept_' Work.whichO2quantity]);
    eval(['argo' num2str(n) 'Struct = argoStruct;']);
end

if mem == 1
    Work.savePlot = 1;
end

%Save drift choice: commented by Thierry Reynaud 12/03/2020
%CORR.whichDrift=Work.whichDrift;
%argo1Struct.Work.whichDrift=Work.whichDrift;
%argo2Struct.Work.whichDrift=Work.whichDrift;
%argo3Struct.Work.whichDrift=Work.whichDrift;


%Save drift choice ->> added by Thierry Reynaud 12/03/2020
CORR.whichDrift=Work.whichDrift;
argo1Struct.Work.whichDrift=Work.whichDrift;
argo2Struct.Work.whichDrift=Work.whichDrift;
argo3Struct.Work.whichDrift=Work.whichDrift;
argo4Struct.Work.whichDrift=Work.whichDrift;


%incoherences
[argo1Struct, argo2Struct,argo3Struct,argo4Struct] = DOXY_corr_apply_main(hFig,...
    CORR, argo1Struct, argo2Struct, argo3Struct, argo3Struct, Work, argoTrajWork);
Work.drift_val=argo2Struct.Work.drift;
