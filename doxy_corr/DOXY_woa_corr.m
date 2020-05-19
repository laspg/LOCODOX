% DOXY_woa_corr do the correction of DOXY argo data using the WOA
% climatology.
%
% SYNTAX
% [argo1Struct,argo2Struct, argo3Struct, argo4Struct, WOA, Work, goProg] = ...
%                   DOXY_inair_corr(CONFIG,WOA,argo1Struct,
%                   argo2Struct, argo3Struct, argo4Struct, argo, argoWork, Work)
%
% DESCRIPTION
% DOXY_woa_corr do the DOXY argo data correction, through different steps :
% * interpolate WOA data to argo levels and position
% * Check data drift between the depth chosen by the user and 6000m
% * Compute correction
% * Apply correction to both ascending and descending profile, and
% propagate the correction computed with the main profile (ascending) to
% the other vertical sampling schemes profile.
% * Apply quality flags to DOXY
% Remarque: For WOA correction, the main profile could be the Secondary
% vertical sampling scheme (the priority), or the Primary Sampling scheme
% if DOXY is not present in the Secondary sampling scheme.
%
% INPUT
%   CONFIG (struct)       Configuration Structure
%   
%   WOA (struct)          Climatology WOA structure.
%                         Example:
%                          WOA = 
%                                init: [1x1 struct]
%                              interp: [1x1 struct]
%                          WOA.init 
%                                dimorder: 'C'
%                                latitude: [1x1 struct]
%                               longitude: [1x1 struct]
%                                   depth: [1x1 struct]
%                                    time: [1x1 struct]
%                                 doxywoa: [1x1 struct]
%                                           ...
%       
%   argo1Struct (struct)    Three data structures organised as the main profile
%   argo2Struct (struct)    for working is in argo1Struct (vertical
%   argo3Struct (struct)    sampling shceme of interest for the kind of
%   argo4Struct (struct)    correction), and the other in the second ones.
%                           The program DOXY_argo_prepare_main.m describes
%                           the choice of the main and other profiles.
%                           Example:
%                             argo1Struct = 
%                                     argo: [1x1 struct]
%                                      Dim: [1x1 struct]
%                                     Work: [1x1 struct]
%                                 argoWork: [1x1 struct]
%                                      VSS: 'Secondary sampling'
%                       
%   argo (structure)      Argo float structure (get directly from data),
%                         for ascending profile and main profile.
%                           Example :
%                               argo =
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
%    Work (structure)     Doxy correction working structure, issued and
%                         computed from argo float data.
%                         Example:
%                             Work = 
%                                   readme: [1x1 struct]
%                                     unit: 'mumol/kg'
%                                      wmo: 1901205
%                                 DOXY_RAW: [85x120 single]
%                                  CAPTEUR: {'Optode'}
%                                  DOXY_QC: [85x120 char]
%
%
% OUTPUT
%    argo1Struct (struct)  argo1Struct, argo2Struct and argo3Struct completed
%    argo2Struct
%    argo3Struct
%    argo4Struct
%
%    WOA (struct)          WOA.interp completed
%
%    Work (struct)         Work completed, in particular with coefficient
%                          of DOXY linear correction.
%
%    goProg (double)       0\1. If this index is 0, LOCODOX stop the
%                          treatment in progress
%
% CALL : 
%   DOXY_interp_WOA_main, DOXY_drift, DOXY_corr_compute_main,
%   DOXY_corr_apply_main, DOXY_interp_REF_main
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 03/05/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%            25/03/2019         Marine GALLIAN, Altran Ouest 
%                               Is reference profil exists, recover it in
%                               order to plot it with corrected data after
%            05/07/19           %Marine GALLIAN, Altran Ouest
%                               Add the possibility to compute the drift on
%                               NCEP data if there is PPOX data
%            12/03/2020         Thierry Reynaud
%                               whichDrift was missing for argoXStruct

function [argo1Struct,argo2Struct,argo3Struct,argo4Struct,WOA, Work,goProg] = DOXY_woa_corr(CONFIG, ...
                WOA, argo1Struct, argo2Struct, argo3Struct,argo4Struct, argo, argoWork, Work)
               
% =========================================================================
% Interpolate WOA data over the argo positions and its vertical
% grid
% =========================================================================
% -------------------------------------------------------------------------
% Step1 : Useful to the correction application, in the folowing.
% -------------------------------------------------------------------------
% Main profile, both ascending and descending profiles
[~, argo1Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo1Struct.argo, WOA.init, ...
    argo1Struct.argoWork, argo1Struct.Work);

% Secondary profiles, both ascending and descending profiles
if argo2Struct.argo.n_prof ~= 0
    [~, argo2Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo2Struct.argo, WOA.init, ...
        argo2Struct.argoWork, argo2Struct.Work);
end
if argo3Struct.argo.n_prof ~= 0
    [~, argo3Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo3Struct.argo, WOA.init, ...
        argo3Struct.argoWork, argo3Struct.Work);
end
if argo4Struct.argo.n_prof ~= 0
    [~, argo4Struct.Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo4Struct.argo, WOA.init, ...
        argo4Struct.argoWork, argo4Struct.Work);
end


%--MG Recover reference profil if exists
if isfield(CONFIG,'REF_ARGO')
    [~, argo1Struct.Work.plotREF] = DOXY_interp_REF_main(CONFIG.bddFile, CONFIG.REF_ARGO, argo1Struct.argo, ...
    argo1Struct.argoWork, argo1Struct.Work);
end

% -------------------------------------------------------------------------
% Step2 : To be able to apply correction
% -------------------------------------------------------------------------
% Main profile, ascending profiles only
[WOA.interp, Work] = DOXY_interp_WOA_main(CONFIG.maskFile, CONFIG.varMask, argo, WOA.init, ...
    argoWork, Work);


% =========================================================================
% Compute drift : Check data drift between the depth chosen by the user and 6000m
% =========================================================================
%marine 20/06/19
fprintf('INFO >>>>>    DRIFT : \n')
if isfield(Work,'NCEP_drift')
    Work.whichDrift=questdlg('Do you want to compute the drift from WOA data or NCEP data ? ',sprintf('DRIFT COMPUTATION - %d',Work.wmo)','WOA','NCEP','WOA');
    if strcmp(Work.whichDrift,'NCEP')
         fprintf('\t The drift is computed on NCEP data \n')   
        Work.NCEP_drift.Work.whichDrift= Work.whichDrift;
        Work.NCEP_drift.Work.isnotINAIRcorr=1;
        Work.NCEP_drift.Work.DOXY_convert = Work.DOXY;
        % Variables necessaires pour appliquer la dérive à DOXY et PSAT 
        Work.NCEP_drift.Work.temp_convert = argoWork.temp_adjusted.data;
        Work.NCEP_drift.Work.psal_convert =argoWork.psal_adjusted.data;
        Work.NCEP_drift.Work.pres_convert =argoWork.pres_adjusted.data; 
        Work.NCEP_drift.Work.units_convert=argoWork.doxy_adjusted.units;
        Work.NCEP_drift.Work.dens_convert=argoWork.an_dens.data;
        [Work.NCEP_drift.Work_NCEP, ~] = DOXY_drift(Work.NCEP_drift.Work,Work.NCEP_drift.argoWork, Work.NCEP_drift.argo, Work.NCEP_drift.argoTrajWork, Work.NCEP_drift.Work.whichCorr);
        
        Work.whichDrift = 'NCEP';
        Work.DODRIFT=Work.NCEP_drift.Work_NCEP.DODRIFT;
        if Work.DODRIFT==1
            Work.DOXY_DRIFTCORR=Work.NCEP_drift.Work_NCEP.DOXY_DRIFTCORR;
            Work.PSAT_DRIFTCORR=Work.NCEP_drift.Work_NCEP.PSAT_DRIFTCORR;
            Work.PPOX_DRIFTCORR_COEF=Work.NCEP_drift.Work_NCEP.PPOX_DRIFTCORR_COEF;
            if isfield(Work.NCEP_drift.Work_NCEP,'ind_drift_stop')
                Work.ind_drift_stop=Work.NCEP_drift.Work_NCEP.ind_drift_stop;
            end
        end
        Work.makePlot=1;
        
    else
        Work.whichDrift = 'WOA';
        fprintf('\t The drift is computed on WOA data \n')   
        [Work, ~] = DOXY_drift(Work, argoWork, argo, [], Work.whichCorr);
    end
    
else
    Work.whichDrift = 'WOA';        
    fprintf('\t The drift is computed on WOA data \n')
    [Work, ~] = DOXY_drift(Work, argoWork, argo, [], Work.whichCorr);
end


% =========================================================================
% Compute correction
% =========================================================================
if Work.savePlot == 1
    Work.savePlot = 0;
    mem = 1;
else
    mem = 0;
end

[Work, CORR, hFig, goProg] = DOXY_corr_compute_main(Work.PresOrDens, Work.whichCorr, Work, argoWork,[], []);

if goProg == 0
    return
end

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
    argoStruct.Work.whichDrift=Work.whichDrift;
    if isfield(Work,'ind_drift_stop')
        argoStruct.Work.ind_drift_stop=Work.ind_drift_stop;
    end
    eval(['argo' num2str(n) 'Struct = argoStruct;']);
end

if mem == 1
    Work.savePlot = 1;
end

%Save drift choice: commented by Thierry Reynaud 12/03/2020
%CORR.whichDrift=Work.whichDrift;

%Save drift choice ->> added by Thierry Reynaud 12/03/2020
CORR.whichDrift=Work.whichDrift;
argo1Struct.Work.whichDrift=Work.whichDrift;
argo2Struct.Work.whichDrift=Work.whichDrift;
argo3Struct.Work.whichDrift=Work.whichDrift;
argo4Struct.Work.whichDrift=Work.whichDrift;

[argo1Struct, argo2Struct,argo3Struct,argo4Struct] = ...
    DOXY_corr_apply_main(hFig, CORR, argo1Struct, argo2Struct, argo3Struct, argo4Struct, Work);
Work.drift_val=argo1Struct.Work.drift;