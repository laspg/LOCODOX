% DOXY_corr_compute_main 
%
% SYNTAX
%[Work, CORR, hFig, goProg] = DOXY_corr_compute_main(PresOrDens, whichCorr,...
%    Work, argoWork, argo, argoTrajWork)
%
% DESCRIPTION
% DOXY_corr_compute cover the correction on DOXY and PSAT, following the
% method ATLN.
%
% INPUT
%
%   PresOrDens (string)  Indicate if the correction is applied using
%                        pressure or potential density. 'PRESSURE' or
%                        'DENSITY'.
%
%   whichCorr (string)  Indicate if the correction is made respect to
%                       climatology (WOA), InSitu reference database
%                       (REF) or In-air measurement (INAIR)
%
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA)
%                           Example:
%                           Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                           where
%                           Work.DOXY_WOA
%                                name: 'DOXY_ADJUSTED'
%                                dim: {'N_PROF'  'N_LEVELS'}
%                               data: [85x120 double]
%                          long_name: 'DISSOLVED OXYGEN adjusted and converted'
%                      standard_name: 'moles_of_oxygen_per_unit_mass_in_sea_water'
%                         FillValue_: 99999
%                              units: 'mumol/kg'
%                          valid_min: 0
%                          valid_max: 650
%                           C_format: '%9.3f'
%                     FORTRAN_format: 'F9.3'
%                         resolution: 0.0010
%                               type: 5
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
%                               argo.pres =
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
%   argowork (structure)  float working structure issued from the argo data
%                         Example:
%                         argoWork =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                                   doxySat: [1x1 struct]
%                                  doxyPSat: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%
%   argoTrajWork (structure)   float trajectory working structure
%
% OUTPUT
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the reference data (WOA or REF),
%                        completed with the annual drift.
%                        Example:
%                        Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                                    DERIVE: 4.8400
%
%   CORR (struct)        Correction structure, with main information and
%                        data from the main profile (ascending only,
%                        vertical schemes used for the correction).
%                        Example :
%                             CORR =
%                          presOrDens: 'Pressure'
%                               level: [164x117 single]
%                                cmpl: '_interpP'
%                           whichCorr: 'WOA'
%                              refCyc: [1x164 double]
%                     whichO2quantity: 'PSAT'
%                            strField: 'psat'
%                                psat: [164x117 single]
%                             psatref: [164x117 double]
%                             psat_ok: [164x117 double]
%                                   ...
%
%   hFig (handle)        handle of the figure with correction plots
%                        (DOXY_PLOT_corr).
%
%   goProg (double)      0\1. If this index is 0, LOCODOX stop the
%                        treatment in progress
%
% CALL : 
%   DOXY_corr_compute
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 11/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v3 26/01/2017    Emilie Brion, Anne Piron, Altran Ouest
%                        argo - ameliorations 2017 phase 1
%       v3.2 10/07/2018  Emilie Brion, Altran Ouest
%                        Manage the case with no valid data in PSAT or DOXY
%           28/02/2019   Marine GALLIAN, Altran Ouest
%                        Add reference datetime to gtime
%       v3.4 10.04.2020  T.Reynaud temporal ==> time


function [Work, CORR, hFig, goProg] = DOXY_corr_compute_main(PresOrDens, whichCorr,...
    Work, argoWork, argo, argoTrajWork)

% =========================================================================
%% Initialisation
% =========================================================================
goProg = 1;

CORR.presOrDens = PresOrDens;
if strcmpi(CORR.presOrDens, 'pressure')
    letterPorD = 'P';
    CORR.level = argoWork.pres_adjusted.data;
elseif strcmpi(CORR.presOrDens, 'density')
    letterPorD = 'D';
    CORR.level = Work.DENS;
end

% Initialise informations
CORR.cmpl = sprintf('_interp%s',letterPorD);
CORR.whichCorr = whichCorr;
okdim = strcmp(argoWork.pres_adjusted.dim,'N_PROF');
dim = size(CORR.level);
nbprof = dim(okdim);
if ~isfield(Work,'refCyc')
    CORR.refCyc = 1:nbprof;
else
    CORR.refCyc = Work.refCyc;
    CORR.nbprof = nbprof; %--MG
    Work.nbprof=nbprof;
end

CORR.presEff = Work.presEff;

% Compute %sat or [O2]
if ismember(Work.whichCorr,{'WOA','REF'})
    string = {'The pressure from CTD sensor is used','=> correction using PSAT is recommanded.',' ','Do the correction with :'};
    if Work.betterUseDoxy == 1
        string = {'The pressure from CTD sensor is used','=> correction using PSAT is recommanded.',...
            'BUT some NaNs in salinity have been detected => LOCODOX suggests using DOXY (instead of PSAT) for the oxygen correction'...
            ' ','Do the correction with :'}
    end
 
    CORR.whichO2quantity = questdlg(string,'PARAMETER FOR COMPUTING THE CORRECTION','DOXY','PSAT','PSAT');
    Work.whichO2quantity = CORR.whichO2quantity;
else
    CORR.whichO2quantity = 'PPOX';
    Work.whichO2quantity = CORR.whichO2quantity;
end

% Initialize fields
switch CORR.whichO2quantity
    case 'DOXY'
        CORR.strField = 'doxy';
    case 'PSAT'
        CORR.strField = 'psat';
    case 'PPOX'
        CORR.strField = 'ppox';
end

field = CORR.strField;
if isfield(Work,'DOXY_DRIFTCORR')
    CORR.(field) = Work.([CORR.whichO2quantity '_DRIFTCORR'])(CORR.refCyc,:);
else
    CORR.(field) = Work.(CORR.whichO2quantity)(CORR.refCyc,:);
end

% In case of PSAT entirely NaN, change to DOXY. If DOXY entirely NaN, WARN
% and stop locodox
if ismember(Work.whichCorr,{'WOA','REF'}) && strcmp(CORR.strField,'psat') && all(all(isnan(CORR.(field))))
    string = {'No valid data in PSAT, probably because of PTS with NaN.',...
              '',...
              'Strategy :'...
              '     - Use DOXY',...
              '     - Try to modify the configuration : QC and data mode taking into account',...
              ''};    
    answer = questdlg(string,'ARGO Float DOXY correction','DOXY','Change Config and Exit LOCODOX','DOXY');
    if strcmp(answer,'DOXY')
        CORR.whichO2quantity = 'DOXY';
        CORR.strField = 'doxy';
        Work.whichO2quantity = CORR.whichO2quantity;
    else
        goProg = 0;
        if ~exist('hFig','var')
            hFig = 0;  % value defined just to have hfig defined (case: makePlot = 0)
        end
        return
    end
    field = CORR.strField;
    if isfield(Work,'DOXY_DRIFTCORR')
        CORR.(field) = Work.([CORR.whichO2quantity '_DRIFTCORR'])(CORR.refCyc,:);
    else
        CORR.(field) = Work.(CORR.whichO2quantity)(CORR.refCyc,:);
    end
    if all(all(isnan(CORR.(field))))
        warndlg('No valid data in DOXY !!','DOXY_corr_compute_main');
        goProg = 0;
        return;
    end       
end

% Initialize work arrays
if strcmp(Work.whichCorr,'WOA') || strcmp(Work.whichCorr,'REF')
    CORR.([field 'ref']) = Work.([CORR.whichO2quantity '_' CORR.whichCorr CORR.cmpl]);
    CORR.([field '_ok']) = NaN(length(CORR.refCyc),dim(~okdim));       % argo data, with bad data put to NaN
    CORR.([field 'ref_ok']) = NaN(length(CORR.refCyc),dim(~okdim));    % WOA data, with bad data put to NaN
    CORR.([field '_ko']) = NaN(length(CORR.refCyc),dim(~okdim));       % rejected argo data
    CORR.([field 'ref_ko']) = NaN(length(CORR.refCyc),dim(~okdim));    % rejected WOA data
    gradO2 = NaN(length(CORR.refCyc),dim(~okdim));                     % d02/dz or dPSAT/dz
elseif strcmp(Work.whichCorr,'INAIR')
    % Keep the common cycle
    isok = ismember(argoTrajWork.cycle_number.data,argo.cycle_number.data);
    
    % Variables
    CORR.presOrDens = Work.PresOrDens;
    if isfield(Work,'PPOX_DRIFTCORR_inair')
        CORR.po2_inwater = Work.PPOX_DRIFTCORR_inwater;
        CORR.po2_inair = Work.PPOX_DRIFTCORR_inair;
    else
        CORR.po2_inwater = Work.surface.argo.pO2_water;
        CORR.po2_inair = Work.surface.argo.pO2_air;
    end
    
    % -- build a po2_inwater the same size than po2_inair, for correction
    % computation     
    CORR.po2_inwater_sizeInair = CORR.po2_inwater;        
    CORR.ncep_pO2_scaled = Work.surface.ncep.pO2;
    % --MG
    if ~Work.NSIAfloat   
        CORR.po2_inair(CORR.po2_inair < 0) = NaN;
        CORR.po2_inwater_gtime = argoTrajWork.juld.data + ...
                datenum(strrep(argoTrajWork.juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');
        CORR.gtime = CORR.po2_inwater_gtime+datenum(strrep(argoTrajWork.juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');       
    else
        CORR.po2_inwater_gtime = Work.surface.argo.juld_water + datenum(strrep(argoTrajWork.juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');
        CORR.gtime = Work.surface.argo.juld_air + datenum(strrep(argoTrajWork.juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');
    end
end

% units
CORR.unit = argoWork.doxy_adjusted.units;

% figure
if Work.makePlot
    if Work.DODRIFT
        cmpl = ', time drift correction applied';
    else
        cmpl = ', time drift correction NOT applied';
    end
    hFig = figure('Name',sprintf('CORRECTION DOXY and PSAT - %d %s',Work.wmo,cmpl),'NumberTitle','off',...
        'unit','normalized','OuterPosition',[0.00 0.05 0.52 0.61]);
end
% =========================================================================
%% Compute correction
% WOA, INSITU REF : Find the linear regression : doxy VS doxy WOA
% INAIR correction : Find the non-linear regression to correct PPOX
% =========================================================================
if ~exist('hFig','var')
    hFig = 0;  % value defined just to have hfig defined (case: makePlot = 0)
end

[Work, CORR, hFig] = DOXY_corr_compute(Work, CORR, hFig);
