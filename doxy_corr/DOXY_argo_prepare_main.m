% DOXY_argo_prepare_main ...
%
% SYNTAX
% [DATA, argo1Struct, argo2Struct, argo3Struct, argoWork, Work] = ...
%                 DOXY_argo_prepare_main(DATA, Work, argo)
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
%     DATA (structure)    argo float "multiprofile" structure gathering the
%                         three kind of vertical sampling scheme. Subfields
%                         are argo (the data) and Dim (the dimension).
%                         Example :
%                             DATA = 
%                                   primary: [1x1 struct]
%                                  nearsurf: [1x1 struct]
%                                 secondary: [1x1 struct]
%                                     other: [1x1 struct]
%
%                             DATA.primary =
%                                 argo: [1x1 struct]
%                                  Dim: [1x1 struct]
%
%                             DATA.primary.argo =
%                                   dimorder: 'C'
%                                  data_type: [1x1 struct]
%                                          ...
%                                       pres: [1x1 struct]
%                                       doxy: [1x1 struct]
%                                    doxy_qc: [1x1 struct]
%                              doxy_adjusted: [1x1 struct]
%                           doxy_adjusted_qc: [1x1 struct]
%                        doxy_adjusted_error: [1x1 struct]
%
%     argo (structure)      Argo float structure (get directly from data), 
%                           selected for correction.
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
% OUTPUT
%     argo (structure)      float data structure (argo) structure chosen
%                           for the Doxy correction.
%
%     argowork (structure)  Float working structure, issued and computed
%                           from argo float data.
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
%     argo1Struct (struct)  Data structure, taken from DATA.(field), to
%     argo2Struct (struct)  identify the main profile (argo1Struct)
%     argo3Struct (struct)  according to the kind of correction.
%     argo4Struct (struct)
%
% CALL : 
%   extract_profile_dim, DOXY_get_PTS, DOXY_argo_prepare
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 26/05/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%              v3 26/01/2017   Anne Piron, Altran Ouest
%                              argo - ameliorations 2017 phase 1
%                 22/02/2019    Marine GALLIAN, Altran Ouest 
%                               Now, when IN AIR correction is used, 
%                               primary sampling and near surface data are
%                               displayed instead of only Near surface data
%

function [argo1Struct, argo2Struct, argo3Struct, argo4Struct, argoWork, Work] = ...
                DOXY_argo_prepare_main(DATA, Work, argo)

% =========================================================================
%% Prepare the argo data, in particular undertakes a quality control of
% doxy data, applying flags.
% =========================================================================
if ~isfield(Work,'is_not_inair_corr')
    fprintf('INFO >>>>>    Prepare data : reference profile is %s\n',argo.VSS);
end
if strcmp(argo.VSS,'Secondary sampling') && ~strcmp(Work.whichCorr,'INAIR')
    h = warndlg({sprintf('Reference profile is Secondary Sampling !'),...
        ' ',...
        'Since several data points may be lost during the conversion from DOXY to PSAT, LOCODOX suggests using the DOXY-based correction.',...
        });
    uiwait(h);
end
% Indicator: PTS data checked if set to 1, not checked elsewhere
Work.checkPTS = 1;
isAprofile  = strcmp(cellstr(DATA.primary.argo.direction.data),'A');
[argoP,Dim] = extract_profile_dim(DATA.primary.argo,DATA.primary.Dim,'N_PROF',isAprofile);
[argoWorkP,Work] = DOXY_get_PTS(argoP,Work,Work.dirLog);

if Work.makePlot==1; Work.makePlot=0; mem_plot=1; else ; mem_plot=0; end

%Compute work structure, do not plot data 
if strcmp(argo.VSS,'Secondary sampling')
    [argoWork,Work] = DOXY_argo_prepare(argo,Work,Work.dirLog,0,argoP,argoWorkP);
else
    [argoWork,Work] = DOXY_argo_prepare(argo,Work,Work.dirLog,0);
end

if mem_plot==1; Work.makePlot=1; end

VSSname = {'primary','secondary','nearsurf','other'};
for iVSS = 1:length(VSSname)
    if iVSS==3 && ~strcmp(Work.whichCorr,'INAIR')
        DATA.(VSSname{iVSS}).Work.checkPTS = Work.checkPTS;
        DATA.(VSSname{iVSS}).Work.makePlot = 0;
    else
        DATA.(VSSname{iVSS}).Work.checkPTS = Work.checkPTS;
        if mem_plot==1
            DATA.(VSSname{iVSS}).Work.makePlot = 1;
        end
    end
end

if Work.savePlot == 1
    Work.savePlot = 0;
    mem = 1;
else
    mem = 0;
end

% Prepare primary data if there is at least one profile
if DATA.primary.argo.n_prof ~= 0
    [DATA.primary.argoWork,DATA.primary.Work] = DOXY_argo_prepare(DATA.primary.argo,...
                    DATA.primary.Work,Work.dirLog,1);
end


% Lignes décommentées : %marine 17/06/19
argoP = DATA.primary.argo;
[argoWorkP,Work] = DOXY_get_PTS(DATA.primary.argo,Work,Work.dirLog);

% Prepare secondary data if there is at least one profile
if strcmp(argo.VSS,'Secondary sampling') && DATA.secondary.argo.n_prof ~= 0
    [DATA.secondary.argoWork,DATA.secondary.Work] = DOXY_argo_prepare(DATA.secondary.argo,...
        DATA.primary.Work,Work.dirLog,2,argoP,argoWorkP);
elseif DATA.secondary.argo.n_prof ~= 0
    [DATA.secondary.argoWork,DATA.secondary.Work] = DOXY_argo_prepare(DATA.secondary.argo,...
        DATA.primary.Work,Work.dirLog,2,argoP,argoWorkP);
end


% Prepare near surface data if there is at least one profile
if DATA.nearsurf.argo.n_prof ~= 0
    if ismember(Work.whichCorr,{'WOA','REF'})
        [DATA.nearsurf.argoWork,DATA.nearsurf.Work] = DOXY_argo_prepare(DATA.nearsurf.argo,...
                        DATA.nearsurf.Work,Work.dirLog,3,argoP,argoWorkP);
    else
        [DATA.nearsurf.argoWork,DATA.nearsurf.Work] = DOXY_argo_prepare(DATA.nearsurf.argo,...
                        DATA.nearsurf.Work,Work.dirLog,3);
    end        
end

%Prepare other data if there is at least one profile
if DATA.other.argo.n_prof ~= 0
    if ismember(Work.whichCorr,{'WOA','REF'})
        %WARNING : nearsurf : if PSAL and PSAL_ADJ is NaN, get the
        %psal/temp/pres from primary profile => la densité calculée (et
        %donc le PSAT) sera biaisée par rapport aux profondeurs du profil
        %Near-Surface.
        [DATA.other.argoWork,DATA.other.Work] = DOXY_argo_prepare(DATA.other.argo,...
                       DATA.other.Work,Work.dirLog,4,argoP,argoWorkP);
    else
        [DATA.other.argoWork,DATA.other.Work] = DOXY_argo_prepare(DATA.other.argo,...
                       DATA.other.Work,Work.dirLog,4);    
    end
end



if mem == 1; Work.savePlot = 1; end

DATA.primary.VSS = DATA.primary.argo.VSS;
DATA.nearsurf.VSS = DATA.nearsurf.argo.VSS;
DATA.secondary.VSS = DATA.secondary.argo.VSS;
DATA.other.VSS = DATA.other.argo.VSS;

for iVSS = 1:length(VSSname)
    DATA.(VSSname{iVSS}).Work.makePlot = 0;
end

% ==================================================================================
%% Identify the main profile, and order the other ones
%  (near-surface and secondary profile)
% For the WOA and REF correction, data of interest are the pumped data. The
% secondary profile is taken prior to the primary.
%       - 1. Secondary profile > Primary profile (if data in secondary)
%       - 2. remaining pumped data
%       - 3. unpumped data (Near-Surface)
% 
% For the INAIR correction, the main profile is the near-surface sampling
% scheme.
%       - 1. Near-Surface (unpumped)
%       - 2. Secondary profile > Primary profile  (if data in secondary)
%       - 3. remaining pumped data
% =========================================================================


if ismember(Work.whichCorr,{'WOA','REF'})
    if any(DATA.primary.argo.hereDoxy.data)
        argo1Struct = DATA.primary;
        argo2Struct = DATA.secondary;
    elseif any(DATA.secondary.argo.hereDoxy.data)
        argo1Struct = DATA.secondary;
        argo2Struct = DATA.primary;
    end
    argo3Struct = DATA.nearsurf;
    argo4Struct = DATA.other;
    
elseif strcmp(Work.whichCorr,'INAIR')
    argo1Struct = DATA.nearsurf;
    if any(DATA.primary.argo.hereDoxy.data)
        argo2Struct = DATA.primary;
        argo3Struct = DATA.secondary;
    elseif any(DATA.secondary.argo.hereDoxy.data)
        argo2Struct = DATA.secondary;
        argo3Struct = DATA.primary;
    end
    argo4Struct = DATA.other;
end

end