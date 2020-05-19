% DOXY_argoTraj_read read the argo trajectory data
%
% SYNTAX
%[argoTraj, argoTrajC, argoTrajWork, Dim] = DOXY_argoTraj_read(CONFIG,
%                                                   wmo)
%
% DESCRIPTION
% DOXY_argoTraj_read        Read the argo trajectory data
%
% INPUT
%     CONFIG (struct)       Configuration structure with data path,
%                           operator choices, ...
%                             CONFIG = 
%                        DataDir: '/home/oo26/coriolis/co05/co0508/dac/coriolis/'
%                    NCEPDataDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/data_input/'
%                    saveDataDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/results/data/'
%                    savePlotDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/results/plots/'
%                              ...
%                              ...
%                        inAirMC: [699 711 799]                         
%                           
%      wmo (double)         WMO of the float 
%
% OUTPUT
%     argoTraj                    float trajectory data structure (argo)
%                                 
%
%     argoTrajWork (structure)    float trajectory intermediate working
%                                 structure
%
%     argoTrajC (structure)       float trajectory CORE data structure
%
%     Dim(structure)              float Dim structure (Dim) choose for
%                                 the Doxy correction.
%
% CALL : 
%   NCR_file, DOXY_get_ppox
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 10/05/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%   v1.2 05/07/2018  Emilie Brion, Altran Ouest
%                    ppox_doxy field in argoTrajWork : get
%                    ppox_doxy_adjusted or ppox_doxy

function [argoTraj, argoTrajC, argoTrajWork, Dim] = ...
    DOXY_argoTraj_read(CONFIG,wmo)

% =========================================================================
%% Read trajectory data
% =========================================================================
trajDir = fullfile(CONFIG.DataDir,num2str(wmo));

bio_file = dir(fullfile(trajDir,'*B*traj.nc'));
[argoTraj,Dim] = NCR_file(fullfile(trajDir,bio_file.name));
argoTraj.cycle_number.data = double(argoTraj.cycle_number.data);
argoTraj.cycle_number.type = 6;
argoTraj = fillValue_2_nan(argoTraj);

core_file = dir(fullfile(trajDir,'*_Rtraj.nc'));
[argoTrajC] = NCR_file(fullfile(trajDir,core_file.name));
argoTrajC.cycle_number.data = double(argoTrajC.cycle_number.data);
argoTrajC.cycle_number.type = 6;
argoTrajC = fillValue_2_nan(argoTrajC);

argoTrajWork.juld = argoTrajC.juld;
argoTrajWork.latitude = argoTrajC.latitude;
argoTrajWork.longitude = argoTrajC.longitude;
argoTrajWork.cycle_number = argoTrajC.cycle_number;
argoTrajWork.measurement_code = argoTrajC.measurement_code;

% Get rid off position_qc >= 3
isnok = argoTrajC.position_qc.data(1:length(argoTraj.longitude.data)) >= '3';
argoTrajWork.latitude.data(isnok) = NaN;
argoTrajWork.longitude.data(isnok) = NaN;
argoTrajWork.juld.data(isnok) = NaN;

% Manage PRES, TEMP and PSAL : find VAR_DOXY prior to VAR_ADJUSTED prior to
% VAR. Look first in the Bio traj file, then in the Core traj file. If the
% VAR_DOXY or the VAR_ADJUSTED is filled with NaN, go down in the priority.
myVar = {'pres','temp','psal'};
argoTrajWork.trajVar.type.(myVar{1}) = 1;
argoTrajWork.trajVar.where.(myVar{1}) = 1;

for v = 1 : length(myVar)
    varadj = sprintf('%s_adjusted',myVar{v});
    if isfield(argoTrajC,varadj) 
        if any(~isnan(argoTrajC.(varadj).data))
            varAdjInCoreTraj = true;
        else
            varAdjInCoreTraj = false;
        end
    end
    if varAdjInCoreTraj    
        argoTrajWork.(varadj) = argoTrajC.(varadj);
        argoTrajWork.trajVar.type.(myVar{v}) = 1;
        argoTrajWork.trajVar.where.(myVar{v}) = 2;
    else
        argoTrajWork.(varadj) = argoTrajC.(myVar{v});
        argoTrajWork.trajVar.type.(myVar{v}) = 2;
        argoTrajWork.trajVar.where.(myVar{v}) = 2;
    end    
end
argoTrajWork.trajVar.type.help = '1: var _adjusted, 2 : raw var';
argoTrajWork.trajVar.where.help = '1 : from traj Bio, 2: from traj Core, 3 : from Near-Surface profile, 4: from primary profile, 5: from secondary profile';

% =========================================================================
%% Choose adjusted field or not
% =========================================================================
%--------------------------------------------------------------------------
% PPOX is carried by MC = CONFIG.inAirMC, and has to be dated by Ascent End
% Time
%
% Which doxy to use ? DOXY or DOXY_ADJUSTED ?
% - any(doxy_adjusted) ~= NaN and/or any(doxy_adjusted_qc) ~= 4 => select
% - any(doxy_adjusted) ~= NaN) but all(doxy_adjusted_qc) == 4   => don't
% select and keep the real time DOXY and DOXY_QC
% - if DOXY and DOXY_ADJUSTED is all nan and/or qc=4, warning : no data !!
%--------------------------------------------------------------------------

%Searching for IN AIR measurements
trajFields = fieldnames(argoTraj);
isok = ~cellfun('isempty',strfind(trajFields,'ppox_doxy_adjusted'));

if any(isok) && ~isempty(str2num(argoTraj.ppox_doxy_adjusted_qc.data')) %#ok<ST2NM>
    doxyField = 'ppox_doxy_adjusted';
else
    isok = ~cellfun('isempty',strfind(trajFields,'ppox_doxy'));
    if any(isok) && ~isempty(str2num(argoTraj.ppox_doxy_qc.data')) %#ok<ST2NM>
        doxyField = 'ppox_doxy';
    else
        warning('No PPOX data, not in real nor in adjusted fields !!!')
        return        
    end
end

argoTrajWork.ppox_doxy_adjusted = argoTraj.(doxyField);
argoTrajWork.ppox_doxy_adjusted_qc = argoTraj.([doxyField '_qc']);
argoTrajWork.ppox_doxy_adjusted.data = NaN(size(argoTrajWork.juld.data));
argoTrajWork.ppox_doxy_adjusted_qc.data = repmat(argoTraj.ppox_doxy_qc.FillValue_,size(argoTrajWork.juld.data));
argoTrajWork.ppox_doxy_adjusted.data(1:length(argoTraj.juld.data)) = argoTraj.(doxyField).data;
argoTrajWork.ppox_doxy_adjusted_qc.data(1:length(argoTraj.juld.data)) = argoTraj.([doxyField '_qc']).data;

% =========================================================================
%% Don't take into account the cycle -1
% If cycle #0 has no measurement_code 1099/1090 (surface drifting measurement)
% and 703 (Ascent End Time), don't take it into account
% =========================================================================
isok = argoTrajWork.cycle_number.data ~= -1;
argoTrajWork.juld.data = argoTrajWork.juld.data(isok);
argoTrajWork.latitude.data = argoTrajWork.latitude.data(isok);
argoTrajWork.longitude.data = argoTrajWork.longitude.data(isok);
argoTrajWork.cycle_number.data = argoTrajWork.cycle_number.data(isok);
argoTrajWork.measurement_code.data = argoTrajWork.measurement_code.data(isok);
argoTrajWork.pres_adjusted.data = argoTrajWork.pres_adjusted.data(isok);
argoTrajWork.temp_adjusted.data = argoTrajWork.temp_adjusted.data(isok);
argoTrajWork.psal_adjusted.data = argoTrajWork.psal_adjusted.data(isok);
argoTrajWork.ppox_doxy_adjusted.data = argoTrajWork.ppox_doxy_adjusted.data(isok);
argoTrajWork.ppox_doxy_adjusted_qc.data = argoTrajWork.ppox_doxy_adjusted_qc.data(isok);

% =========================================================================
%% Don't get the cycle 0 data : possibly on board, ... data not confident
% =========================================================================
isok = argoTrajWork.cycle_number.data ~= 0;
argoTrajWork.juld.data = argoTrajWork.juld.data(isok);
argoTrajWork.latitude.data = argoTrajWork.latitude.data(isok);
argoTrajWork.longitude.data = argoTrajWork.longitude.data(isok);
argoTrajWork.cycle_number.data = argoTrajWork.cycle_number.data(isok);
argoTrajWork.measurement_code.data = argoTrajWork.measurement_code.data(isok);
argoTrajWork.pres_adjusted.data = argoTrajWork.pres_adjusted.data(isok);
argoTrajWork.temp_adjusted.data = argoTrajWork.temp_adjusted.data(isok);
argoTrajWork.psal_adjusted.data = argoTrajWork.psal_adjusted.data(isok);
argoTrajWork.ppox_doxy_adjusted.data = argoTrajWork.ppox_doxy_adjusted.data(isok);
argoTrajWork.ppox_doxy_adjusted_qc.data = argoTrajWork.ppox_doxy_adjusted_qc.data(isok);

% =========================================================================
%% Keep only data with QC allowed by operator
% (QC defined in the configuration file)
% =========================================================================
tabdoxyqc = argoTrajWork.ppox_doxy_adjusted_qc.data;
tabdoxyqc(tabdoxyqc == ' ') = '9';
tabdoxyqc = str2num(tabdoxyqc'); 

QCp = 1:1:9;       % all possible QCs
QCnan = NaN(1,length(QCp)-length(CONFIG.QC_O));
cmpt = 0;
for i = 1:length(QCp)
    isok = find(QCp(i) == CONFIG.QC_O, 1);
    if isempty(isok)
        cmpt = cmpt + 1;
        QCnan(cmpt) = QCp(i);
    end
    clear isok
end
for i = 1:length(QCnan)
    argoTrajWork.ppox_doxy_adjusted.data(tabdoxyqc==num2str(QCnan(i))) = NaN;
end

% % =========================================================================
% %% Compute doxy sat, psat and ppox
% % Compute doxy psat from ppox
% % WARNING : at this level, the psat is uncorrect as the psal is not the one
% % from the primary profile. it is just to initialisz
% % =========================================================================
% % argoTrajWork.psat.data = O2ctoO2s(argoTrajWork.doxy_adjusted.data,...
% %     argoTrajWork.temp_adjusted.data,argoTrajWork.psal_adjusted.data);
% % 
% % argoTrajWork.ppox.data = O2stoO2p(argoTrajWork.psat.data,...
% %     argoTrajWork.temp_adjusted.data,argoTrajWork.psal_adjusted.data);
% 
% argoTrajWork.psat.data = O2ptoO2s(argoTrajWork.ppox_doxy_adjusted.data,...
%     argoTrajWork.temp_adjusted.data,argoTrajWork.psal_adjusted.data);

% =========================================================================
%% Select : juld, position and doxy of PPOX
% =========================================================================
% Ascent end time kept if no explicite datation
% Lat/lon from the first argos localisation available
% DOXY : last point of the doxy profile kept as PPOX. WARNING : evolution :
% all the unpumped DOXY data will be kept as PPOX in the traj file, dated
% with AET ?
% [argoTrajWork] = DOXY_get_ppox(argoTrajWork,Dim,CONFIG,argo1Struct,argo2Struct,argo3Struct, argo4Struct);
[argoTrajWork] = DOXY_get_ppox(argoTrajWork,Dim,CONFIG);

% % =========================================================================
% %% Get Psal in the primary profile
% % the CTD pump is switched off in the near-surface and surface phases.
% % Moreover, no PSAL could be measured in the "in-air" phase.
% % Therefore, PSAL is get as the first correct value (QC ok) in the primary
% % profile between Pcutoff and a maximum depth defined in the configuration
% % file (ex: 20m).
% %
% % A VOIR : si jamais les profileurs font des mesures de salinité dans leur
% % phase de dérive de subsurface (avant coup de pompe = emersion puis mesure
% % dans l'air "vraie"), valeurs à prendre ? je ne pense pas car pompe
% % arrêtée...
% % =========================================================================
% [argoTrajWork, goProg] = DOXY_get_psal_for_traj(CONFIG, REF_ARGO,...
%                         argoTrajWork, argo1Struct,argo2Struct,argo3Struct, argo4Struct);
% if goProg == 0
%     return
% end

% Compute doxy psat
%--EB 20180423
% PSAT should be computed from data at the same time in water or at
% the same time in air. Here, salinity is always measured in water, so the
% psat should be computed only for the in-water part.
% => the psat computation has no sense here for the moment
% for iCyc = 1:length(argoTrajWork.psal_adjusted.data)
%     if ~isempty(argoTrajWork.psal_adjusted.data{iCyc})
%         argoTrajWork.psat.data{iCyc} = O2ptoO2s(argoTrajWork.ppox_doxy_adjusted.data{iCyc},...
%             argoTrajWork.temp_adjusted.data{iCyc},argoTrajWork.psal_adjusted.data{iCyc});
%     end
% end
%--EB