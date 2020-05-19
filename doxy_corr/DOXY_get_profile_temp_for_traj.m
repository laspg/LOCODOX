% DOXY_get_profile_temp_for_traj gets the TEMP data in the main profile
% usesul for inAir correction from trajectory data.
%
% SYNTAX
% [argoTrajWork] = DOXY_get_profile_temp_for_traj(CONFIG, Work,...
%                        argoTrajWork,argoWork,argo)
%
% DESCRIPTION
% DOXY_get_profile_temp_for_traj gets the first correct value of TEMP in
% the main profile to set it up in the trajectory information, for the
% INAIR correction computation. The correct value is defined as the first
% value with allowed QC and in the shallower between Pcutoff and a deepest
% allowed pressure. The allowed QC and the maximum pressure are defined in
% DOXY_config.m.
%
% INPUT
%     CONFIG (struct)       Configuration structure with data path,
%                           operator choices, ...
%                             CONFIG =
%                        DataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/data_coriolis/INAIR/'
%                    NCEPDataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04\Data_NCEP'
%                    saveDataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/data/'
%                    savePlotDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/plots/'
%                              ...
%                              ...
%                    inAirFormat: 0
%                        inAirMC: 1090
%                              ...
%                              ...
%     Work (structure)       Doxy correction working structure, issued and
%                            computed from argo float data.
%                            Example:
%                             Work =
%                                   readme: [1x1 struct]
%                                     unit: 'mumol/kg'
%                                      wmo: 1901205
%                                 DOXY_RAW: [85x120 single]
%                                  CAPTEUR: {'Optode'}
%                                  DOXY_QC: [85x120 char]
%
%    argoTrajWork (structure)    float trajectory intermediate working
%                                 data
% 
%    argowork (structure)     Float working structure, issued and computed
%                             from argo float data.
%                             Example:
%                             argoWork = 
%                               pres_adjusted: [1x1 struct]
%                               temp_adjusted: [1x1 struct]
%                               psal_adjusted: [1x1 struct]
%
%   argo (structure)        Argo float structure (get directly from data).
%                         Example :
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
% OUTPUT
%     argoTrajWork (structure)    float trajectory intermediate working
%                                 data completed
%     goProg (double)             0\1. If this index is 0, LOCODOX stop the
%                                 treatment in progress
% CALL :
%
% SEE ALSO
%   DOXY_argoTraj_read, DOXY_corr_main

% HISTORY
%   $created: 19/04/2018 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:

function [argoTrajWork, goProg] = DOXY_get_profile_temp_for_traj(CONFIG, Work,...
    argoTrajWork,argoWork,argo)

% =========================================================================
%% Initialisation
% =========================================================================
goProg = 1;

maxPres = CONFIG.inAirMaxPresForTS;

temp_adjusted.data = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
temp_adjusted.cycle_number = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
temp_adjusted.pres = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
noDataShallow = false(1,size(argoTrajWork.profile.temp_adjusted.data,1));
temp_adjusted.from = 'main profile';

% =========================================================================
%% Get temp from main profile
% TEMP is get as the first correct value (QC ok) in the main profile
% between Pcutoff and a maximum depth defined in the configuration file
% (ex: 20m).
% =========================================================================
inlinePres = 1;
while inlinePres == 1
     logId = fopen(fullfile(Work.dirLog,sprintf('DOXY_get_nearsurf_PTS_for_traj_TEMP_%d.log',Work.wmo)),'w');
     fprintf(logId,'WMO; PRES; TEMP; TEMP_QC; WARN\n');
    
    % ---------------------------------------------------------------------
    % find shallowed correct value (QC allowed (CONFIG.QC_T), max pres
    % allowed (CONFIG.inAirMaxPresForTS).
    % If no value find, warn user and offer to change the threshold pressure
    % masPres.
    % ---------------------------------------------------------------------
    for icycle = 1:length(argoTrajWork.ppox_doxy_adjusted.data)
        iCyc = argo.cycle_number.data == argoTrajWork.cycle_number.data(icycle);
        if any(iCyc)
            isokT = ~isnan(argoWork.temp_adjusted.data(iCyc,:));
            isokP = argoWork.pres_adjusted.data(iCyc,:) < maxPres;
            tmpQC = argoWork.temp_adjusted_qc.data(iCyc,:);
            tmpQC(tmpQC == ' ') = '9';
            isokQC = ismember(str2num(tmpQC'),CONFIG.QC_T)';
            isok = isokP & isokT & isokQC;
            iShallower = find(isok,1,'first');
            
            if ~isempty(iShallower)
                % valid data (not NaN, allowed QC, pres over maxPres)
                % in argoWork (inair correction => Near-Surface
                % profile, bearing DOXY)
                pres = argoWork.pres_adjusted.data(iCyc,iShallower);
                temp = argoWork.temp_adjusted.data(iCyc,iShallower);
                temp_qc = argoWork.temp_adjusted_qc.data(iCyc,iShallower);
                fprintf(logId,sprintf('%d, %2.2f, %2.2f %d %s\n',Work.wmo,pres,temp,temp_qc,'OK'));
            else
                pres = NaN;
                temp = NaN;
                temp_qc = '9';
                fprintf(logId,sprintf('%d, %2.2f, %2.2f %d %s\n',Work.wmo,pres,temp,temp_qc,'KO'));
                noDataShallow(iCyc) = true;
            end
            
            temp_adjusted.cycle_number(icycle) = argoTrajWork.cycle_number.data(icycle);
            temp_adjusted.data(icycle) = temp;
            temp_adjusted.pres(icycle) = pres;
            %--EB
            % WARNING : if data all nan, LOCODOX should take the first shallower
            % data in the profile not bearing doxy (near-surface "core"), then
            % in the profile bearing doxy ('secondary' or 'primary'), then in
            % the profile core.
            %--EB
        end
    end
    
    fclose(logId);
    if any(~isnan(temp_adjusted.data))       
        % -----------------------------------------------------------------
        % if more than 30% of the cycles has no valid data (QC KO or temp at NaN or no pres between
        % surface and maxPres), ask user to modify maxPres in order to look
        % for valid TEMP below maxPres.
        % -----------------------------------------------------------------
        badDataPerc = (sum(noDataShallow)*100)/length(noDataShallow);
        if badDataPerc < 30
            inlinePres = 0;
            argoTrajWork.profile.temp_adjusted.data = temp_adjusted.data;
            argoTrajWork.profile.temp_adjusted.cycle_number = temp_adjusted.cycle_number;
            argoTrajWork.profile.temp_adjusted.pres = temp_adjusted.pres;
            argoTrajWork.profile.temp_adjusted.from = 'Main profile';

            clear temp_adjusted;            
        elseif badDataPerc > 30 && badDataPerc < 100
            prompt = {'More than 30% of the cycles has no valid TEMP data in the main profil',...
                sprintf('over maxPRES = %2.2f.',maxPres),...
                'For these cycles, the TEMP values has been be set to NaN.',...
                'If you want to look deeper, you had to change the maxPres',...
                'directly or in the configuration file.',...
                'Do you want to change directly the maxPres for all cycles ?'};
            
            input = questdlg(prompt,'DOXY_argoTraj_read : NOT ENOUGH TEMP',...
                'NO','CHANGE','CHANGE');
            if strcmp(input,'CHANGE')
                maxPres = inputdlg('change the maxPres value:','MODIFY MAXPRES',1,{num2str(maxPres)});
                maxPres = str2num(maxPres{1});
            else
                inlinePres = 0;
                goProg = 0;
            end
        end
    else
        % no valid data over all the cycles : change maxPres or change
        % the temp recovery strategy
        % WARNING : if data all nan, LOCODOX should take the first shallower
        % data in the profile not bearing doxy (near-surface "core"), then
        % in the profile bearing doxy ('secondary' or 'primary'), then in
        % the profile core.
        % NOT IMPLEMENTED
        prompt = {sprintf('No valid data over all the cycles, below %dm',maxPres),...
            '',...
            '!!!change maxPres or change the temp recovery strategy!!!',...
            'if data all nan, LOCODOX should take the first shallower',...
            'data in the profile not bearing doxy (near-surface "core"),',...
            'then in the profile bearing doxy (''secondary'' or ''primary''),',...
            'then in the profile core.',...
            '!!!NOT YET IMPLEMENTED!!!'};
        warndlg(prompt,'DOXY_argoTraj_read : NO TEMP','modal');
        goProg = 0;
        inlinePres = 0;
    end
    inlinePres = 0;
end

% =========================================================================
% repeat the temp profile data as in argoTrajWork.temp_adjusted.data
% =========================================================================
argoTrajWork.profile.cellShape.temp_adjusted.data = cell(1,length(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.cellShape.temp_adjusted.pres = cell(1,length(argoTrajWork.ppox_doxy_adjusted.data));
for icycle = 1:length(argoTrajWork.cycle_number.data)
    if ~isempty(argoTrajWork.temp_adjusted.data{icycle})
        % compare cycle_number
        okCyc = argoTrajWork.profile.temp_adjusted.cycle_number == argoTrajWork.cycle_number.data(icycle);
        % set psal_adjusted to the primary profil value
        if any(okCyc)
            argoTrajWork.profile.cellShape.temp_adjusted.data{icycle} = repmat(argoTrajWork.profile.temp_adjusted.data(okCyc),size(argoTrajWork.temp_adjusted.data{icycle}));
            argoTrajWork.profile.cellShape.temp_adjusted.pres{icycle} = repmat(argoTrajWork.profile.temp_adjusted.pres(okCyc),size(argoTrajWork.temp_adjusted.data{icycle}));
        end
    end
end