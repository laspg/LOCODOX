% DOXY_get_primary_PTS_for_traj gets the TEMP and PSAL data in the primary
% profile, for inAir correction.
%
% SYNTAX
% [argoTrajWork, goProg] = DOXY_get_primary_PTS_for_traj(CONFIG, Work,...
%                          argoTrajWork,argo1Struct,argo2Struct,argo3Struct,...
%                          argo4Struct)
% DESCRIPTION
% DOXY_get_primary_PTS_for_traj gets the first correct value of PSAL and
% TEMP in the primary profile to set it up in the trajectory information,
% for the INAIR correction computation. The correct value is defined as the
% first value with allowed QC and in the shallower between Pcutoff and a
% deepest allowed pressure. The allowed QC and the maximum pressure are
% defined in DOXY_config.m.
% The Pressure retained is that associated with the temperature.
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
%     argoTrajWork (structure)    float trajectory intermediate working
%                                 data
%
%     argo1Struct (struct)    Three data structures organised with the main profile
%     argo2Struct (struct)    for working is in argo1Struct (vertical
%     argo3Struct (struct)    sampling shceme of interest for the kind of
%     argo4Struct (struct)    correction), and the others are in the second ones.
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
%                         with:
%                         argo1Struct.argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                                         ...
%
%                         argo1Struct.argoWork =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%                              the         sat: [1x1 struct]
%                                      psat: [1x1 struct]
%
%                          argo1Struct.Work =
%                                    readme: [1x1 struct]
%                                      unit: 'mumol/kg'
%                                       wmo: 1901205
%                                    sensor: 'Optode'
%                                 whichCorr: 'WOA'
%                                  DOXY_RAW: [165x117 single]
%                                     timar: [165x20 char]
%                                     datat: [165x1 double]
%                                   DENSITY: [165x117 single]
%                                     DEPTH: [165x117 double]
%                                   CAPTEUR: 'Optode'
%                                inlinePres   DOXY_QC: [165x117 char]
%                                      DENS: [165x117 single]
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
%   $created: 01/07/2018 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:

function [argoTrajWork, goProg] = DOXY_get_primary_PTS_for_traj(CONFIG, Work,...
    argoTrajWork,argo1Struct,argo2Struct,argo3Struct, argo4Struct) 


% =========================================================================
%% Initialisation
% =========================================================================
goProg = 1;
maxPres = CONFIG.inAirMaxPresForTS;

% Find the primary sampling profile
for is = 1:4
    vss = eval(sprintf('argo%dStruct.VSS',is));
%     if contains(vss,'Primary')
    if ~isempty(strfind(vss,'Primary'))
        break
    end
end
primaryStruct = eval(sprintf('argo%dStruct',is));
isA = primaryStruct.argo.direction.data == 'A';

argoTrajWork.profile.psal_adjusted.data = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.psal_adjusted.cycle_number = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.psal_adjusted.pres = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.temp_adjusted.data = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.temp_adjusted.cycle_number = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.temp_adjusted.pres = NaN(size(argoTrajWork.ppox_doxy_adjusted.data));
noDataShallow.temp_adjusted = false(1,size(argoTrajWork.profile.temp_adjusted.data,2));
noDataShallow.psal_adjusted = false(1,size(argoTrajWork.profile.psal_adjusted.data,2));
argoTrajWork.profile.psal_adjusted.from = 'primary profile';
argoTrajWork.profile.temp_adjusted.from = 'primary profile';

% =========================================================================
%% Get Psal and Temp in the primary profile
% PSAL
% the CTD pump is switched off in the near-surface and surface phases.
% Moreover, no PSAL could be measured in the "in-air" phase.
% Therefore, PSAL is get as the first correct value (QC ok) in the primary
% profile between Pcutoff and a maximum depth defined in the configuration
% file (ex: 20m).
%
% TEMP
% TEMP is get as the first correct value (QC ok) in the primary profile
% between Pcutoff and a maximum depth defined in the configuration file
% (ex: 20m).
% =========================================================================

inlinePres = 1;
while inlinePres == 1
    % logfile
    logIdS = fopen(fullfile(Work.dirLog,sprintf('DOXY_get_primary_PTS_for_traj_PSAL_%d.log',Work.wmo)),'w');
    fprintf(logIdS,'WMO; PRES; PSAL; PSAL_QC; WARN\n');
    logIdT = fopen(fullfile(Work.dirLog,sprintf('DOXY_get_primary_PTS_for_traj_TEMP_%d.log',Work.wmo)),'w');
    fprintf(logIdT,'WMO; PRES; TEMP; TEMP_QC; WARN\n');
    
    % ---------------------------------------------------------------------
    % find shallowed correct value (QC allowed (CONFIG.QC_S), max pres
    % allowed (CONFIG.inAirMaxPresForTS).
    % If no value find, warn user and offer to change the threshold pressure
    % masPres.
    % ---------------------------------------------------------------------
    for icycle = 1:length(argoTrajWork.ppox_doxy_adjusted.data)
        iCyc = primaryStruct.argo.cycle_number.data == argoTrajWork.cycle_number.data(icycle) & isA;
        if any(iCyc)
            % PSAL
            isokS = ~isnan(primaryStruct.argoWork.psal_adjusted.data(iCyc,:));
            isokP = primaryStruct.argoWork.pres_adjusted.data(iCyc,:) < maxPres;
            tmpQC = primaryStruct.argoWork.psal_adjusted_qc.data(iCyc,:);
            tmpQC(tmpQC == ' ') = '9';
            isokQC = ismember(str2num(tmpQC'),CONFIG.QC_S)';
            isok = isokP & isokS & isokQC;
            iShallower = find(isok,1,'first');
            
            if ~isempty(iShallower)
                % valid data (not NaN, allowed QC, pres over maxPres)
                % in argoWork (inair correction => Near-Surface
                % profile, bearing DOXY)
                pres = primaryStruct.argoWork.pres_adjusted.data(iCyc,iShallower);
                psal = primaryStruct.argoWork.psal_adjusted.data(iCyc,iShallower);
                psal_qc = primaryStruct.argoWork.psal_adjusted_qc.data(iCyc,iShallower);
                fprintf(logIdS,sprintf('%d, %2.2f, %2.2f %d %s\n',Work.wmo,pres,psal,str2num(psal_qc),'OK'));
            else
                pres = NaN;
                psal = NaN;
                psal_qc = '9';
                fprintf(logIdS,sprintf('%d, %2.2f, %2.2f %d %s\n',Work.wmo,pres,psal,psal_qc,'KO'));
                noDataShallow.psal_adjusted(iCyc) = true;
            end
            argoTrajWork.profile.psal_adjusted.cycle_number(icycle) = argoTrajWork.cycle_number.data(icycle);
            argoTrajWork.profile.psal_adjusted.data(icycle) = psal;
            argoTrajWork.profile.psal_adjusted.pres(icycle) = pres;
            
            % TEMP
            isokT = ~isnan(primaryStruct.argoWork.temp_adjusted.data(iCyc,:));
            isokP = primaryStruct.argoWork.pres_adjusted.data(iCyc,:) < maxPres;
            tmpQC = primaryStruct.argoWork.temp_adjusted_qc.data(iCyc,:);
            tmpQC(tmpQC == ' ') = '9';
            isokQC = ismember(str2num(tmpQC'),CONFIG.QC_T)';
            isok = isokP & isokT & isokQC;
            iShallower = find(isok,1,'first');
            
            if ~isempty(iShallower)
                % valid data (not NaN, allowed QC, pres over maxPres)
                % in argoWork (inair correction => Near-Surface
                % profile, bearing DOXY)
                pres = primaryStruct.argoWork.pres_adjusted.data(iCyc,iShallower);
                temp = primaryStruct.argoWork.temp_adjusted.data(iCyc,iShallower);
                temp_qc = primaryStruct.argoWork.temp_adjusted_qc.data(iCyc,iShallower);
                fprintf(logIdT,sprintf('%d, %2.2f, %2.2f %d %s\n',Work.wmo,pres,temp,str2num(temp_qc),'OK'));
            else
                pres = NaN;
                temp = NaN;
                temp_qc = '9';
                fprintf(logIdS,sprintf('%d, %2.2f, %2.2f %d %s\n',Work.wmo,pres,temp,temp_qc,'KO'));
                noDataShallow.temp_adjusted(iCyc) = true;
            end
            argoTrajWork.profile.temp_adjusted.cycle_number(icycle) = argoTrajWork.cycle_number.data(icycle);
            argoTrajWork.profile.temp_adjusted.data(icycle) = temp;
            argoTrajWork.profile.temp_adjusted.pres(icycle) = pres;
        end
    end
    
    % ---------------------------------------------------------------------
    % if more than 30% of the cycles has no valid data (QC KO or temp at
    % NaN or no pres between surface and maxPres), ask user to modify
    % maxPres in order to look for valid TEMP below maxPres.
    % ---------------------------------------------------------------------    
    var = {'temp','psal'};
    for iv = 1:2
        myvar = var{iv};
        badDataPerc = (sum(noDataShallow.([myvar '_adjusted']))*100)/length(noDataShallow.([myvar '_adjusted']));
        if badDataPerc > 30 && badDataPerc < 100
            prompt = {'More than 30% of the cycles has no valid TEMP data in the main profil',...
                sprintf('over maxPRES = %2.2f.',maxPres),...
                sprintf('For these cycles, the %s values has been be set to NaN.',upper(myvar)),...
                'If you want to look deeper, you had to change the maxPres',...
                'directly or in the configuration file.',...
                'Do you want to change directly the maxPres for all cycles ?'};
            
            input = questdlg(prompt,sprintf('DOXY_argoTraj_read : NOT ENOUGH %s',upper(myvar)),...
                'NO','CHANGE','CHANGE');
            if strcmp(input,'CHANGE')
                maxPres = inputdlg('change the maxPres value:','MODIFY MAXPRES',1,{num2str(maxPres)});
                maxPres = str2num(maxPres{1});
            else
                inlinePres = 0;
                goProg = 0;
            end
        elseif badDataPerc == 100
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
            warndlg(prompt,sprintf('DOXY_argoTraj_read : NO %s',upper(myvar)),'modal');
            goProg = 0;
            inlinePres = 0;
        end
    end
    inlinePres = 0;
end

% =========================================================================
% repeat the psal profile data as in argoTrajWork.psal_adjusted.data
% repeat the temp profile data as in argoTrajWork.temp_adjusted.data
% =========================================================================
% TEMP
argoTrajWork.profile.cellShape.temp_adjusted.data = cell(1,length(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.cellShape.temp_adjusted.pres = cell(1,length(argoTrajWork.ppox_doxy_adjusted.data));
for icycle = 1:length(argoTrajWork.cycle_number.data)
     if ~isempty(argoTrajWork.temp_adjusted.data{icycle})
        % compare cycle_number
        okCyc = argoTrajWork.profile.temp_adjusted.cycle_number == argoTrajWork.cycle_number.data(icycle);
        % set temp_adjusted to the primary profil value
        if any(okCyc)
            argoTrajWork.profile.cellShape.temp_adjusted.data{icycle} = repmat(argoTrajWork.profile.temp_adjusted.data(okCyc),size(argoTrajWork.temp_adjusted.data{icycle}));
            argoTrajWork.profile.cellShape.temp_adjusted.pres{icycle} = repmat(argoTrajWork.profile.temp_adjusted.pres(okCyc),size(argoTrajWork.temp_adjusted.data{icycle}));
        end
     end
end
% PSAL
argoTrajWork.profile.cellShape.psal_adjusted.data = cell(1,length(argoTrajWork.ppox_doxy_adjusted.data));
argoTrajWork.profile.cellShape.psal_adjusted.pres = cell(1,length(argoTrajWork.ppox_doxy_adjusted.data));
for icycle = 1:length(argoTrajWork.psal_adjusted.data)
    if ~isempty(argoTrajWork.psal_adjusted.data{icycle})
        % compare cycle_number
        okCyc = argoTrajWork.profile.psal_adjusted.cycle_number == argoTrajWork.cycle_number.data(icycle);
        % set psal_adjusted to the primary profil value
        if any(okCyc)
            argoTrajWork.profile.cellShape.psal_adjusted.data{icycle} = repmat(argoTrajWork.profile.psal_adjusted.data(okCyc),size(argoTrajWork.psal_adjusted.data{icycle}));
            argoTrajWork.profile.cellShape.psal_adjusted.pres{icycle} = repmat(argoTrajWork.profile.psal_adjusted.pres(okCyc),size(argoTrajWork.psal_adjusted.data{icycle}));
        end
    end
end
