% DOXY_get_PTS select the PTS fields used to correct the DOXY data in argo
%
% SYNTAX
% [Sout,Work] = DOXY_get_PTS(argo,Work,LogDir)
%
% DESCRIPTION
% DOXY_get_PTS get the PTS fields used to correct the DOXY data in argo
% float. The PTS data are selected:
%       - by default: data adjusted if exist, raw data if not
%       - if it's the user's choice: raw data (even if the adjusted data
%                                              exist)
%       - the only exception (Whatever the user choice): If VSS is Near-Surface sampling (unpumped),
%         get the psal instead of the psal_adjusted, because unpumped data get artificially a QC = 4
%         (so data = NaN), not justified for DOXY calculation
%
% INPUT
%     argo (structure)       float structure
%                           Example :
%                               argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%
%                               argo.pres =
%                                   name: 'PRES'
%                                    dim: {'N_PROF'  'N_LEVELS'}
%                                   data: [85x120 single]
%                              long_name: 'SEA PRESSURE'
%                          standard_name: 'sea_water_pressure'
%                             FillValue_: 99999
%                                  units: 'decibar'
%                              valid_min: 0
%                              valid_max: 12000
%                               C_format: '%7.1f'
%                         FORTRAN_format: 'F7.1'
%                             resolution: 0.1000
%                                   axis: 'Z'
%             coordinate_reference_frame: 'urn:ogc:crs:EPSG::5113'
%                                   type: 5
%
%   Work (struct)        float working structure, issued and computed from
%                        the argo data and the reference data (WOA or REF).
%                        Carries plot informations among other things.
%                        Example:
%                         Work =
%
%                                  readme: [1x1 struct]
%                                    unit: 'mumol/kg'
%                              R2treshold: 0.8000
%                                     wmo: 1901205
%                                  sensor: 'Optode'
%                                makePlot: 1
%                                savePlot: 1
%
%   LogDir             .....
%
% OUTPUT
%     Sout (structure)      float structure on the Sin model, with the PTS
%                           adjusted field
%                           Example :
%                             Sout =
%                                 pres_adjusted: [1x1 struct]
%                                 temp_adjusted: [1x1 struct]
%                                 psal_adjusted: [1x1 struct]
%
%     Work (structure)       float structure Work
%
% CALL : 
%   DOXY_PLOT_PTS_QCD
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 25/11/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%         v2   26/01/2017 Anne Piron, Emilie Brion, Altran Ouest
%                         argo - ameliorations 2017 phase 1
%         v2   08/11/2017 Emilie Brion, Altran Ouest
%                         If no QC > 2, the plot is still done, but no
%                         message dialog box is displayed.
%         v2.1 10/07/2018 Emilie Brion, Altran Ouest
%                         Add comments, clear the programm.

function [Sout,Work] = DOXY_get_PTS(argo,Work,LogDir)

% =========================================================================
%% Initialisation
% =========================================================================
myDim = argo.pres_adjusted.dim;
okDim = find(strcmp(myDim,'N_PROF'));
Work.betterUseDoxy = 0;

% =========================================================================
%% QC selection: replaces the data with NaN for all unwanted QCs
% =========================================================================
% All possible QCs
QC_po = 1:1:9;

% -------------------------------------------------------------------------
% 1. PRES
% -------------------------------------------------------------------------
QCnan = NaN(1,length(QC_po)-length(Work.QC_P));
cmpt = 0;
for iq = 1:length(QC_po)
    i = find(QC_po(iq) == Work.QC_P, 1);
    if isempty(i)
        cmpt = cmpt + 1;
        QCnan(cmpt) = QC_po(iq);
    end
    clear i
end
for iq = 1:length(QCnan)
    if Work.Wdata.DM_pres
        argo.pres_adjusted.data(argo.pres_adjusted_qc.data==num2str(QCnan(iq))) = NaN;
    else
        argo.pres.data(argo.pres_qc.data==num2str(QCnan(iq))) = NaN;
    end
end

% -------------------------------------------------------------------------
% 2. TEMP
% -------------------------------------------------------------------------
QCnan = NaN(1,length(QC_po)-length(Work.QC_T));
cmpt = 0;
for iq = 1:length(QC_po)
    i = find(QC_po(iq) == Work.QC_T, 1);
    if isempty(i)
        cmpt = cmpt + 1;
        QCnan(cmpt) = QC_po(iq);
    end
    clear i
end
for iq = 1:length(QCnan)
    if Work.Wdata.DM_temp
        argo.temp_adjusted.data(argo.temp_adjusted_qc.data==num2str(QCnan(iq))) = NaN;
    else
        argo.temp.data(argo.temp_qc.data==num2str(QCnan(iq))) = NaN;
    end
end

% -------------------------------------------------------------------------
% 3. PSAL
% -------------------------------------------------------------------------
QCnan = NaN(1,length(QC_po)-length(Work.QC_S));
cmpt = 0;
for iq = 1:length(QC_po)
    i = find(QC_po(iq) == Work.QC_S, 1);
    if isempty(i)
        cmpt = cmpt + 1;
        QCnan(cmpt) = QC_po(iq);
    end
    clear i
end
for iq = 1:length(QCnan)
    if Work.Wdata.DM_psal
        argo.psal_adjusted.data(argo.psal_adjusted_qc.data==num2str(QCnan(iq))) = NaN;
    else
        argo.psal.data(argo.psal_qc.data==num2str(QCnan(iq))) = NaN;
    end
end

% =========================================================================
%% Manage the combo QC and Data Mode
% Prepare the log file to retrieve a posteriori the choice made
% =========================================================================
% Vector (logical) getR = 1 if Core File 'R' (0 if not,i.e. Core File 'D')
% by default: adjusted data (=> 0)
getR_pres = false(1,size(argo.pres_adjusted.data,okDim));
getR_temp = false(1,size(argo.pres_adjusted.data,okDim));
getR_psal = false(1,size(argo.pres_adjusted.data,okDim));

cmptR = 0; cmptD = 0;

for ip = 1:length(getR_pres)    
    % Indices of the cycle numbers for Core Files R and D
    if strcmp(argo.data_mode_CF.data(ip),'R')
        cmptR = cmptR + 1; indR(cmptR) = ip;
    elseif strcmp(argo.data_mode_CF.data(ip),'D')
        cmptD = cmptD + 1; indD(cmptD) = ip;
    end
    
    % Case: 'D' File does not exist => LOCODOX take the 'R' File data
    if strcmp(argo.data_mode_CF.data(ip),'R')
        getR_pres(ip) = true;
        getR_temp(ip) = true;
        getR_psal(ip) = true;
    end
end

% User's choice: force the use of the 'Real Time' fields even if the
% 'Delayed Mode' fields exist (DM_fields set to 0)
if Work.Wdata.DM_pres == 0
    getR_pres(:) = true;
end
if Work.Wdata.DM_temp == 0
    getR_temp(:) = true;
end
if Work.Wdata.DM_psal == 0
    getR_psal(:) = true;
end

% cycle numbers of the Core Files
if exist('indR','var')
    cycn_CF_R = argo.cycle_number.data(indR);
end
if exist('indD','var')
    cycn_CF_D = argo.cycle_number.data(indD);
end

% =========================================================================
%% Filled the Sout "adjusted" field.
% If VSS is Near-Surface sampling (unpumped), get the psal instead of
% the psal_adjusted, because unpumped data get artificially a QC = 4
% (so data = NaN), not justified for DOXY calculation
% =========================================================================
VarIn = {'pres','temp','psal'};
for v = 1:length(VarIn)
    if strcmp(VarIn{v},'psal') && ~isempty(strfind(argo.VSS,'Near')) 
        Sout.([VarIn{v} '_adjusted']) = argo.([VarIn{v}]);
        Sout.([VarIn{v} '_adjusted_qc']) = argo.([VarIn{v} '_qc']);
        Sout.([VarIn{v} '_adjusted']).name = upper([VarIn{v} '_adjusted']);
        getR_psal_log = true(1,length(getR_psal));
    else
        if strcmp(VarIn{v},'pres')
            getR_pres_log = getR_pres; getR_log = getR_pres_log;
        elseif strcmp(VarIn{v},'temp')
            getR_temp_log = getR_temp; getR_log = getR_temp_log;
        elseif strcmp(VarIn{v},'psal') && isempty(strfind(argo.VSS,'Near'))
            getR_psal_log = getR_psal; getR_log = getR_psal_log;
        end
        Sout.([VarIn{v} '_adjusted']) = argo.([VarIn{v} '_adjusted']);
        Sout.([VarIn{v} '_adjusted_qc']) = argo.([VarIn{v} '_adjusted_qc']);
        Sout.([VarIn{v} '_adjusted']).data(getR_log,:) = argo.(VarIn{v}).data(getR_log,:);
        Sout.([VarIn{v} '_adjusted_qc']).data(getR_log,:) = argo.([VarIn{v} '_qc']).data(getR_log,:);
        Sout.([VarIn{v} '_adjusted']).name = upper([VarIn{v} '_adjusted']);
    end
    clear getR_log
end

% =========================================================================
%% Writing to a Log File: which Data used to correct the DOXY data
% =========================================================================
if Work.checkPTS == 1
%     Work.checkPTS = 0;
    if exist('logId','var')
        clear logId
    end
    logId = fopen(fullfile(LogDir,sprintf('DOXY_PTS_data_used_for_doxycorr_%d.log',Work.wmo)),'w');
    fprintf(logId,sprintf('%d float ; PTS of the reference profile (%s) selected by LOCODOX for DOXY correction \nCore Files: \n',Work.wmo,argo.VSS));
    if exist('indR','var') && exist('indD','var') %!!!!!!!!!!!!!!!!!!! PAS TESTE !!!!!!!!!!!!!!!!!!!
        fprintf(logId,' - Delayed Mode: cycle numbers ');
        for i = 1:length(cycn_CF_D)
            fprintf(logId,sprintf(' %d ;',cycn_CF_D(i)));
            if i == length(cycn_CF_D); fprintf(logId,'\n'); end
        end
        fprintf(logId,' - Real Time: cycle numbers ');
        for i = 1:length(cycn_CF_R)
            fprintf(logId,sprintf(' %d ;',cycn_CF_R(i)));
            if i == length(cycn_CF_R); fprintf(logId,'\n'); end
        end
    elseif exist('indR','var') && ~exist('indD','var') %!!!!!!!!!!!!!!!!!!! PAS TESTE !!!!!!!!!!!!!!!!!!!
        fprintf(logId,' - Delayed Mode: no cycle numbers\n - Real Time: cycle numbers ');
        for i = 1:length(cycn_CF_R)
            fprintf(logId,sprintf(' %d ;',cycn_CF_R(i)));
            if i == length(cycn_CF_R); fprintf(logId,'\n'); end
        end
    elseif ~exist('indR','var') && exist('indD','var')
        fprintf(logId,' - Delayed Mode: cycle numbers');
        for i = 1:length(cycn_CF_D)
            fprintf(logId,sprintf(' %d ;',cycn_CF_D(i)));
            if i == length(cycn_CF_D); fprintf(logId,'\n'); end
        end
        fprintf(logId,' - Real Time: no cycle numbers \n');
    end
    fprintf(logId,'=========================================================================\n');
    VarName = {'Pressure','Temperature','Salinity'};
    for ic = 1:length(VarIn)
        if strcmp(VarIn{ic},'pres')
            getR_log = getR_pres_log;
        elseif strcmp(VarIn{ic},'temp')
            getR_log = getR_temp_log;
        elseif strcmp(VarIn{ic},'psal')
            getR_log = getR_psal_log;
        end
        if length(unique(getR_log)) == 1
            if unique(getR_log) == 0
                fprintf(logId,sprintf('%s data: \n - %s_ADJUSTED fields: all cycle numbers \n',VarName{ic},upper(VarIn{ic})));
            elseif unique(getR_log) == 1
                fprintf(logId,sprintf('%s data: \n - %s fields: all cycle numbers \n',VarName{ic},upper(VarIn{ic})));
            end
        else
            indR = find(getR_log==1);
            fprintf(logId,sprintf('%s data: \n - %s fields for cycle numbers:',VarName{ic},upper(VarIn{ic})));
            for i = 1:length(indR)
                fprintf(logId,sprintf(' %d ;',argo.cycle_number.data(indR(i))));
                if i == length(indR); fprintf(logId,'\n'); end
            end; clear indR
            indD = find(getR_log==0);
            fprintf(logId,sprintf(' - %s_ADJUSTED fields for cycle numbers:',upper(VarIn{ic})));
            for i = 1:length(indD)
                fprintf(logId,sprintf(' %d ;',argo.cycle_number.data(indD(i))));
                if i == length(indD); fprintf(logId,'\n'); end
            end; clear indD
        end
        fprintf(logId,'=========================================================================\n');
        clear getR_log
    end
    fclose(logId);
end

% =========================================================================
%% Check the salinity data quality: Plot QC flag if selected data (PTS)
% are not of good quality (QC > 2)
% =========================================================================
if Work.checkPTS == 1
    Work.checkPTS = 0;
    % Control Figure:
    if Work.makePlot
        hFig = figure('Name',sprintf('PTS data QC - %d',Work.wmo),'NumberTitle','off',...
            'unit','normalized','Outerposition',[0.70 0.26 0.30 0.74]);
        DOXY_PLOT_PTS_QC(hFig,...
            argo.cycle_number.data,...
            Sout.pres_adjusted.data,...
            Sout.pres_adjusted_qc.data,...
            Sout.temp_adjusted_qc.data,...
            Sout.psal_adjusted_qc.data,...
            Work);
    end
    
    % Warning Message:
    if ~isempty(find(str2numQC(Sout.pres_adjusted_qc.data) > 2, 1))
        warn_P = '- PRESSURE';
    else
        warn_P = ' ';
    end
    if ~isempty(find(str2numQC(Sout.temp_adjusted_qc.data) > 2, 1))
        warn_T = '- TEMPERATURE';
    else
        warn_T = ' ';
    end
    if ~isempty(find(str2numQC(Sout.psal_adjusted_qc.data) > 2, 1)) && ~strcmp(Work.whichCorr,'INAIR')
        warn_S = '- SALINITY => LOCODOX suggests using DOXY (instead of PSAT) for the oxygen correction';
        Work.betterUseDoxy = 1;
    elseif ~isempty(find(str2numQC(Sout.psal_adjusted_qc.data) > 2, 1)) && strcmp(Work.whichCorr,'INAIR')
        warn_S = '- SALINITY';
    else
        warn_S = ' ';
    end
    
    if (~isempty(find(str2numQC(Sout.pres_adjusted_qc.data) > 2, 1)) ||...
            ~isempty(find(str2numQC(Sout.temp_adjusted_qc.data) > 2, 1)) ||...
            ~isempty(find(str2numQC(Sout.psal_adjusted_qc.data) > 2, 1))) && ~isfield(Work,'is_not_inair_corr')
        rep = warndlg({['WARNING : Reference profile (' argo.VSS ') => some PTS data are not of good quality (QC > 2). It concerns the following fields:']...
            warn_P warn_T warn_S...
            '     ' 'Press OK to continue' },...
            '!! WARNING QC !!');
        
        log_summ = Work.log_summ;
        message=sprintf('WARNING : Reference profile (%s) => some PTS data are not of good quality (QC > 2).',argo.VSS );%
          fprintf(log_summ,'----------------------------------------------------\n');
          fprintf(log_summ,[message '\n']); clear message
          fprintf(log_summ,'It concerns the following fields: \n');
          fprintf(log_summ,'%s\n',warn_P);
          fprintf(log_summ,'%s\n',warn_T);
          fprintf(log_summ,'%s\n',warn_S);
          fprintf(log_summ,'----------------------------------------------------\n');  
          fprintf(log_summ,'\n'); 
        waitfor(rep);
    end
end