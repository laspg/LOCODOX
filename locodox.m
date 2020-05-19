% locodox apply correction on the Argo Oxygen data
%
% SYNTAX
% [] = locodox
%
% DESCRIPTION
% locodox apply correction on the Argo Oxygen data, using the ATLN
% method.
%
% INPUT
%   config_file     (Optional) the full path of the configuration file.
%
% OUTPUT
%   Adjusted file written
%
% CALL
%   DOXY_NCR_WOA, DOXY_argo_read, DOXY_argo_prep, DOXY_MAP, DOXY_PLOT_QC,
%   DOXY_interp_WOA_main, grep -ni doxy_drift *.m, DOXY_corr_compute,
%   DOXY_corr_apply_main,DOXY_interp_REF_main,DOXY_update_fields,
%
% SEE ALSO
%

% HISTORY
%   $created: //2009 $author: Mathieu Le Steun, Thomas Bouinot, LPO
%   $Revision: version $Date: $author:
%       v1.2  /
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       ergonomics
%       v3 25/05/2016   Emilie Brion, ALTRAN OUEST
%                       argo 3.1 adjustement, Vertical Sampling Scheme
%                       management
%       v4 26/01/2017   Anne Piron, Emilie Brion, Altran Ouest
%                       argo - ameliorations 2017 phase 1
%       v4 21/02/2019   Marine GALLIAN, Altran Ouest
%                       Take into account when there is no PPOX data, the 
%                       IN AIR correction is not possible
%                       Add pressure effect correction and display
%                       corrected data 
%                       Write configuration and choices of the user in the 
%                       summary file
%      v5 07/02/2020    Access to M_MAP added by T. Reynaud
%      V5 12/03/2020    Bug corrected by Thierry Reynaud for NC/MAT directories for REF case
%                       (dirout)
%      v5 05.04.2020    Modifications by T. Reynaud for plotting and
%                       generating corrected oxy profiles figures without   
%                       a REF profile
%      v6 24.04.2020    DOXY_corr_main.m renamed as locodox.m    


function [] = locodox(config_prog)

close all
% =========================================================================
%% *Initialisation*
% =========================================================================
goProg = 1;

% Read the configuration file, addpath and load files
if nargin == 0
    config_prog = '/Users/thierry_reynaud/IFREMER/MATLAB/LOCODOX/LOCODOX3.4/locodox_config';
end
[config_dir,config_prog] = fileparts(config_prog);
cd(config_dir);
CONFIG = feval(config_prog);

% Prepare results directories 
CONFIG.saveDataDir=fullfile(CONFIG.resultsDir,'data',filesep);
CONFIG.savePlotDir = fullfile(CONFIG.resultsDir,'plots',filesep);
CONFIG.saveLogDir = fullfile(CONFIG.resultsDir,'logs',filesep);
% CONFIG.NCEPDataDir ==> move in DOXY_config.m by TR 09.04.2020
% CONFIG.NCEPDataDir = fullfile(CONFIG.LocodoxMainDir,'data_input',filesep);

% Initialise Work Structure
Work = load('readme');
Work.R2threshold = CONFIG.R2threshold;
Work.adjusted_error_abs = CONFIG.adjusted_error_abs;
Work.adjusted_error_rel = CONFIG.adjusted_error_rel;

% Units
% -------------------------------------------------
% Units shortname = {'mumol/kg','mumol/L','mL/L'};
% Units longname = {'micromole per kilogram','micromole per liter','milliliter per liter'};
Work.unit = 'mumol/kg';
Work.refUnit = CONFIG.refUnit;

% General
Work.Wdata.DM_psal = CONFIG.DM_psal;
Work.Wdata.DM_temp = CONFIG.DM_temp;
Work.Wdata.DM_pres = CONFIG.DM_pres;
Work.presEff = CONFIG.presEff;
Work.QC_O = CONFIG.QC_O;
Work.QC_P = CONFIG.QC_P;
Work.QC_T = CONFIG.QC_T;
Work.QC_S = CONFIG.QC_S;
Work.makePlot = 1;
Work.savePlot = CONFIG.savePlot;
Work.formattype = CONFIG.formattype;
Work.drift_fitPolynomialDegree = CONFIG.drift_fitPolynomialDegree;
Work.drift_spec = CONFIG.drift_spec;
if CONFIG.drift_spec == 1
    Work.drift_fittype = CONFIG.drift_fittype;
end
Work.fontsize = CONFIG.fontsize;
Work.resol = sprintf('-r%d',CONFIG.resolution);
Work.history.history_software = CONFIG.history_software;
Work.history.history_reference = CONFIG.history_reference;
Work.history.history_software_release=CONFIG.history_software_release;
Work.isokC=CONFIG.isokC;
initialWork = Work;
Work.min_drift_depth=CONFIG.min_drift_depth;

% =====================================================================
%% *Read the WOA and convert in the usefull unit*
% =====================================================================
WOA.init = DOXY_WOA_read(CONFIG.WOAfile, Work.unit);
if CONFIG.presEff == 1
    [WOA.init.psatwoa] = WOA.init.psatwoa_preswoa;
    WOA.init = rmfield(WOA.init,'psatwoa_preswoa');
else
    WOA.init = rmfield(WOA.init,'psatwoa_preswoa');
end


% =====================================================================
%% *Doxy correction*
% =====================================================================
while goProg
    options.WindowStyle = 'normal';
    wmoIn = inputdlg('WMO of the float: ','ARGO Float DOXY correction',1,cellstr(''),options);
    if isempty(wmoIn)
        msgbox({'  Good Bye !  ';''},sprintf('%s',CONFIG.history_software),'custom',imread(CONFIG.logo));
        return
    end
    Work.wmo=str2double(wmoIn);
    Work.replist=fullfile(cell2mat(wmoIn),'/profiles/');
    close all

    % =================================================================
    %% Chose the correction option : WOA, REF or INAIR
    % =================================================================
    CONFIG.corrTypeDesc = {'WOA climatology','Reference In-Situ database','In-Air correction'};
    CONFIG.corrTypeShort = {'WOA','REF','INAIR'};
    nbr_tmp = listdlg('Name','METHOD FOR THE DOXY GAIN CALCULATION',...
    'PromptString','Which method do you want to use for your DOXY Gain calculation (slope/offset) ?',...
    'SelectionMode','single',...
    'Listsize',[500,150],...
    'ListString',CONFIG.corrTypeDesc);
    
    if isempty(nbr_tmp) 
        Work.whichCorr = questdlg({'No correction method has been choosen.';...
            'Do you want to restart ?'},...
            sprintf('%s: Cancel method correction choice',CONFIG.history_software'),...
            'YES','NO','NO');
        if strcmp(Work.whichCorr,'YES')
            continue
        else
            msgbox({'  Good Bye !  ';''},sprintf('%s',CONFIG.history_software),'custom',imread(CONFIG.logo));
            return
        end        
    end
    Work.whichCorr = CONFIG.corrTypeShort{nbr_tmp};
    clear nbr_tmp
    
    
    fprintf('INFO >>>>>    Apply correction from %s\n',Work.whichCorr);
    
    % Create directory for plots and logs
    if ~exist(fullfile(CONFIG.savePlotDir,Work.whichCorr,num2str(Work.wmo)),'dir')
        mkdir(fullfile(CONFIG.savePlotDir,Work.whichCorr,num2str(Work.wmo)));
    end
    Work.dirPlot = fullfile(CONFIG.savePlotDir,Work.whichCorr,num2str(Work.wmo));
    
    if ~exist(fullfile(CONFIG.saveLogDir,Work.whichCorr,num2str(Work.wmo)),'dir')
        mkdir(fullfile(CONFIG.saveLogDir,Work.whichCorr,num2str(Work.wmo)));
    end
    Work.dirLog = fullfile(CONFIG.saveLogDir,Work.whichCorr,num2str(Work.wmo));
    
    %% --------   SUMMARY FILE   --------------
    % CREATE THE DOC SUMMARIZING ALL THE OPTION CONSIDERED BY THE USER
    summary_file=fullfile(Work.dirLog,sprintf('Summary_correction_options_%d_%s.txt',Work.wmo,datestr(now,'yyyymmdd_HHMM')));
    log_summ= fopen(summary_file,'w');
    DOXY_print_config_in_logfile(log_summ,CONFIG,Work.wmo);
    Work.log_summ= log_summ;

    fprintf(log_summ,'================================================================================\n');
    fprintf(log_summ,'---------------------    PART 2 : CORRECTIONS APPLIED     -----------------------\n');
    fprintf(log_summ,'================================================================================\n');
    fprintf(log_summ,'\n');
    fprintf(log_summ,sprintf('The float correction is based on the %s method\n',Work.whichCorr));
    fprintf(log_summ,'\n');
    fprintf(log_summ,'----------------------------------------------------\n');
    fprintf(log_summ,'\n');
    
    % =====================================================================
    %% *Read the Reference In-Situ data if the correction used is REF*
    % =====================================================================
    if strcmp(Work.whichCorr,'REF')
        load(CONFIG.bddFile);
        ff=fopen(CONFIG.RefArgoFile);
        linewmo=textscan(ff,'%s %s %s %s');  
        ii = find(strcmp(linewmo{1},wmoIn{1}));
        if isempty(ii)
            warndlg({sprintf('The wmo %d does not exist in the reference list %s.',...
                str2double(wmoIn{1}),CONFIG.RefArgoFile); '';...
                '=> Choose another one or use WOA or IN AIR correction'},...
                sprintf('%s: Uncorrect WMO',CONFIG.history_software'));
            continue
        else
            %Creation de la structure REF_ARGO
            REF_ARGO.wmo   = str2double(cell2mat((linewmo{1}(ii))));
            REF_ARGO.cycle = str2double(cell2mat(linewmo{2}(ii)));
            REF_ARGO.refId = [cell2mat(linewmo{3}(ii)) '_' cell2mat(linewmo{4}(ii))];
        end
    elseif exist(CONFIG.RefArgoFile,'file')==2
        load(CONFIG.bddFile);
        ff=fopen(CONFIG.RefArgoFile);
        linewmo=textscan(ff,'%s %s %s %s');
        ii = find(strcmp(linewmo{1},wmoIn{1}));
        if ~isempty(ii)
            REF_ARGO.wmo   = str2double(cell2mat((linewmo{1}(ii))));
            REF_ARGO.cycle = str2double(cell2mat(linewmo{2}(ii)));
            REF_ARGO.refId = [cell2mat(linewmo{3}(ii)) '_' cell2mat(linewmo{4}(ii))];      
            CONFIG.REF_ARGO=REF_ARGO;
            %--MG ajouter recuperation profil de reference
        end
    end
    
    
    % =================================================================
    %% Open and prepare the argo data
    % =================================================================
    % ******************************************************
    % check if files exist. If not, warn, and reload LOCODOX to the wmo
    % asking for.
    % ******************************************************
    biofiles = dir(fullfile(CONFIG.DataDir,Work.replist,'B*nc'));
    core_files = dir(fullfile(CONFIG.DataDir,Work.replist,'/*nc'));
    core_files = {core_files(:).name}';
    iscore = ~cellfun('isempty',regexp(core_files,'^R[0-9]*')) | ~cellfun('isempty',regexp(core_files,'^D[0-9]*'));
    core_files = core_files(iscore);
    
    if isempty(core_files) && isempty(core_files)
        warndlg(sprintf('No Bio nor Core file in the directory : %s\n \n=> Check CONFIG.DataDir \n',fullfile(CONFIG.DataDir,Work.replist)),'Warning');
        continue
    elseif isempty(biofiles)
        warndlg(sprintf('No Bio file in the directory : %s\n \n=> Check CONFIG.DataDir \n',fullfile(CONFIG.DataDir,Work.replist)),'Warning');
        continue
    elseif isempty(core_files)
        warndlg(sprintf('No Core file in the directory : %s\n \n=> Check CONFIG.DataDir \n',fullfile(CONFIG.DataDir,Work.replist)),'Warning');
        continue
    end
    
    % *****************************************************************
    % Read profile data : if no profil data, warn and stop
    % *****************************************************************
    [DATA.primary.argo, DATA.nearsurf.argo, DATA.secondary.argo,DATA.other.argo, ...
        DATA.primary.Dim, DATA.nearsurf.Dim, DATA.secondary.Dim,DATA.other.Dim,...
        argo, ~] = DOXY_argo_read(CONFIG.DataDir,Work.wmo,Work.replist,Work.whichCorr,Work.dirLog,1);
    if isfield(DATA.primary.argo,'no_inair_data') && strcmp(Work.whichCorr,'INAIR')
        warndlg({sprintf('The wmo %d can''t be corrected with IN AIR method, ',...
            str2double(wmoIn{1})); '';...
            '=> Choose another one or try WOA or REF correction'},...
            'IN AIR correction not possible');
        continue
    end
    
    %Prepare data to compute drift on NCEP %marine 20/06/19
    if CONFIG.ok_inair_drift==1 && ~strcmp(Work.whichCorr,'INAIR')
        Work_NCEP=Work;
        Work_NCEP.whichCorr='INAIR';          
        Work_NCEP.is_not_inair_corr=1;
        [DATA_NCEP.primary.argo, DATA_NCEP.nearsurf.argo, DATA_NCEP.secondary.argo,DATA_NCEP.other.argo, ...
        DATA_NCEP.primary.Dim, DATA_NCEP.nearsurf.Dim, DATA_NCEP.secondary.Dim,DATA_NCEP.other.Dim,...
        argo_NCEP, ~] = DOXY_argo_read(CONFIG.DataDir,Work_NCEP.wmo,Work_NCEP.replist,Work_NCEP.whichCorr,Work_NCEP.dirLog,0);  
 
        if isfield(DATA_NCEP.primary.argo,'no_inair_data')
            no_inair_drift=1;
            h = warndlg({sprintf('No PPOX data or in-air measurement for %d !!!',Work.wmo),...
                ' ',...
                'Impossible to compute drift with NCEP',...
                '=> Drift will be compute on WOA'});
            uiwait(h);        
            clear DATA_NCEP
        end
    end   
    
    
    
    % *****************************************************************
    % Apply a pressure effect correction 
    % *****************************************************************    

    a=questdlg('Do you want to apply a pressure effect correction on this float ?',sprintf('PRESSURE EFFECT CORRECTION - %d',Work.wmo),'Yes','No','No');   
    if strcmp(a,'Yes')
        %figure
        m_outerpos = [0.52 0.15 0.24 0.50];
        fig=figure('unit','normalized','Outerposition',m_outerpos,...
            'Name', sprintf('PRESSURE EFFECT CORRECTION - %d',Work.wmo),'NumberTitle','off','visible','off');
        hold on; grid on; hg1=hggroup;
        if strcmp(argo.VSS,'Near-surface sampling')
            if isempty(DATA.primary.argo.doxy.data)
                plot(DATA.secondary.argo.doxy.data',DATA.secondary.argo.pres.data','b','Parent',hg1)                
            else
                plot(DATA.primary.argo.doxy.data',DATA.primary.argo.pres.data','b','Parent',hg1)
            end
        else
            plot(argo.doxy.data',argo.pres.data','b','Parent',hg1)
        end

        %Apply the pressure effect correction to the data
        Work.anteriorcorr=1;
        Work.coeff_corr = cell2mat(inputdlg(sprintf('Enter coefficient value for pressure effect correction with the formula : \n    "DOXY_corr=DOXY*(1+coeff*PRES/1000)" \n '),sprintf('PRESSURE EFFECT CORRECTION - %d',Work.wmo)));
        if isempty(Work.coeff_corr)
            Work.anteriorcorr=0;
            fprintf('INFO >>>>>   No pressure effect applied \n')
            % ----------    SUMMARY FILE     -------------
            fprintf(log_summ,'No presure effect applied.\n');
            fprintf(log_summ,'\n');
            fprintf(log_summ,'----------------------------------------------------\n');
            fprintf(log_summ,'\n');
        else
            fprintf('INFO >>>>>   A pressure effect is applied \n')
            %apply presure effect
            DATA.primary.argo.doxy.data=DATA.primary.argo.doxy.data.*(1+DATA.primary.argo.pres.data.*str2double(Work.coeff_corr)/1000);
            DATA.secondary.argo.doxy.data=DATA.secondary.argo.doxy.data.*(1+DATA.secondary.argo.pres.data.*str2double(Work.coeff_corr)/1000);
            DATA.nearsurf.argo.doxy.data=DATA.nearsurf.argo.doxy.data.*(1+DATA.nearsurf.argo.pres.data.*str2double(Work.coeff_corr)/1000);
            DATA.other.argo.doxy.data=DATA.other.argo.doxy.data.*(1+DATA.other.argo.pres.data.*str2double(Work.coeff_corr)/1000);
            argo.doxy.data=argo.doxy.data.*(1+argo.pres.data.*str2double(Work.coeff_corr)/1000);
            
            %Add corrected data to the figure
            hg2=hggroup;
            if strcmp(argo.VSS,'Near-surface sampling')
                if isempty(DATA.primary.argo.doxy.data)
                    plot(DATA.secondary.argo.doxy.data',DATA.secondary.argo.pres.data','r','Parent',hg2)
                else
                    plot(DATA.primary.argo.doxy.data',DATA.primary.argo.pres.data','r','Parent',hg2)
                end
            else
                plot(argo.doxy.data',argo.pres.data','r','Parent',hg2)
            end
           
          %  plot(DATA.primary.argo.doxy.data',DATA.primary.argo.pres.data','r','Parent',hg2)
            set(get(get(hg1,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
            set(get(get(hg2,'Annotation'),'LegendInformation'),...hD
                'IconDisplayStyle','on');
            legend('Raw data','Raw data corrected from pressure effect','Location','Southwest')
            title(['Pressure effect correction with coeff=' num2str(Work.coeff_corr)])
            set(gca,'YDir','reverse');
            set(gca,'fontweight','bold');
            xlabel('[O2] (umol/kg)','fontweight','bold','fontsize',Work.fontsize+2);
            ylabel('Pres (dbar)','fontweight','bold','fontsize',Work.fontsize+2);
            set(fig,'visible','on')
            
            warning('From now on, every "raw" data presented will be corrected of the pressure effect');
            
            % ----------    SUMMARY FILE     -------------
            fprintf(log_summ,sprintf('A presure effect was applied : \n  DOXY_corr=DOXY*(1+coeff*PRES/1000)  with coeff=%s\n',num2str(Work.coeff_corr)));
            fprintf(log_summ,'\n');            
            fprintf(log_summ,'----------------------------------------------------\n');
            fprintf(log_summ,'\n');
            
            if Work.savePlot == 1
                hFig=gcf;
                saveFile = fullfile(Work.dirPlot,sprintf('DOXY_pressure_effect_coeff%s_%d',num2str(Work.coeff_corr),Work.wmo));
                DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);       
            end
        end
    else
       Work.anteriorcorr=0;  
       fprintf('INFO >>>>>    No pressure effect applied \n')   
       
       % ----------    SUMMARY FILE     -------------  
       fprintf(log_summ,'No presure effect applied.\n');   
       fprintf(log_summ,'\n');         
       fprintf(log_summ,'----------------------------------------------------\n');     
       fprintf(log_summ,'\n');  
    end
    
      DATA.primary.Work = Work;
      DATA.nearsurf.Work = Work;
      DATA.secondary.Work = Work;
      DATA.other.Work = Work;
    
    if argo.n_prof == 0
        warnstr = {sprintf('%d : the main profile is empty (here : %s). ',REF_ARGO.wmo,argo.VSS);...
            ' ';...
            '=> No correction could be applied.';' '};
        h = warndlg(warnstr,'DOXY_argo_read WARNING');
        uiwait(h);
        continue
    end
    
    
    [argo1Struct, argo2Struct, argo3Struct, argo4Struct, ...
        argoWork, Work] = DOXY_argo_prepare_main(DATA, Work, argo);

    % *****************************************************************
    % Read trajectory data : if no in-air measurement, warn and stop
    % *****************************************************************
    noIAmeas = 0;
    if strcmp(Work.whichCorr,'INAIR')
        [~, ~, argoTrajWork] = DOXY_argoTraj_read(CONFIG,Work.wmo);
        Work.trajVar = argoTrajWork.trajVar;
        argoTrajWork = rmfield(argoTrajWork,'trajVar');
        
        if ~isfield(argoTrajWork,'ppox_doxy_adjusted')
            noIAmeas = 1;        
            h = warndlg({sprintf('No PPOX data for %d !!!',Work.wmo),...
                ' ',...
                'Impossible to use the INAIR correction',...
                '=> Try with WOA or REF correction'},'DOXY_argoTraj_read : no PPOX data');
            uiwait(h);        
            continue
        
        elseif all(cellfun('isempty',argoTrajWork.juld.data))
            noIAmeas = 1;
            h = warndlg({sprintf('No Inair Measurement for %d !!!',Work.wmo),...
                ' ',...
                sprintf('INFO: Measurement code for in-air measurement get in [%s]',num2str(CONFIG.inAirMC)),...
                '',...
                'Impossible to use the INAIR correction',...
                '=> Try with WOA or REF correction'},'DOXY_argoTraj_read : no Inair Measurement');
            uiwait(h);
            continue
        end
    elseif CONFIG.ok_inair_drift==1 && ~exist('no_inair_drift','var') %Prepare data to compute drift on NCEP %marine 20/06/19
        [~, ~, argoTrajWork_NCEP] = DOXY_argoTraj_read(CONFIG,Work_NCEP.wmo);
        Work_NCEP.trajVar = argoTrajWork_NCEP.trajVar;
        argoTrajWork_NCEP = rmfield(argoTrajWork_NCEP,'trajVar');
        
        if ~isfield(argoTrajWork_NCEP,'ppox_doxy_adjusted') || all(cellfun('isempty',argoTrajWork_NCEP.juld.data))
            noIAmeas_NCEP = 1;        
            h = warndlg({sprintf('No PPOX data or in-air measurement for %d !!!',Work.wmo),...
                ' ',...
                'Impossible to compute drift with NCEP',...
                '=> Drift will be compute on WOA'});
            uiwait(h);        
            clear argoTrajWork_NCEP Work_NCEP noIAmeas_NCEP
            no_inair_drift=1;
        else
            Work_NCEP.makePlot=0;  
            Work_NCEP.is_not_inair_corr=1;
            DATA_NCEP.primary.Work = Work_NCEP;
            DATA_NCEP.nearsurf.Work = Work_NCEP;
            DATA_NCEP.secondary.Work = Work_NCEP;
            DATA_NCEP.other.Work = Work_NCEP;
            [argo1Struct_NCEP, argo2Struct_NCEP, argo3Struct_NCEP, argo4Struct_NCEP, ...
        argoWork_NCEP, Work_NCEP] = DOXY_argo_prepare_main(DATA_NCEP, Work_NCEP, argo_NCEP);
        end
    end
    % =====================================================================
    %% Plot argo trajectory
    % =====================================================================
    if CONFIG.M_MAP_ACTIVE==0
        if strcmp(Work.whichCorr,'WOA') || strcmp(Work.whichCorr,'INAIR')
            DOXY_MAP(argo,Work,1)
        else
            iok = strcmp(REF.id,REF_ARGO.refId);
            DOXY_MAP(argo,Work,3,REF,iok)
        end
    end
    
    if CONFIG.M_MAP_ACTIVE>0
        if strcmp(Work.whichCorr,'WOA') || strcmp(Work.whichCorr,'INAIR')
            DOXY_MAP2(CONFIG,argo,Work,1)
        else
            iok = strcmp(REF.id,REF_ARGO.refId);
            DOXY_MAP2(CONFIG,argo,Work,3,REF,iok)
        end
    end
    
    % =====================================================================
    %% Choose the vertical scale : pressure or density ?
    % =====================================================================
    Work.PresOrDens = 'Pressure';

    % =====================================================================
    %% Correction
    % =====================================================================
    % Probably : the DOXY_corr_compute and DOXY_corr_apply_main is the
    % same for woa and ref : keep it outside, and replace DOXY_woa_corr
    % and DOXY_ref_corr by DOXY_woa_corr_prepare and
    % DOXY_ref_corr_prepare
    switch Work.whichCorr
        case 'WOA'
            if CONFIG.ok_inair_drift==1 && ~exist('no_inair_drift','var') 
                NCEP.init = DOXY_NCEP_read(CONFIG,CONFIG.ncepYears);
                [Work_NCEP] = DOXY_inair_prepare_for_drift(CONFIG, WOA,...
    NCEP, argo1Struct_NCEP, argo2Struct_NCEP, argo3Struct_NCEP, argo4Struct_NCEP, argo_NCEP, argoWork_NCEP, Work_NCEP, argoTrajWork_NCEP);
                clear argo1Struct_NCEP argo2Struct_NCEP argo3Struct_NCEP argo4Struct_NCEP 
                Work.NCEP_drift.Work=Work_NCEP;
                Work.NCEP_drift.argo=argo_NCEP;
                Work.NCEP_drift.argoWork=argoWork_NCEP;
                Work.NCEP_drift.argoTrajWork=argoTrajWork_NCEP;
                clear Work_NCEP argo_NCEP argoWork_NCEP argoTrajWork_NCEP
            end
            [argo1Struct,argo2Struct,argo3Struct,argo4Struct,WOA,Work,goProg] = DOXY_woa_corr(CONFIG, ...
                WOA, argo1Struct, argo2Struct, argo3Struct, argo4Struct, argo, argoWork, Work);
        case 'REF'
            if CONFIG.ok_inair_drift==1 && ~exist('no_inair_drift','var') 
                NCEP.init = DOXY_NCEP_read(CONFIG,CONFIG.ncepYears);
                [Work_NCEP] = DOXY_inair_prepare_for_drift(CONFIG, WOA,...
                    NCEP, argo1Struct_NCEP, argo2Struct_NCEP, argo3Struct_NCEP, argo4Struct_NCEP, argo_NCEP, argoWork_NCEP, Work_NCEP, argoTrajWork_NCEP);
                clear argo1Struct_NCEP argo2Struct_NCEP argo3Struct_NCEP argo4Struct_NCEP
                Work.NCEP_drift.Work=Work_NCEP;
                Work.NCEP_drift.argo=argo_NCEP;
                Work.NCEP_drift.argoWork=argoWork_NCEP;
                Work.NCEP_drift.argoTrajWork=argoTrajWork_NCEP;
                clear Work_NCEP argo_NCEP argoWork_NCEP argoTrajWork_NCEP
            end
            
            [argo1Struct,argo2Struct,argo3Struct,argo4Struct,WOA,Work,goProg] = DOXY_ref_corr(CONFIG, ...
                WOA, REF_ARGO, argo1Struct, argo2Struct, argo3Struct, argo4Struct,argo, argoWork, Work);
        case 'INAIR'
            % Read the NCEP data
            NCEP.init = DOXY_NCEP_read(CONFIG,CONFIG.ncepYears);
            [argo1Struct,argo2Struct,argo3Struct,argo4Struct,WOA,NCEP,Work,goProg] = DOXY_inair_corr(CONFIG,...
                WOA, NCEP, argo1Struct, argo2Struct, argo3Struct, argo4Struct, argo, argoWork, Work, argoTrajWork);
    end
    if goProg == 0
        continue;
    end
    
    % =====================================================================
    %% Save structures for wmo
    % =====================================================================
    %marine
    if Work.presEff, presEffStr = 'okpreseff'; else, presEffStr = 'nopreseff'; end
    if Work.DODRIFT
        if strcmp(Work.whichCorr,'INAIR') && strcmp(Work.whichDrift,'NCEP')
            driftStr = 'driftonNCEP';
        elseif strcmp(Work.whichCorr,'REF') && strcmp(Work.whichDrift,'NCEP')
        % Next lines added by Thierry Reynaud 12/03/2020
            driftStr = 'driftonNCEP';
        elseif strcmp(Work.whichCorr,'WOA') && strcmp(Work.whichDrift,'NCEP')
        % Next lines added by Thierry Reynaud 13/03/2020
            driftStr = 'driftonNCEP';
        else
            driftStr = 'driftonWOA';
        end
    else
        driftStr = 'nodrift';
    end
    if Work.(['FIT_intercept_' Work.whichO2quantity]) == 0, offsetStr = 'nooffset';
    else, offsetStr = 'okoffset';
    end
    if isfield(argo1Struct.Work,['OFFSET_' argo1Struct.Work.whichO2quantity])
        corr_type=[Work.whichCorr 'constant'];
    else
        corr_type=Work.whichCorr;
    end
    Work.dirout=sprintf('%s_%s_%s_%s_%s',corr_type,presEffStr,driftStr,offsetStr,Work.whichO2quantity);
    
    if ~exist(fullfile(CONFIG.saveDataDir,'MAT',Work.dirout),'dir')
        mkdir(fullfile(CONFIG.saveDataDir,'MAT',Work.dirout));
    end 
    save(fullfile(CONFIG.saveDataDir,'MAT',Work.dirout,...
        sprintf('%d_corr',Work.wmo)),'argo','argoWork','Work',...
        'argo1Struct','argo2Struct','argo3Struct','argo4Struct');
    
    % =====================================================================
    %% Write NetCDF corrected
    % =====================================================================
    % Prepare the data to be writen
    [argo1Struct, argo2Struct, argo3Struct, argo4Struct, CONFIG] = DOXY_write_prepare(argo1Struct, ...
        argo2Struct, argo3Struct, argo4Struct, CONFIG, Work);
      
    % Write the new monoprofiles
    if strcmp(Work.whichCorr,'REF')
        DOXY_argo_write(CONFIG.DataDir,Work.replist,CONFIG.saveDataDir,...
        CONFIG.prefix,argo1Struct, argo2Struct, argo3Struct, argo4Struct, Work,REF_ARGO);  
    else
        DOXY_argo_write(CONFIG.DataDir,Work.replist,CONFIG.saveDataDir,...
        CONFIG.prefix,argo1Struct, argo2Struct, argo3Struct, argo4Struct, Work);  
    end
    
    % --------   SUMMARY FILE   --------------
    fprintf(log_summ,'\n');
    fprintf(log_summ,'             CORRECTION COMPLETED\n');
    
    
    % close the summary file
    fclose('all');

    newWMO = questdlg('Do you want to manage a new wmo ?',...
        sprintf('%s: Continue the correction', CONFIG.history_software),'YES',...
        sprintf('NO, STOP %s',CONFIG.history_software),sprintf('NO, STOP %s',CONFIG.history_software));
    if goProg == 0, 
        goProg =1; 
    end
    if strcmp(newWMO,sprintf('NO, STOP %s',CONFIG.history_software))
        goProg = 0;
        %msgbox({'  Good Bye !  ';''},sprintf('%s',CONFIG.history_software),'custom',imread(CONFIG.logo));
    else
        if strcmp(Work.whichCorr,'INAIR') && ~noIAmeas
            NCEP = rmfield(NCEP,'coloc');
        end
        Work = initialWork;
        if length(fieldnames(WOA)) == 2
            WOA = rmfield(WOA,'interp');
            clear DATA argo1Struct argo2Struct argo3Struct argo argoWork;
        end
    end
end


