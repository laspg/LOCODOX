% DOXY_corr_apply_main apply correction on data.
%
% SYNTAX
% [argo1Struct, argo2Struct, argo3Struct] = DOXY_corr_apply_main(hFig,
% CORR,argo1Struct, argo2Struct, argo3Struct, argo4Struct, Work, argoTrajWork)
%
% DESCRIPTION
% DOXY_corr_apply_main prepares data for correction and apply correction.
% The correction applied is linear correction, and Constant Correction if
% needed.
%
% INPUT
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
%   argo1Struct (struct)    Three data structures organised with the main profile
%   argo2Struct (struct)    for working is in argo1Struct (vertical
%   argo3Struct (struct)    sampling shceme of interest for the kind of
%   argo4Struct (struct)    correction), and the others are in the second ones.
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
% 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%                                       sat: [1x1 struct]
%                                      psat: [1x1 struct]
%
%                          argo1Struct.Work =
% 
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
%                                   DOXY_QC: [165x117 char]
%                                      DENS: [165x117 single]
%                                           ...
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
%   argoTrajWork (structure)     float data structure (argoTrajWork) choose for the Doxy
%                               correction.
%                       Example: 
%                           argoTrajWork =
%                          float_serial_no: [1x1 struct]
%                            wmo_inst_type: [1x1 struct]
%                                     juld: [1x1 struct]
%                                 latitude: [1x1 struct]
%                                longitude: [1x1 struct]
%
%
% OUTPUT
%   argo1Struct (struct)    Three data structures organised with the main profile
%   argo2Struct (struct)    for working is in argo1Struct (vertical
%   argo3Struct (struct)    sampling shceme of interest for the kind of
%   argo4Struct (struct)    correction), and the others are in the second ones.
%                           The DOXY data are corrected with the correction
%                           calculated in Work. The correction has been
%                           computed using main profile, ascending
%                           profiles. The correction is applied to both
%                           ascending and descending profiles, and is
%                           propagated to every vertical sampling scheme.
%
% CALL : 
%   DOXY_corr_prepare, DOXY_corr_apply, DOXY_corr_cst, DOXY_PLOT_corr
%
% SEE ALSO
%   DOXY_corr_compute, DOXY_corr_main
%
% HISTORY
%   $created: 14/01/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v2 25/05/2016   Emilie Brion, ALTRAN OUEST
%                       takes into account three Vertical Sampling Scheme
%                       (primary, secondary and nearsurface).
%          25/03/2019   Marine GALLIAN, ALTRAN OUEST 
%                       Add another figure representing corrected data and
%                       reference profile
%                       Fill the summary file

function [argo1Struct, argo2Struct,argo3Struct,argo4Struct] = ...
    DOXY_corr_apply_main(hFig, CORR, argo1Struct, argo2Struct,argo3Struct,argo4Struct, Work, argoTrajWork)


% =========================================================================
%% Initialisation
% =========================================================================
% arguments
minArgs = 7;
maxArgs = 8;
narginchk(minArgs,maxArgs);

if Work.savePlot == 1
    Work.savePlot = 0;
    mem = 1;
else
    mem = 0;
end

% =========================================================================
%% Apply the Linear regression : correction of DOXY data
% =========================================================================
Work.isargo1Struct=1;
% Get Ascendant and descendant profiles, main vertical scheme
if ismember(Work.whichCorr,{'WOA','REF'})
    [CORR1] = DOXY_corr_prepare(CORR, Work, argo1Struct);
else
    [CORR1] = DOXY_corr_prepare(CORR, Work, argo1Struct, argoTrajWork);
end
Work=rmfield(Work,'isargo1Struct');

if ismember(Work.whichCorr,{'WOA','REF'})
    CORR1.LRcoef = CORR.LRcoef;
elseif ismember(Work.whichCorr,{'INAIR'})
    CORR1.NLRcoef = CORR.NLRcoef;
    CORR1.NLRinfo = CORR.NLRinfo;
end

% Apply correction on the ascendant and descendant profiles, main vertical
% scheme
[CORR1] = DOXY_corr_apply(CORR1,argo1Struct.argoWork,Work.(['FIT_slope_' CORR.whichO2quantity]),...
    Work.(['FIT_intercept_' CORR.whichO2quantity]));
argo1Struct.Work.(['PSAT_LINCORR_' CORR.whichO2quantity]) = CORR1.psat_corr;
argo1Struct.Work.(['DOXY_LINCORR_' CORR.whichO2quantity]) = CORR1.doxy_corr;
% carry the correction on the doxy_adjusted field
argo1Struct.argoWork.doxy_adjusted.data = argo1Struct.Work.(['DOXY_LINCORR_' CORR1.whichO2quantity]);

% Get ascendant and descendant profiles and apply correction, secondary
% vertical schemes
for n=2:4
    eval(sprintf('argoStruct = argo%dStruct;',n));
    if argoStruct.argo.n_prof ~= 0
        % propagate the correction to the secondary profile
        if ismember(Work.whichCorr,{'WOA','REF'})
            [CORRn] = DOXY_corr_prepare(CORR, Work,argoStruct);
        else
            [CORRn] = DOXY_corr_prepare(CORR, Work,argoStruct, argoTrajWork);
        end
        if ismember(Work.whichCorr,{'WOA','REF'})
            CORRn.LRcoef = CORR.LRcoef;
        elseif ismember(Work.whichCorr,{'INAIR'})
            CORRn.NLRcoef = CORR.NLRcoef;
            CORRn.NLRinfo = CORR.NLRinfo;
            CORRn.NLRstreq = CORR.NLRstreq;
        end
        [CORRn] = DOXY_corr_apply(CORRn,argoStruct.argoWork,Work.(['FIT_slope_' CORR.whichO2quantity]),...
                Work.(['FIT_intercept_' CORR.whichO2quantity]));
        argoStruct.Work.(['PSAT_LINCORR_' CORRn.whichO2quantity]) = CORRn.psat_corr;
        argoStruct.Work.(['DOXY_LINCORR_' CORRn.whichO2quantity]) = CORRn.doxy_corr;
        argoStruct.argoWork.doxy_adjusted.data = argoStruct.Work.(['DOXY_LINCORR_' CORRn.whichO2quantity]);
        eval(sprintf('CORR%d = CORRn;',n));    
        eval(sprintf('argo%dStruct = argoStruct;',n));    
    end
end

% Get the correct levels (Pres or Dens)
if strcmpi(CORR.presOrDens, 'pressure')
    CORR1.level = argo1Struct.argoWork.pres_adjusted.data;
    if exist('CORR2','var')
        CORR2.level = argo2Struct.argoWork.pres_adjusted.data;
    end
    if exist('CORR3','var')
        CORR3.level = argo3Struct.argoWork.pres_adjusted.data;
    end
    if exist('CORR4','var')
        CORR4.level = argo4Struct.argoWork.pres_adjusted.data;
    end
elseif strcmpi(CORR.presOrDens, 'density')
    CORR1.level = argo1Struct.Work.DENS;
    if isfield(CORR2,'argoWork')
        CORR2.level = argo2Struct.Work.DENS;
    end
    if isfield(CORR3, 'argoWork')
        CORR3.level = argo3Struct.Work.DENS;
    end
    if isfield(CORR4, 'argoWork')
        CORR4.level = argo4Struct.Work.DENS;
    end
end

% Control Plots
if Work.makePlot
    if ismember(Work.whichCorr,{'WOA','REF'})
        % Nb of day after deployment
        dayjul = argo1Struct.argo.juld.data - argo1Struct.argo.juld.data(1);
        % Plots
        argo1Struct.Work.savePlot=0;
        DOXY_PLOT_corr(hFig, 4, CORR1, argo1Struct.Work, dayjul);
        DOXY_PLOT_corr(hFig, 5, CORR1, argo1Struct.Work);
        DOXY_PLOT_corr(hFig, 6, CORR1, argo1Struct.Work);
        
    elseif ismember(Work.whichCorr,{'INAIR'})
        % keep primary sampling for plots
        if argo2Struct.argo.n_prof ~= 0
            argoStruct = argo2Struct;
            CORRtmp = CORR2;
        elseif argo2Struct.argo.n_prof == 0 && argo3Struct.argo.n_prof ~= 0
            argoStruct = argo3Struct;
            CORRtmp = CORR3;
        end            
   
        % Nb of day after deployment
        dayjul = argoStruct.argo.juld.data - argoStruct.argo.juld.data(1);
                
        % Plots
        % plot primary argo raw profile, primary argo corrected profile,
        % and WOA profile (PSAT and DOXY)        
        % keep primary sampling for plots
        if mem == 1
            argoStruct.Work.savePlot = 0;
        end
        DOXY_PLOT_corr(hFig, 16, CORRtmp, argoStruct.Work);
        DOXY_PLOT_corr(hFig, 4, CORRtmp, argoStruct.Work, dayjul,7);
        DOXY_PLOT_corr(hFig, 5, CORRtmp, argoStruct.Work,8);
        if mem == 1
            Work.savePlot = 1;
            argoStruct.Work.savePlot = 1;
        end
        DOXY_PLOT_corr(hFig, 6, CORRtmp, argoStruct.Work,9);        
    end
end

% 
% =========================================================================
%% Constant Correction for R² to small : offset with WOA
% =========================================================================
if ismember(Work.whichCorr,{'WOA','REF'})
    if Work.(['FIT_rsquare_' CORR.whichO2quantity]) < Work.R2threshold
        input = questdlg({sprintf('The R2 of linear correction is under the threshold defined in the configuration file : R2 = %2.2f  <  %2.2f',Work.(['FIT_rsquare_' CORR.whichO2quantity]), Work.R2threshold);...
            '  ';'=> LOCODOX suggests to apply the constant correction instead of the correcton based on the linear correction.';'   '},...
            'APPLY CONSTANT CORRECTION ?','YES','NO','NO');
        
        % -----------------      SUMMARY FILE    ---------------------------
        fprintf(Work.log_summ,'\n');
        if strcmp(input,'NO')
            fprintf(Work.log_summ,...
                sprintf('The R-squared indicator is equal to %s and is lower than the threshold define by the user in the configuration file (=%s) \n --> THE USER CHOSE TO STILL APPLY THE LINEAR REGRESSION \n',...
                num2str(Work.(['FIT_rsquare_' CORR.whichO2quantity])),num2str(Work.R2threshold)));
        else
            fprintf(Work.log_summ,...
                sprintf('The R-squared indicator is equal to %s and is lower than the threshold define by the user in the configuration file (=%s) \n --> THE USER CHOSE TO APPLY A CONSTANT CORRECTION INSTEAD\n',...
                num2str(Work.(['FIT_rsquare_' CORR.whichO2quantity])),num2str(Work.R2threshold)));
        end
    else
        % -----------------      SUMMARY FILE    ---------------------------
        fprintf(Work.log_summ,...
            sprintf('The R-squared indicator is equal to %s and is greater than the threshold define by the user in the configuration file (=%s) \n --> THE LINEAR REGRESSION IS APPLIED\n',...
            num2str(Work.(['FIT_rsquare_' CORR.whichO2quantity])),num2str(Work.R2threshold)));
        input = 'NO';
    end
    
    if strcmp(input,'YES')
        fprintf('INFO >>>>>    Constant correction is applied : R² = %2.2f, < %2.2f\n',Work.(['FIT_rsquare_' CORR.whichO2quantity]), Work.R2threshold)
        % -----------------------------------------------------------------
        % Offset computation between argo and WOA
        % -----------------------------------------------------------------
        diff = CORR.([CORR.strField 'ref']) - CORR.(CORR.strField);
        dec = nanmean(diff(:));
        
        % -----------------------------------------------------------------
        % Offset correction
        % -----------------------------------------------------------------
        % Apply offset correction on the ascendant and descendant profiles,
        % main vertical scheme
        [argo1Struct.Work] = DOXY_corr_cst(CORR1,argo1Struct.Work,argo1Struct.argoWork,dec);
        if argo1Struct.Work.savePlot==1
            argo1Struct.Work.savePlot=0;
            mem=1;
        else
            mem=0;
        end
        
        if Work.makePlot
            DOXY_PLOT_corr(hFig, 7,CORR, argo1Struct.Work, argo1Struct.argoWork);
        end
        
        % carry the correction on the doxy_adjusted field
        argo1Struct.argoWork.doxy_adjusted.data = argo1Struct.Work.(['DOXY_OFFSETCORR_' CORR.whichO2quantity]);
        argo1Struct.Work.drift=CORR1.drift;
        
        % Apply offset correction on the ascendant and descendant profiles,
        % other vertical shemes
        if argo2Struct.argo.n_prof ~= 0
            CORR2.no_log=1;
            [argo2Struct.Work] = DOXY_corr_cst(CORR2,argo2Struct.Work,argo2Struct.argoWork,dec);
            % carry the correction on the doxy_adjusted field
            argo2Struct.argoWork.doxy_adjusted.data = argo2Struct.Work.(['DOXY_OFFSETCORR_' CORR.whichO2quantity]);
        end
        if argo3Struct.argo.n_prof ~= 0
            CORR3.no_log=1;
            [argo3Struct.Work] = DOXY_corr_cst(CORR3,argo3Struct.Work,argo3Struct.argoWork,dec);
            % carry the correction on the doxy_adjusted field
            argo3Struct.argoWork.doxy_adjusted.data = argo3Struct.Work.(['DOXY_OFFSETCORR_' CORR.whichO2quantity]);
        end
        if argo4Struct.argo.n_prof ~= 0
            CORR4.no_log=1;
            [argo4Struct.Work] = DOXY_corr_cst(CORR4,argo4Struct.Work,argo4Struct.argoWork,dec);
            % carry the correction on the doxy_adjusted field
            argo4Struct.argoWork.doxy_adjusted.data = argo4Struct.Work.(['DOXY_OFFSETCORR_' CORR.whichO2quantity]);
        end
        
        % -----------------------------------------------------------------
        % Control plot about Constante correction
        % -----------------------------------------------------------------
        if Work.makePlot
            DOXY_PLOT_corr(hFig, 8,CORR1, argo1Struct.Work, argo1Struct.argoWork );
            if mem == 1
                Work.savePlot = 1;
                argo1Struct.Work.savePlot = 1;
            end
            DOXY_PLOT_corr(hFig, 9,CORR1, argo1Struct.Work, argo1Struct.argoWork );
        end
    % if no constant correction 
    else
        argo1Struct.Work.drift=CORR1.drift;
        fprintf('\t No need to apply constant correction : R² = %2.3f, >= %2.2f\n',...
            Work.(['FIT_rsquare_' CORR.whichO2quantity]),Work.R2threshold);
        if Work.makePlot
            DOXY_PLOT_corr(hFig, 7,CORR1, Work, 0);
            DOXY_PLOT_corr(hFig, 8,CORR1, Work, 0);
            if mem == 1
                Work.savePlot = 1;
            end
            DOXY_PLOT_corr(hFig, 9,CORR1, Work, 0);
        end
    end
    
    % -----------------      SUMMARY FILE    ---------------------------
    fprintf(Work.log_summ,'\n');
    fprintf(Work.log_summ,'-----------------------------------------------------------\n');
    
    
    % -----------------------------------------------------------------
    % Plot corrected data with the constant correction and reference profile
    % -----------------------------------------------------------------
    %Create another fig for plotting corrected data --MG
    if strcmp(Work.whichCorr,'REF') ||  isfield(argo1Struct.Work,'plotREF')
        m_outerpos = [0.52 0.15 0.24 0.50];
        m_Fig = figure('unit','normalized','Outerposition',m_outerpos,...
            'Name', sprintf('DATA CORRECTED AND REFERENCE PROFILE - %d',Work.wmo),'NumberTitle','off');
        if mem==1
            argo1Struct.Work.savePlot=1;
        else
            argo1Struct.Work.savePlot=0;
        end
        DOXY_PLOT_corr(m_Fig, 10, CORR1, argo1Struct.Work);
    else
        % Added by T Reynaud 06.04.2020
         m_outerpos = [0.52 0.15 0.24 0.50];
        m_Fig = figure('unit','normalized','Outerposition',m_outerpos,...
            'Name', sprintf('DATA CORRECTED without REFERENCE PROFILE - %d',Work.wmo),'NumberTitle','off');
        if mem==1
            argo1Struct.Work.savePlot=1;
        else
            argo1Struct.Work.savePlot=0;
        end
        DOXY_PLOT_corr(m_Fig, 10, CORR1, argo1Struct.Work);
    end
    
else 
    argo2Struct.Work.drift=CORR2.drift;       
    if isfield(argoStruct.Work,'plotREF')
        m_outerpos = [0.52 0.15 0.24 0.50];
        m_Fig = figure('unit','normalized','Outerposition',m_outerpos,...
            'Name', sprintf('DATA CORRECTED AND REFERENCE PROFILE - %d',Work.wmo),'NumberTitle','off');
        if mem==1
            argoStruct.Work.savePlot=1;
        else
            argoStruct.Work.savePlot=0;
        end
        DOXY_PLOT_corr(m_Fig, 10, CORRtmp, argoStruct.Work,9);
    else
        % added by Thierry Reynaud 10.04.2020
        m_outerpos = [0.52 0.15 0.24 0.50];
        m_Fig = figure('unit','normalized','Outerposition',m_outerpos,...
            'Name', sprintf('DATA CORRECTED WITHOUT REFERENCE PROFILE - %d',Work.wmo),'NumberTitle','off');
        if mem==1
            argoStruct.Work.savePlot=1;
        else
            argoStruct.Work.savePlot=0;
        end
        DOXY_PLOT_corr(m_Fig, 10, CORRtmp, argoStruct.Work,9);        
    end
end
