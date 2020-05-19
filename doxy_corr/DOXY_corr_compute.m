% DOXY_corr_compute compute the correction over DOXY and PSAT.
%
% SYNTAX
% [Work, CORR, hFig] = DOXY_corr_compute(Work, CORR, hFig)
%
% DESCRIPTION
% DOXY_corr_compute cover the correction on DOXY and PSAT, following the
% method ATLN.
%
% INPUT
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
%                           Work.DOXY_WOA =
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
%  hFig (handle)        handle of the figure with correction plots
%                       (DOXY_PLOT_corr).
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
%                            whichO2quantity: 'PSAT'
%                            strField: 'psat'
%                                psat: [164x117 single]
%                             psatref: [164x117 double]
%                             psat_ok: [164x117 double]
%                                   ...
%
%   hFig (handle)        handle of the figure with correction plots
%                        (DOXY_PLOT_corr).
%
% CALL : 
%   DOXY_PLOT_corr
%
% SEE ALSO
%   DOXY_corr_main, DOXY_corr_linreg

% HISTORY
%   $created: 11/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v2,   30/05/2016, Emilie Brion, Altran Ouest
%              break in DOXY_corr_compute_main and DOXY_corr_compute, with
%              optional choice
%       v3    26/01/2017  Emilie Brion, Anne Piron, Altran Ouest
%                         argo - ameliorations 2017 phase 1
%       v3.2  12/09/2017  Emilie Brion, Altran Ouest
%                         If the doxy (or psat) carries too much NaN (more
%                         than 20%), the level scale is turned to PRES
%                         instead of DEPTH to avoid gradO2 filled with too
%                         much NaN.
%             25/03/2019  Marine GALLIAN, Altran ouest 
%                         Take in acount when the Carry over parameter is
%                         set to 0 in the configuration file (isokC=0)
%                         Write information in the log file
%                         Now we use only the PRES variable to compute the
%                         grad02

function [Work, CORR, hFig] = DOXY_corr_compute(Work, CORR, hFig)
%Messages for log file
msg=strings;
cpt_msg=0;

if ismember(Work.whichCorr,{'WOA','REF'})
    % =====================================================================
    %% WOA or INSITU REF correction
    % Find the linear regression : doxy VS doxy WOA
    % =====================================================================
    field = CORR.strField;
    
    % ---------------------------------------------------------------------
    % Keep data with slope dO2/dz < 0.2
    % ---------------------------------------------------------------------
    % Check the gradient : if dVAR/dZ > 0.2, bad data are put to NaN. Bad
    % data values are kept in the table :
    % ([field '_ko']) / ([field 'ref_ko']).
    % ---------------------------------------------------------------------
    for z = 2:size(CORR.level,2)
        for i = 1:length(CORR.refCyc)
            % Test the gradient of O2 or PSAT over depth : dVAR/dZ >0.2 (Z
            % is pres) --MG
            gradO2(i,z) = ((CORR.(field)(i,z) - CORR.(field)(i,z-1)))./((Work.PRES(CORR.refCyc(i),z)-Work.PRES(CORR.refCyc(i),z-1)));
            if abs(gradO2(i,z)) > 0.2
                CORR.([field '_ok'])(i,z) = NaN;
                CORR.([field 'ref_ok'])(i,z) = NaN;
                CORR.([field '_ko'])(i,z) = CORR.(field)(i,z);
                CORR.([field 'ref_ko'])(i,z) = CORR.([field 'ref'])(i,z);
            elseif abs(gradO2(i,z)) <= 0.2 || isnan(gradO2(i,z)) %marine 19/06/19
                CORR.([field '_ok'])(i,z) = CORR.(field)(i,z);
                CORR.([field 'ref_ok'])(i,z)= CORR.([field 'ref'])(i,z);
            end
        end
    end
    
    % Choose the kind of linear regression : y ~ a*x +b or y ~ a*x
    % The second case assumes that the optode is calibrated so that it
    % mesures practically 0 in a no-DOXY solution.
    %   - Make a plot
    %   - the user choose the regression
    %   - make the correction
    %
    % Linear regression : the linear regression formula is  y ~ a*x + b
    % However, against the kind of optode calibration (multipoint, ...),
    % the operator should choose to force b to 0 (linear regression with
    % zero intercept) or to keep it free. LOCODOX ask for it.    
    myVar = CORR.(field);
    myVarRef = CORR.([field 'ref']);
    isok = ~isnan(myVar)  &  ~isnan(myVarRef);
    
    linReg.equation{1} = 'y ~ a*x';
    linReg.equation{2} = 'y ~ a*x + b';
    linReg.equation_for_plot{1} = sprintf('%s_adjusted ~ slope*%s',field,field); %--MG
    linReg.equation_for_plot{2} = sprintf('%s_adjusted ~ slope*%s + offset',field,field);    
    linReg.intercept(1) = false;
    linReg.intercept(2) = true;     
    [CORR.LRfit, CORR.LRFullFit, CORR.LRcoef] = DOXY_corr_linreg(myVar,myVarRef,isok,linReg);
    
    % Control plot
    if Work.makePlot
        fprintf('INFO >>>>>    Correction :\n') %2/07/19
        DOXY_PLOT_corr(hFig, 0, CORR, Work,linReg,1);
    end
    
    
    [whichRegCoef,answer] = listdlg('PromptString',...
        {sprintf('A linear regression between the vertical profile of %s measured by the float and the reference was estimated.',field);...
        'Now, you can choose between:';...
        sprintf('   # %s_adjusted ~ slope*%s : regression forced to the origin (0,0)',field,field);...
        sprintf('   # %s_adjusted ~ slope*%s +offset : keep the y intercept estimated by the linear regression',field,field);...
        ''},...
        'SelectionMode','single',...
        'ListString',linReg.equation_for_plot,'ListSize',[500 120],'Name','DOXY CORRECTION');
    
    %prepare msg for log file --MG
    cpt_msg=cpt_msg+1;
    msg(cpt_msg)=sprintf('A linear regression between the vertical profile of %s measured by the float and the reference was estimated.  ',field);
    if whichRegCoef==1
        cpt_msg=cpt_msg+1;
        msg(cpt_msg)= sprintf('The user choose to compute a regression forced to the origin (0,0) :  %s_adjusted ~ slope*%s \n',field,field);
    else
        cpt_msg=cpt_msg+1;
        msg(cpt_msg)= sprintf('The user choose to keep the y intercept estimated by the linear regression : %s_adjusted ~ slope*%s +offset \n',field,field);
    end
    
    
    if answer == 0
        return
    end
    
    clf;
    
    [CORR.LRfit, CORR.LRFullFit,CORR.LRcoef] = DOXY_corr_linreg(myVar,myVarRef,isok,...
                                      linReg,0,1:length(linReg.intercept));
    CORR.LRfit = CORR.LRfit(whichRegCoef,:);
    CORR.LRFullFit = CORR.LRFullFit(whichRegCoef,:,:);
    CORR.LRcoef = CORR.LRcoef(whichRegCoef,:);
    
    % model II linear regression for data whose gradient < 0.2
    myVar = CORR.([field '_ok']);
    myVarRef = CORR.([field 'ref_ok']);
    isok = ~isnan(myVar)  &  ~isnan(myVarRef);
    [LRfit, LRFullFit, LRcoef] = DOXY_corr_linreg(myVar,myVarRef,isok,linReg,0,whichRegCoef);
    CORR.LRfit(2,:) = LRfit;
    CORR.LRFullFit(2,:,:) = LRFullFit;
    CORR.LRcoef(2,:) = LRcoef;    
    CORR.linReg = linReg.equation{whichRegCoef};
    % Control plot
    if Work.makePlot
        DOXY_PLOT_corr(hFig, 1, CORR, Work);
    end

    % ---------------------------------------------------------------------
    % Manage the data too much different
    % ---------------------------------------------------------------------
    % Difference between WOA and the corrected float
    tmp = squeeze(CORR.LRFullFit(2,:,:));
    if size(tmp,1) == size(CORR.([field 'ref_ok']),2)
        tmp = tmp'; 
    end
    Diff = CORR.([field 'ref_ok']) - tmp;
%     Diff = CORR.([field 'ref_ok']) - squeeze(CORR.LRFullFit(2,:,:));
    
    isnok = isnan(Diff);
    tmpStd = nanmean(nanstd(Diff(~isnok)));
    avg = nanmean(nanmean(Diff));
    CORR.Diff = Diff;
    CORR.([field '_std']) = tmpStd;
    CORR.([field '_avg']) = avg;
    
    if Work.makePlot
        DOXY_PLOT_corr(hFig, 2, CORR, Work);
    end
    
    % Eliminate data too much different
    tt = 0;
    while isempty(find(abs(Diff) > avg + 2.8*tmpStd, 1)) == 0
        tt = tt + 1;
        tmpStd = nanmean(nanstd(Diff)); avg = nanmean(nanmean(Diff)); Diff(abs(Diff) > avg + 2.8*tmpStd) = NaN;
        if tt > length(Diff)
            cprintf('err', 'WARNING : the loop is too long!\n');
            break
        end
    end
    
    CORR = rmfield(CORR,([field '_std']));
    CORR = rmfield(CORR,([field '_avg']));
    
    % ---------------------------------------------------------------------
    % Linear Regression over good data
    % ---------------------------------------------------------------------
    isnok = logical(isnan(Diff) - isnok);
    myVar = CORR.([field '_ok'])(~isnok);
    myVarRef = CORR.([field 'ref_ok'])(~isnok);
    isok = ~isnan(myVar)  &  ~isnan(myVarRef);
    
    [LRfit, LRFullFit, LRcoef] = DOXY_corr_linreg(myVar,myVarRef,isok,linReg,0,whichRegCoef);
    CORR = rmfield(CORR,'LRFullFit');
    CORR.LRfit(3,:) = LRfit;
    CORR.LRcoef(3,:) = LRcoef;
        
    % Control plot
    if Work.makePlot
        DOXY_PLOT_corr(hFig, 3,CORR,Work,isnok);
    end
    
    % ---------------------------------------------------------------------
    % Keep information on Linear Regression
    % ---------------------------------------------------------------------
    Work.(['FIT_slope_' CORR.whichO2quantity]) = CORR.LRcoef(3,1);% m;
    Work.(['FIT_intercept_' CORR.whichO2quantity]) = CORR.LRcoef(3,2); %b;
    Work.(['FIT_rsquare_' CORR.whichO2quantity]) = CORR.LRcoef(3,3); %r*r;
    
    % -------------- SUMMARY FILE -----------------%--MG
    fprintf(Work.log_summ,'\n');
    if strcmp(Work.whichO2quantity,'DOXY')
        fprintf(Work.log_summ,'Correction computed on DOXY\n');
    else 
        fprintf(Work.log_summ,'Correction computed on PSAT\n');    
    end
    for i_msg=1:cpt_msg
        fprintf(Work.log_summ,msg(i_msg));
    end
elseif ismember(Work.whichCorr,{'INAIR'})
    % =====================================================================
    %% INAIR correction : Find the non-linear regression to correct PPOX
    % ppox inflated as function of ppox deflated and ppox NCEP :
    % PO2inflated = c*PO2deflated + [(1-c)/m]PO2air
    % Bittig, 2015
    % =====================================================================
    % Manage the license
    [status,errmsg] = license('checkout','statistics_toolbox');
    if status == 0
        warndlg({'Statistics Toolbox unavailable !!','=> INAIR correction impossible',...
                 'Try the WOA or the REF correction'},'DOXY_corr_compute')
    end
    
    isok = ~isnan(CORR.ncep_pO2_scaled);
    if Work.isokC==1
        slopefun=@(beta,X)beta(2).*X(:,1)+(1-beta(2))./(1+beta(1)/100).*X(:,2);
        strfun = 'PO2_inair ~ c*PO2_inwater + ((1-c)/m)PO2_ncep';
        % compute the coefficient c and m
        [beta,mdl] = nlinfitci([CORR.po2_inwater_sizeInair(isok)' ...
        CORR.ncep_pO2_scaled(isok)'],CORR.po2_inair(isok)',slopefun,[0;0.23]);
    
    elseif Work.isokC==0
        slopefun=@(beta,X)(1)./(1+beta/100).*X; %--MG
        strfun = 'PO2_inair ~ (1/m)PO2_ncep'
        % compute the coefficient c and m
        [beta,mdl] = nlinfitci(CORR.ncep_pO2_scaled(isok)',CORR.po2_inair(isok)',slopefun,[0]); %--MG    
    end
    
    
    CORR.NLRcoef = beta;
    CORR.NLRcoef(:,3) = mdl.RMSE;
    CORR.NLRinfo = mdl;
    CORR.NLRstreq = strfun;
    
    if Work.savePlot == 1
        Work.savePlot = 0;
        mem=1;
    else
        mem=0;
    end
    
    if Work.makePlot
        DOXY_PLOT_corr(hFig,11,CORR,Work);
        DOXY_PLOT_corr(hFig,13,CORR,Work);
        DOXY_PLOT_corr(hFig,14,CORR,Work);
    end
    
    % ---------------------------------------------------------------------
    % Keep information on Linear Regression
    % ---------------------------------------------------------------------
    Work.(['FIT_slope_' CORR.whichO2quantity]) = CORR.NLRcoef(1,1)/100 + 1;
    Work.(['FIT_intercept_' CORR.whichO2quantity]) = 0;
    Work.(['FIT_rsquare_' CORR.whichO2quantity]) = CORR.NLRcoef(1,3);
    Work.NLRcoef = CORR.NLRcoef;
    Work.NLRinfo = CORR.NLRinfo;
    Work.NLRstreq = CORR.NLRstreq;
    
    %--MG
    if mem==1
        Work.savePlot=1;
    else
        Work.savePlot=0;
    end
    
end