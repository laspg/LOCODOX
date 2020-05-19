% DOXY_PLOT_corr do the control plots for the correction step of LOCODOX
%
% SYNTAX
% [hFig] = DOXY_PLOT_corr(varargin)
%
% DESCRIPTION
% DOXY_PLOT_corr do the control plots for the correction step of LOCODOX.
% The control plot slighlty differ from one correction type to another one.
%
% INPUT
% the four first argument are the following:
%         hFig = varargin{1};
%         plotTyp = varargin{2};
%         CORR = varargin{3};
%         Work = varargin{4};
%
%   hFig (handle)       figure handles of the plot.
%
%   plotTyp (double)    make the plot following the case
%                       "plotTyp".
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
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA)
%                        Example:
%                        Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                                       wmo: 1901205
%
% Optionnal argin should be call for different cases:
%      linReg           regression equation description, and Z intercept
%                        information. linReg.intercept = false means that
%                        the regression equation is forced to origin.
%                        Example: linReg = {'y ~ a*x + b','y ~ a*x'};
%                         linReg =
%                             equation: {'y ~ a*x + b'  'y ~ a*x'}
%                             intercept: [1 0]
%     fitted0 (logical)  Use the minimum in data (0) or the
%                        value O (1) to build the fitted data. Default
%                        = 0 (use the minimum in data).
%
% OUTPUT
%   hFig (handle)       figure handles of the plot.
%
% CALL :
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 14/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v2 25/05/2016   Emilie Brion, ALTRAN OUEST
%                       adapted to INAIR correction, debugging
%       v3 26/01/2017   Emilie Brion, Altran Ouest
%                       argo - ameliorations 2017 phase 1
%       v3.1 09/11/2017 Emilie Brion, Altran Ouest
%                       take into account po2_air has becomed an array (no
%                       more a vector). Use Work.DODRIFT to now make the
%                       plot name when saving.
%           22/02/2019  Marine GALLIAN, Altran Ouest
%                       Take in account QC of data : Now for data corrected
%                       only points with QC<4 are plot
%                       Improve figure for in air correction : plot mean of
%                       data and std, display stats ...
%          13/03/2020   Thierry Reynaud 13/03/2020
%                       Work.whichDrift introduced for filename
%                       construction line 581
%          05.04.2020   Modifications by T. Reynaud for plotting and
%                       generating corrected oxy profiles figures without   
%                       a REF profile

%
function [hFig] = DOXY_PLOT_corr(varargin)

% =========================================================================
%% Initialisation
% =========================================================================
hFig = varargin{1};
plotTyp = varargin{2};
CORR = varargin{3};
Work = varargin{4};
if length(varargin) == 5
    linReg = varargin{5};
end
if length(varargin) == 6
    linReg = varargin{5};
    fitted0 = varargin{6};
end


if ~isempty(hFig)
    set(0,'currentFigure',hFig);
end

if strcmpi(CORR.presOrDens, 'pressure')
    levelLabel = 'Pres (dbar)';
elseif strcmpi(CORR.presOrDens, 'density')
    levelLabel = 'Density';
end

% =========================================================================
%% Plots
% =========================================================================

switch plotTyp
    case 0
        % -----------------------------------------------------------------
        % DOXY or PSAT from the reference data VS from float : choose the
        % kind of linear regression
        % -----------------------------------------------------------------
        if isempty(hFig)
            hFig = figure('unit','normalized','OuterPosition',[0 0 1 1],...
                'Name',sprintf('CORRECTION DOXY and PSAT - %d',Work.wmo),'NumberTitle','off');
        end
        
        hold on;
        grid on;
        plot(CORR.(CORR.strField),CORR.([CORR.strField 'ref']),'.k');
        plotColor = {'b','r','m'};
        strColor = {'blue','red','magenta'};
        sizeText = 0.9;
        xlabel('PSAT from float');
        ylabel('PSAT REF');
        if fitted0 == 0
            minFit = nanmin(nanmin(CORR.(CORR.strField)));
        else
            minFit = 0;
        end
        for i = 1:length(linReg.intercept)
            plot([minFit nanmax(nanmax(CORR.(CORR.strField)))],CORR.LRfit(i,:),'color',plotColor{i});
            if linReg.intercept(i) == false
                cmpl = 'free z intercept';
            else
                cmpl = 'zero intercept';
            end
            fprintf('\t %s, in %s %s => a=%2.3f, b=%2.3f, R²=%2.3f\n',cmpl,strColor{i},...
                linReg.equation_for_plot{i},CORR.LRcoef(i,1),CORR.LRcoef(i,2),CORR.LRcoef(i,3));
            text(0.05,sizeText,sprintf('%s => a=%2.3f, b=%2.3f, R²=%2.3f\n',...
                linReg.equation_for_plot{i},CORR.LRcoef(i,1),CORR.LRcoef(i,2),CORR.LRcoef(i,3)),'color',plotColor{i},'units','normalized','fontsize',Work.fontsize,'interp','none')
            sizeText = sizeText-0.05;
        end
        hold off;
        drawnow;
        
        
    case 1
        % -----------------------------------------------------------------
        % PSAT from the reference data VS from float (resp. O2)
        % -----------------------------------------------------------------
        if isempty(hFig)
            hFig = figure('unit','normalized','OuterPosition',[0 0 1 1],...
                'Name',sprintf('CORRECTION DOXY and PSAT - %d',Work.wmo),'NumberTitle','off');
        else
            subplot(3,3,1);
        end
        
        hold on;
        grid on;
        plot(CORR.([CORR.strField '_ko']),CORR.([CORR.strField 'ref_ko']),'ok','MarkerSize', 5,'MarkerFaceColor','w');
        plot(CORR.([CORR.strField '_ok']),CORR.([CORR.strField 'ref_ok']),'ok','MarkerSize', 5,'MarkerFaceColor','k');
        plot([nanmin(nanmin(CORR.(CORR.strField))) nanmax(nanmax(CORR.(CORR.strField)))],CORR.LRfit(1,:),'b');
        plot([nanmin(nanmin(CORR.([CORR.strField '_ok']))) nanmax(nanmax(CORR.([CORR.strField '_ok'])))],CORR.LRfit(2,:),'r');
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        xlabel(sprintf('%s from float',CORR.whichO2quantity));
        ylabel(sprintf('%s %s',CORR.whichO2quantity, CORR.whichCorr));
        %fprintf('all data, in blue %s => a=%2.3f, b=%2.3f, R²=%2.3f\n',CORR.linReg,CORR.LRcoef(1,1),CORR.LRcoef(1,2),CORR.LRcoef(1,3));
        %fprintf('grad<0.2, in red %s => a=%2.3f, b=%2.3f, R²=%2.3f\n',CORR.linReg,CORR.LRcoef(2,1),CORR.LRcoef(2,2),CORR.LRcoef(2,3));
        text(0.05,0.9,sprintf('all: %s => a=%2.3f, b=%2.3f, R²=%2.3f\n',...
            CORR.linReg,CORR.LRcoef(1,1),CORR.LRcoef(1,2),CORR.LRcoef(1,3)),'color','b','units','normalized','fontsize',Work.fontsize)
        text(0.05,0.85,sprintf('grad<0.2: %s => a=%2.3f, b=%2.3f, R²=%2.3f\n',...
            CORR.linReg,CORR.LRcoef(2,1),CORR.LRcoef(2,2),CORR.LRcoef(2,3)),'color','r','units','normalized','fontsize',Work.fontsize)
        hold off;
        drawnow;
        
    case 2
        % -----------------------------------------------------------------
        % PSAT reference - PSAT float VS PSAT float (resp. O2)
        % -----------------------------------------------------------------
        if isempty(hFig)
            figure('unit','normalized','Position',[0 0 0.48 0.55]);
        else
            subplot(3,3,2);
        end
        hold on;
        grid on;
        plot(CORR.(CORR.strField), CORR.Diff,'k*')
        plot(CORR.(CORR.strField), CORR.Diff - CORR.Diff + ...
            CORR.([CORR.strField '_avg']) + 2.8 * CORR.([CORR.strField '_std']),'+m')
        plot(CORR.(CORR.strField), CORR.Diff - CORR.Diff + ...
            CORR.([CORR.strField '_avg']) - 2.8 * CORR.([CORR.strField '_std']),'+m')
        hold off
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        xlabel(sprintf('%s from float',CORR.whichO2quantity));
        ylabel(sprintf('%s %s - %s float',CORR.whichO2quantity, CORR.whichCorr, CORR.whichO2quantity));
        title(sprintf(' %d - Correction applied on %s',Work.wmo,CORR.whichO2quantity));
        drawnow
        
    case 3
        % -----------------------------------------------------------------
        % Final PSAT From ref VS PSAT from float (resp. O2)
        % -----------------------------------------------------------------
        % isnok : cycles for which not to plot
        isnok = varargin{5};
        myVar = CORR.([CORR.strField '_ok'])(~isnok);
        
        if isempty(hFig)
            figure('unit','normalized','Position',[0.71 0.37 0.29 0.33]);
        else
            subplot(3,3,3);
        end
        hold on;
        grid on;
        plot(CORR.([CORR.strField '_ko']),CORR.([CORR.strField 'ref_ko']),...
            'ok','MarkerSize', 5,'MarkerFaceColor','w')                % rejected by dO2/dz %marine 19/06/19
        plot(CORR.([CORR.strField '_ok'])(isnok),CORR.([CORR.strField 'ref_ok'])(isnok),...
            'ok','MarkerSize', 5,'MarkerFaceColor','y')                % rejected by avg+std
        plot(CORR.([CORR.strField '_ok'])(~isnok),CORR.([CORR.strField 'ref_ok'])(~isnok),...
            'ok','MarkerSize', 5,'MarkerFaceColor','k')                % the good ones
        plot([nanmin(nanmin(myVar)) nanmax(nanmax(myVar))],CORR.LRfit(3,:),'r','LineWidth',2);
        set(gca,'fontweight','bold','fontsize',Work.fontsize)
        fprintf('\t aX+b <0.2, avg ± 2.8 std => a=%2.3f, b=%2.3f, R²=%2.3f\n',...
            CORR.LRcoef(3,1),CORR.LRcoef(3,2),CORR.LRcoef(3,3));
        text(0.05,0.9,'aX+b<0.2, avg+-2.8 std =>','color','r','units','normalized','fontsize',Work.fontsize)
        text(0.05,0.8,sprintf('a=%2.3f, b=%2.3f, R²=%2.3f\n',...
            CORR.LRcoef(3,1),CORR.LRcoef(3,2),CORR.LRcoef(3,3)),'color','r','units','normalized','fontsize',Work.fontsize)
        hold off
        xlabel(sprintf('%s from float',CORR.whichO2quantity));
        ylabel(sprintf('%s %s', CORR.whichO2quantity, CORR.whichCorr));
        drawnow
        
    case 4
        % -----------------------------------------------------------------
        % PSAT VS number of day since launch (resp. O2), for REF (blue),
        % float (black), and float corrected (red)
        % -----------------------------------------------------------------
        dayjul = varargin{5};
        
        if nargin == 6
            subplotPlace = varargin{6};
        else
            subplotPlace = 4;
        end
        if isempty(hFig)
            figure('unit','normalized','Position',[0.29 0 0.29 0.33])
        else
            subplot(3,3,subplotPlace);
        end
        hold on;
        grid on;
        
        isPos=and(Work.PSAT(:,1:3)>20,Work.PSAT(:,1:3)<170); %marine 26/06/19  %marine 18/06/19
        NS_data_raw=Work.PSAT(:,1:3);
        NS_data_raw(~isPos)=NaN;
        h1 = plot(dayjul,NS_data_raw,'-*k');
        isPos=and(CORR.psat_corr(:,1:3)>50,CORR.psat_corr(:,1:3)<150); %marine 11/06/19  %marine 18/06/19
        NS_datacorr=CORR.psat_corr(:,1:3);
        NS_datacorr(~isPos)=NaN;
        h2 = plot(dayjul,NS_datacorr,'-*r');
        
        if isfield(Work,['PSAT_WOA' CORR.cmpl])
            h3 = plot(dayjul,...
                Work.(['PSAT_WOA' CORR.cmpl])(:,1:3),'-*b');
        end
        
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        abscisseLim = get(gca,'XLim');
        set(gca,'XLim',[-10 abscisseLim(2)]);
        xlabel('Julian day from deployement','fontsize',Work.fontsize);
        ylabel('PSAT (%)','fontsize',Work.fontsize);
        legend;
        hDGroup = hggroup;
        hSGroup = hggroup;
        hCGroup = hggroup;
        set(h1,'Parent',hDGroup);
        set(h2,'Parent',hSGroup);
        if exist('h3','var')
            set(h3,'Parent',hCGroup);
        end
        set(get(get(hDGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        if exist('h3')
            set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
        end
        if exist('h3')
            legendstr = {sprintf('%s float','PSAT'), ...
                sprintf('%s float corrected','PSAT'), ...
                sprintf('%s %s','PSAT', 'WOA')};
        else
            legendstr = {sprintf('%s float','PSAT'), ...
                sprintf('%s float corrected','PSAT')};
        end
        legend(legendstr,'location','SouthEast');
        hold off;
        drawnow;
        
    case 5
        % -----------------------------------------------------------------
        % PSAT VS PRES, for REF (blue), float (black), and float corrected
        % (red).
        % -----------------------------------------------------------------
        if nargin == 5
            subplotPlace = varargin{5};
        else
            subplotPlace = 5;
        end
        if isempty(hFig)
            figure('unit','normalized','Position',[0.58 0 0.29 0.33])
        else
            subplot(3,3,subplotPlace);
        end
        
        hold on;
        grid on;
        h1=hggroup;
        plot(Work.PSAT, CORR.level,'k*','Parent',h1);
        
        %Plotting only data with good QC
        nbr_prof = size(Work.PSAT,1); %marine 11/06/19
        QC=[1,2,3];
        h2=hggroup;
        for i_plot = 1:nbr_prof
            for k=1:length(QC)
                isok = strfind(Work.DOXY_QC(i_plot,:),num2str(QC(k)));
                plot(Work.(['PSAT_LINCORR_' CORR.whichO2quantity])(i_plot,isok), CORR.level(i_plot,isok),'r*','Parent',h2)
            end
        end
        
        if ismember(Work.whichCorr,{'WOA'})
            refField = Work.whichCorr;
            ind_ref=CORR.refCyc;
        elseif ismember(Work.whichCorr,{'REF'})
            refField = Work.whichCorr;
            if (CORR.nbprof+1)==size(CORR.level,1) && CORR.refCyc~=1
                ind_ref=CORR.refCyc+1;
            else
                ind_ref=CORR.refCyc;
            end
        elseif ismember(Work.whichCorr,{'INAIR'})
            refField = 'WOA';
            ind_ref=CORR.refCyc;
        end
        
        h3=hggroup;
        plot(Work.(['PSAT_' refField CORR.cmpl])(1:length(ind_ref),:), CORR.level(ind_ref,:),'b*','Parent',h3);
        set(gca,'YDir','reverse');
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        xlabel('PSAT (%)');
        ylabel(levelLabel);
        legend;
        set(get(get(h1,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(h2,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(h3,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        legend(sprintf('%s float','PSAT'), ...
            sprintf('%s float corrected','PSAT'), ...
            sprintf('%s %s','PSAT', refField),'location','SouthEast');
        drawnow;
        
    case 6
        % -----------------------------------------------------------------
        % O2 VS PRES, for REF (blue), float (black), and float corrected
        % (red)
        % -----------------------------------------------------------------
        if nargin == 5
            subplotPlace = varargin{5};
        else
            subplotPlace = 6;
        end
        if isempty(hFig)
            figure('unit','normalized','Position',[0.71 0 0.29 0.33]);
        else
            subplot(3,3,subplotPlace);
        end
        hold on;
        grid on;
        h1=hggroup;
        plot(Work.DOXY, CORR.level,'k*','parent',h1);
        
        %Plotting only data with good QC
        nbr_prof = size(CORR.([CORR.strField]),1);
        QC=[1,2,3];
        h2=hggroup;
        for i_plot = 1:nbr_prof
            for k=1:length(QC)
                isok = strfind(Work.DOXY_QC(i_plot,:),num2str(QC(k)));
                plot(Work.(['DOXY_LINCORR_' CORR.whichO2quantity])(i_plot,isok), CORR.level(i_plot,isok),'r*','parent',h2)
                
            end
        end
        
        if ismember(Work.whichCorr,{'WOA'})
            refField = Work.whichCorr;
            ind_ref=CORR.refCyc;
        elseif ismember(Work.whichCorr,{'REF'})
            refField = Work.whichCorr;
            if (CORR.nbprof+1)==size(CORR.level,1) && CORR.refCyc~=1
                ind_ref=CORR.refCyc+1;
            else
                ind_ref=CORR.refCyc;
            end
        elseif ismember(Work.whichCorr,{'INAIR'})
            refField = 'WOA';
            ind_ref=CORR.refCyc;
        end
        h3=hggroup;
        plot(Work.(['DOXY_' refField CORR.cmpl])(1:length(ind_ref),:), CORR.level(ind_ref,:),'b*','parent',h3);
        set(gca,'YDir','reverse');
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        xlabel(sprintf('DOXY (%s)',Work.unit));
        ylabel(levelLabel);
        legend
        
        set(get(get(h1,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(h2,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(h3,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        legend('DOXY float',sprintf('DOXY float corrected (%s)',CORR.whichCorr),...
            sprintf('DOXY %s',refField),...
            'location','NorthEast');
        drawnow;
        
    case 7
        % -----------------------------------------------------------------
        % Constant correction
        % -----------------------------------------------------------------
        argoWork = varargin{5};
        
        if isempty(hFig)
            figure('unit','normalized','Position',[0 0.56 0.29 0.33]);
        else
            subplot(3,3,7);
        end
        hold on;
        grid on;
        
        if ~isstruct(argoWork)
            plot(0,0);
            text(0,0,{sprintf('R² = %2.3f',CORR.LRcoef(3,3));'=> No constant correction'},'fontsize',Work.fontsize);
            set(gca,'fontweight','bold','fontsize',Work.fontsize);
            drawnow;
            return
        end
        
        if strcmp(CORR.whichO2quantity,'PSAT')
            raw = argoWork.(CORR.strField).data;
        elseif strcmp(CORR.whichO2quantity,'DOXY')
            raw = Work.DOXY;
        end
        
        plot(nanmean(raw,1),nanmean(CORR.level,1),'k--','LineWidth',1);
        plot(nanmean(Work.([CORR.whichO2quantity '_OFFSETCORR_' upper(CORR.strField)]),1),nanmean(CORR.level,1),'g-','LineWidth',2);
        plot(nanmean(Work.(sprintf('%s_%s%s',upper(CORR.strField),CORR.whichCorr,CORR.cmpl)),1),nanmean(CORR.level(CORR.refCyc,:),1),'b-','LineWidth',2);
        hold off;
        set(gca,'YDir','reverse');
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        xlabel('PSAT (%)');
        ylabel(levelLabel);
        set(gca,'fontweight','bold','fontsize',Work.fontsize)
        legend(sprintf('averaged %s float',CORR.whichO2quantity), ...
            sprintf('averaged %s float offset',CORR.whichO2quantity), ...
            sprintf('averaged %s %s',CORR.whichO2quantity, CORR.whichCorr),'location','NorthEast');
        drawnow;
        
    case 8
        % -----------------------------------------------------------------
        % Control plot
        % -----------------------------------------------------------------
        argoWork = varargin{5};
        if isempty(hFig)
            figure('unit','normalized','Position',[0.29 0.56 0.29 0.33]);
        else
            subplot(3,3,8);
        end
        hold on;
        grid on;
        
        if ~isstruct(argoWork)
            plot(0,0);
            text(0,0,{sprintf('R² = %2.3f',CORR.LRcoef(3,3));'=> No constant correction'},'fontsize',Work.fontsize);
            set(gca,'fontweight','bold','fontsize',Work.fontsize);
            drawnow;
            
            if Work.savePlot
                if Work.(['FIT_intercept_' Work.whichO2quantity]) == 0, offset = 'nooffset';
                else, offset = 'okoffset';
                end
                if Work.presEff, presEff = 'okpreseff'; else, presEff = 'nopreseff'; end
                if Work.DODRIFT
                    saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_corr_%d_%s_%s_okdrift_on%s_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,Work.whichDrift,offset,presEff));
                else
                    saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_corr_%d_%s_%s_nodrift_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,offset,presEff));
                end
                [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);
                
            end
            return
        end
        h1 = plot(Work.PSAT, CORR.level,'k*');
        h2 = plot(Work.(sprintf('PSAT_LINCORR_%s',CORR.whichO2quantity)),CORR.level,'r*');
        h3 = plot(Work.(sprintf('PSAT_OFFSETCORR_%s',CORR.whichO2quantity)),CORR.level,'g*');
        h4 = plot(Work.(['PSAT_' CORR.whichCorr CORR.cmpl]),CORR.level(CORR.refCyc,:),'b*');
        set(gca,'YDir','reverse');
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        xlabel('PSAT (%)');
        ylabel(levelLabel);
        title(sprintf('%d',Work.wmo))
        legend
        hDGroup = hggroup;
        hSGroup = hggroup;
        hCGroup = hggroup;
        hFGroup = hggroup;
        set(h1,'Parent',hDGroup)
        set(h2,'Parent',hSGroup)
        set(h3,'Parent',hCGroup)
        set(h4,'Parent',hFGroup)
        set(get(get(hDGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(hFGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        legend('PSAT float',sprintf('PSAT float corrected (%s)',CORR.whichCorr),...
            'PSAT float offset',sprintf('PSAT %s',CORR.whichCorr),...
            'location','NorthEast');
        hold off;
        drawnow;
        
    case 9
        % -----------------------------------------------------------------
        % Control plot
        % -----------------------------------------------------------------
        argoWork = varargin{5};
        if isempty(hFig)
            figure('unit','normalized','Position',[0.29 0.56 0.29 0.33]);
        else
            subplot(3,3,9);
        end
        hold on;
        grid on;
        
        if ~isstruct(argoWork)
            plot(0,0);
            text(0,0,{sprintf('R² = %2.3f',CORR.LRcoef(3,3));'=> No constant correction'},'fontsize',Work.fontsize);
            set(gca,'fontweight','bold','fontsize',Work.fontsize);
            drawnow;
            
            if Work.savePlot
                if Work.(['FIT_intercept_' Work.whichO2quantity]) == 0, offset = 'nooffset';
                else, offset = 'okoffset';
                end
                if Work.presEff, presEff = 'okpreseff'; else, presEff = 'nopreseff'; end
                if Work.DODRIFT
                    if strcmp(Work.whichCorr,'REF') || strcmp(Work.whichCorr,'WOA')
                        %whichDrift='WOA'; % Commented by Thierry Reynaud
                        %13/03/2020
                        whichDrift=Work.whichDrift;
                    end
                    
                    saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_corr_%d_%s_%s_okdrift_on%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,Work.whichDrift,offset,presEff));
                else
                    saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_corr_%d_%s_%s_nodrift_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,offset,presEff));
                end
                [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);
                
            end
            return
        end
        h1 = plot(Work.DOXY,CORR.level,'k*');
        h2 = plot(Work.(sprintf('DOXY_LINCORR_%s',CORR.whichO2quantity)),CORR.level,'r*');
        h3 = plot(Work.(sprintf('DOXY_OFFSETCORR_%s',CORR.whichO2quantity)),CORR.level,'g*');
        h4 = plot(Work.(['DOXY_' CORR.whichCorr CORR.cmpl]),CORR.level(CORR.refCyc,:),'b*');
        set(gca,'YDir','reverse');
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        xlabel(sprintf('DOXY (%s)',Work.unit));
        ylabel(levelLabel);
        legend
        hDGroup = hggroup;
        hSGroup = hggroup;
        hCGroup = hggroup;
        hFGroup = hggroup;
        set(h1,'Parent',hDGroup)
        set(h2,'Parent',hSGroup)
        set(h3,'Parent',hCGroup)
        set(h4,'Parent',hFGroup)
        set(get(get(hDGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        set(get(get(hFGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        legend('DOXY float',sprintf('DOXY float corrected (%s)',CORR.whichCorr),...
            'DOXY float offset',sprintf('DOXY %s',CORR.whichCorr),...
            'location','NorthEast');
        hold off;
        drawnow;
        
    case 10 %--MG
        % -----------------------------------------------------------------
        % DATA corrected with reference profile
        % -----------------------------------------------------------------
        
        if ismember(Work.whichCorr,{'WOA'})
            refField = Work.whichCorr;
        elseif ismember(Work.whichCorr,{'REF'})
            refField = Work.whichCorr;
        elseif ismember(Work.whichCorr,{'INAIR'})
            refField = 'WOA';
        end
        
        %For constante correction (only with WOA or REF correction)
        if isfield(Work,['OFFSET_' CORR.whichO2quantity])
            set(0,'CurrentFigure',hFig);
            hold on;
            nbr_prof = size(CORR.psat,1);
            map = colormap(jet(nbr_prof));
            QC=[1,2,3];
            hg_mar=hggroup;
            for i_mar = 1:nbr_prof
                if i_mar~=2
                    for k=1:length(QC)
                        isok = strfind(Work.DOXY_QC(i_mar,:),num2str(QC(k)));
                        plot(Work.(['DOXY_OFFSETCORR_' CORR.whichO2quantity])(i_mar,isok),CORR.level(i_mar,isok),'Color',map(i_mar,:),'Parent',hg_mar);
                    end
                end
            end
            hg_mar_r=hggroup;
            if strcmp(Work.whichCorr,'REF')
                %Ignore 1D profile
                if (CORR.nbprof+1)==size(CORR.level,1) && CORR.refCyc~=1
                    ind_ref=CORR.refCyc+1;
                else
                    ind_ref=CORR.refCyc;
                end
                plot(Work.(['DOXY_LINCORR_' CORR.whichO2quantity])(ind_ref,:), CORR.level(ind_ref,:),'r','LineWidth',3,'Parent',hg_mar_r);
                hg_mar_c=hggroup;
                plot(Work.(['DOXY_' refField CORR.cmpl]), CORR.level(ind_ref,:),'k','LineWidth',3,'Parent',hg_mar_c);
                ref_leg=CORR.refCyc;
            else
                %Ignore 1D profile
                if Work.plotREF.refCycle~=1
                    ind_ref=Work.plotREF.refCycle+1;
                else
                    ind_ref=Work.plotREF.refCycle;
                end
                plot(Work.(['DOXY_LINCORR_' CORR.whichO2quantity])(ind_ref,:), CORR.level(ind_ref,:),'r','LineWidth',3,'Parent',hg_mar_r);
                hg_mar_c=hggroup;
                plot(Work.plotREF.(['DOXY_REF' CORR.cmpl]), CORR.level(ind_ref,:),'k','LineWidth',3,'Parent',hg_mar_c);
                ref_leg=Work.plotREF.refCycle;
            end
            
            grid on;
            legend;
            set(get(get(hg_mar,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
            set(get(get(hg_mar_r,'Annotation'),'LegendInformation'),...hD
                'IconDisplayStyle','on');
            set(get(get(hg_mar_c,'Annotation'),'LegendInformation'),...hD
                'IconDisplayStyle','on');
            
            legend('O2 float',sprintf('O2 float, cycle %d',ref_leg),'O2 REF','Location','Southwest');
            set(gca,'YDir','reverse');
            set(gca,'fontweight','bold');
            xlabel('[O2] (umol/kg)','fontweight','bold','fontsize',Work.fontsize);
            ylabel('Pres (dbar)','fontweight','bold','fontsize',Work.fontsize);
            title(sprintf('CORRECTED DATA \n %d',Work.wmo),...
                'fontweight','bold','fontsize',Work.fontsize);
            ylim([0,nanmax(nanmax(Work.PRES))]);
            hold off;
            drawnow;
            
            %For normal correction
        else
            set(0,'CurrentFigure',hFig);
            hold on;
            nbr_prof = size(Work.PSAT,1); %marine 11/06/19
            map = colormap(jet(nbr_prof));
            QC=[1,2,3];
            hg_mar=hggroup;
            for i_mar = 1:nbr_prof
                if i_mar~=2
                    for k=1:length(QC)
                        isok = strfind(Work.DOXY_QC(i_mar,:),num2str(QC(k)));
                        plot(Work.(['DOXY_LINCORR_' CORR.whichO2quantity])(i_mar,isok),CORR.level(i_mar,isok),'Color',map(i_mar,:),'Parent',hg_mar);
                    end
                end
            end
            hg_mar_r=hggroup;
            
            if strcmp(Work.whichCorr,'REF')
                %Ignore 1D profile
                if (CORR.nbprof+1)==size(CORR.level,1) && CORR.refCyc~=1
                    ind_ref=CORR.refCyc+1;
                else
                    ind_ref=CORR.refCyc;
                end
                plot(Work.(['DOXY_LINCORR_' CORR.whichO2quantity])(ind_ref,:), CORR.level(ind_ref,:),'r','LineWidth',3,'Parent',hg_mar_r);
                hg_mar_c=hggroup;
                plot(Work.(['DOXY_' refField CORR.cmpl]), CORR.level(ind_ref,:),'k','LineWidth',3,'Parent',hg_mar_c);
                ref_leg=CORR.refCyc;
            else
                %Ignore 1D profile
                if exist('Work.plotREF.refCycle)')  % Added byTR 05.04.2020
                    if Work.plotREF.refCycle~=1
                        ind_ref=Work.plotREF.refCycle+1;
                    else
                        ind_ref=Work.plotREF.refCycle;
                    end
                    plot(Work.(['DOXY_LINCORR_' CORR.whichO2quantity])(ind_ref,:), CORR.level(ind_ref,:),'r','LineWidth',3,'Parent',hg_mar_r);
                    hg_mar_c=hggroup;
                    plot(Work.plotREF.(['DOXY_REF' CORR.cmpl]), CORR.level(ind_ref,:),'k','LineWidth',3,'Parent',hg_mar_c);
                    ref_leg=Work.plotREF.refCycle;
                end
            end
            
            grid on;
            legend;
            if exist('hg_mar')
            set(get(get(hg_mar,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
            end
            if exist('hg_mar_r')
            set(get(get(hg_mar_r,'Annotation'),'LegendInformation'),...hD
                'IconDisplayStyle','on');
            end
            if exist('hg_mar_c')
            set(get(get(hg_mar_c,'Annotation'),'LegendInformation'),...hD
                'IconDisplayStyle','on');
            end
            
            % Modified by. TR 05.04.2020
            if exist('ref_leg')
                legend('O2 float',sprintf('O2 float, cycle %d',ref_leg),'O2 REF','Location','Southwest');
            else
                ref_leg=NaN;
                legend('O2 float',sprintf('O2 float, cycle %d',ref_leg),'O2 REF','Location','Southwest');  
            end
            set(gca,'YDir','reverse');
            set(gca,'fontweight','bold');
            xlabel('[O2] (umol/kg)','fontweight','bold','fontsize',Work.fontsize);
            ylabel('Pres (dbar)','fontweight','bold','fontsize',Work.fontsize);
            title(sprintf('CORRECTED DATA \n %d',Work.wmo),...
                'fontweight','bold','fontsize',Work.fontsize);
            ylim([0,nanmax(nanmax(Work.PRES))]);
            hold off;
            drawnow;
        end
        
        
        
    case 11
        % -----------------------------------------------------------------
        % PPOX timeSerie measurement, inflated, deflated and atmospheric
        % -----------------------------------------------------------------
        if isempty(hFig)
            hFig = figure('unit','normalized','OuterPosition',[0 0 1 1],...
                'Name',sprintf('CORRECTION DOXY and PSAT - %d',Work.wmo),'NumberTitle','off');
        else
            subplot(3,3,[1 2]);
        end
        
        largeMarker = 10;
        hold on;
        grid on;
        
        % -----------------------------------------------------------------
        %--MG Calculation of average value per cycle (figure more readable)
        % -----------------------------------------------------------------
        k1_air=1;
        k2_air=0;
        k1_water=1;
        k2_water=0;
        for i=1:length(Work.nbrMeas.inAir)
            k2_air=Work.nbrMeas.inAir(i);
            k2_water=Work.nbrMeas.inWater(i);
            %inwater
            po2_inwater(i)=mean(CORR.po2_inwater(k1_water:k1_water+k2_water-1));
            std_inwater(i)=nanstd(CORR.po2_inwater(k1_water:k1_water+k2_water-1));
            po2_time_water(i)=mean(CORR.po2_inwater_gtime(k1_water:k1_water+k2_water-1));
            %inair
            isokDim = size(CORR.po2_inair) == length(CORR.gtime);
            gtime_for_plot = repmat(CORR.gtime,size(CORR.po2_inair,find(~isokDim)),1);
            po2_time_air(i)=mean(gtime_for_plot(k1_air:k1_air+k2_air-1));
            po2_inair(i)=mean(CORR.po2_inair(k1_air:k1_air+k2_air-1));
            std_inair(i)=nanstd(CORR.po2_inair(k1_air:k1_air+k2_air-1));
            ncep_po2(i)=mean(CORR.ncep_pO2_scaled(k1_air:k1_air+k2_air-1));
            ncep_time(i)=mean(CORR.gtime(k1_air:k1_air+k2_air-1));
            k1_air=k1_air+k2_air;
            k1_water=k1_water+k2_water;
        end
        
        h1 = errorbar(po2_time_water,po2_inwater,std_inwater,'vb','MarkerSize',5);
        h2 = errorbar(po2_time_air,po2_inair,std_inair,'^r','MarkerSize',5);
        h3 = plot(ncep_time,ncep_po2,'.k','MarkerSize',8);
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        ylabel('\itp\rmO_2 / mbar')
        lh = legend([h1 h2(1),h3],{'argo water';'argo air';'ncep'});
        set(lh,'box','off','orientation','horizontal','location','SouthOutside')
        title(sprintf('In-Air measurements time series, WMO = %d',Work.wmo))
        dynamicDateTicks(gca);
        set(gca,'XTickLabelRotation',45);
        
    case 13
        % -----------------------------------------------------------------
        % PPOX timeSerie measurement, inflated, deflated and atmospheric
        % -----------------------------------------------------------------
        if isempty(hFig)
            hFig = figure('unit','normalized','OuterPosition',[0 0 1 1],...
                'Name',sprintf('CORRECTION DOXY and PSAT - %d',Work.wmo),'NumberTitle','off');
        else
            subplot(3,3,3);
        end
        
        hold on;
        grid on;
        plot(CORR.po2_inwater_sizeInair - CORR.ncep_pO2_scaled,CORR.po2_inair - CORR.ncep_pO2_scaled,'r.');
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        axis auto
        xlabel('argo water - ncep (mbar)');
        ylabel('argo air - ncep (mbar)');
        title({'Argo in air vs Argo in water';'referenced to NCEP values'})
        textpos = [];
        
        if Work.isokC==1 %--MG
            plotleg = ['c = ' num2str(CORR.NLRcoef(2,1),'%.1f') ...
                '\pm' num2str(CORR.NLRcoef(2,2),'%.1f') ' mbar' char(10) ...
                'm = ' num2str(1+CORR.NLRcoef(1,1)/100,'%.1f')...
                '\pm' num2str(CORR.NLRcoef(1,2),'%.1f') '%' char(10) ...
                'RMSE = ' num2str(CORR.NLRcoef(1,3),'%.2f') ' mbar'];
        elseif Work.isokC==0
            plotleg = ['m = ' num2str(1+CORR.NLRcoef(1)/100,'%.1f')...
                '\pm' num2str(CORR.NLRcoef(2),'%.1f') '%' char(10) ...
                'RMSE = ' num2str(CORR.NLRcoef(3),'%.2f') ' mbar'];
        end
        
        legend(plotleg,'location','best');
        clear plotleg
        hold off;
        
    case 14
        % -----------------------------------------------------------------
        % PPOX timeSerie measurement, inflated, deflated and atmospheric
        % -----------------------------------------------------------------
        if isempty(hFig)
            hFig = figure('unit','normalized','OuterPosition',[0 0 1 1],...
                'Name',sprintf('CORRECTION DOXY and PSAT - %d',Work.wmo),'NumberTitle','off');
        else
            subplot(3,3,[4:5]);
        end
        largeMarker = 10;
        hold on;
        grid on;
        
        % -----------------------------------------------------------------
        %--MG Calculation of average value per cycle (figure more readable)
        % -----------------------------------------------------------------
        k1_air=1;
        k2_air=0;
        k1_water=1;
        k2_water=0;
        for i=1:length(Work.nbrMeas.inAir)
            k2_air=Work.nbrMeas.inAir(i);
            k2_water=Work.nbrMeas.inWater(i);
            %inwater
            po2_inwater(i)=mean(CORR.po2_inwater(k1_water:k1_water+k2_water-1));
            std_inwater(i)=nanstd(CORR.po2_inwater(k1_water:k1_water+k2_water-1));
            po2_time_water(i)=mean(CORR.po2_inwater_gtime(k1_water:k1_water+k2_water-1));
            %inair
            isokDim = size(CORR.po2_inair) == length(CORR.gtime);
            gtime_for_plot = repmat(CORR.gtime,size(CORR.po2_inair,find(~isokDim)),1);
            po2_time_air(i)=mean(gtime_for_plot(k1_air:k1_air+k2_air-1));
            po2_inair(i)=mean(CORR.po2_inair(k1_air:k1_air+k2_air-1));
            std_inair(i)=nanstd(CORR.po2_inair(k1_air:k1_air+k2_air-1));
            %ncep
            ncep_time(i)=mean(CORR.gtime(k1_air:k1_air+k2_air-1));
            ncep_po2(i)=mean(CORR.ncep_pO2_scaled(k1_air:k1_air+k2_air-1));
            k1_air=k1_air+k2_air;
            k1_water=k1_water+k2_water;
        end
        
        h1 = plot(po2_time_water,po2_inwater,'vb','MarkerSize',5);
        h2 = plot(po2_time_air,po2_inair,'^r','MarkerSize',5);
        h4 = plot(po2_time_water,po2_inwater*(1+CORR.NLRcoef(1)/100),'vc', 'MarkerFaceColor','c','MarkerSize',5);
        h5 = plot(po2_time_air,po2_inair*(1+CORR.NLRcoef(1)/100),'^m','MarkerSize',5,'MarkerFaceColor','m');
        h3 = plot(ncep_time,ncep_po2,'.k','MarkerSize',8);
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        %  set(h3,'markersize',largeMarker)
        ylabel('\itp\rmO_2 / mbar')
        lh = legend([h1,h2(1),h3,h4,h5(1)],{'argo water';'argo air';'ncep';'argo water corrected';'argo air corrected'});
        set(lh,'box','off','orientation','horizontal','location','SouthOutside')
        title(sprintf('In-Air measurements time series, WMO = %d',Work.wmo))
        dynamicDateTicks(gca);
        set(gca,'XTickLabelRotation',45);
        
    case 16
        % -----------------------------------------------------------------
        % PPOX timeSerie measurement, inflated, deflated and atmospheric
        % -----------------------------------------------------------------
        if isempty(hFig)
            hFig = figure('unit','normalized','OuterPosition',[0 0 1 1],...
                'Name',sprintf('CORRECTION DOXY and PSAT - %d',Work.wmo),'NumberTitle','off');
        else
            subplot(3,3,6);
        end
        
        %--MG
        if Work.isokC==1
            infoStat = {'Equation:';...
                sprintf('%s\n',CORR.NLRstreq);...
                'Coefficients:';...
                sprintf('m= %2.2f +/- %2.2f %%', 1+CORR.NLRcoef(1,1)/100, CORR.NLRcoef(1,2));...
                sprintf('c = %2.2f +/- %2.2f mBar\n', CORR.NLRcoef(2,1), CORR.NLRcoef(2,2));...
                'Stats:';...
                sprintf('pValues = %2.2e ', CORR.NLRinfo.Coefficients.pValue(2));...
                sprintf('c pValues = %2.2e ', CORR.NLRinfo.Coefficients.pValue(1));...
                sprintf('Number of obs : %d',CORR.NLRinfo.NumObservations);...
                sprintf('RMSE : %2.3f',CORR.NLRinfo.RMSE);...
                sprintf('R-Squared : %2.3f',CORR.NLRinfo.Rsquared.Ordinary);...
                sprintf('Adjusted R-Squared : %2.3f',CORR.NLRinfo.Rsquared.Adjusted);...
                };
        else
            infoStat = {'Equation:';...
                sprintf('%s\n',CORR.NLRstreq);...
                'Coefficients:';...
                sprintf('m = %2.2f +/- %2.2f %%\n', 1+CORR.NLRcoef(1,1)/100, CORR.NLRcoef(1,2));...
                'Stats:';...
                sprintf('pValues = %2.2e ', CORR.NLRinfo.Coefficients.pValue(1));...
                sprintf('Number of obs : %d',CORR.NLRinfo.NumObservations);...
                sprintf('RMSE : %2.3f',CORR.NLRinfo.RMSE);...
                sprintf('R-Squared : %2.3f',CORR.NLRinfo.Rsquared.Ordinary);...
                sprintf('Adjusted R-Squared : %2.3f',CORR.NLRinfo.Rsquared.Adjusted);...
                };
        end
        
        
        plot(0,0);
        text(-0.7,0.1,infoStat,'interpreter','none');
        set(gca,'fontweight','bold','fontsize',Work.fontsize);
        set(gca,'Visible','Off');
        drawnow;
        if Work.savePlot
            if Work.(['FIT_intercept_' Work.whichO2quantity]) == 0, offset = 'nooffset';
            else, offset = 'okoffset';
            end
            if Work.presEff, presEff = 'okpreseff'; else, presEff = 'nopreseff'; end
            if Work.DODRIFT
                saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_corr_%d_%s_%s_okdrift_on%s_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,Work.whichDrift,offset,presEff));
            else
                saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_corr_%d_%s_%s_nodrift_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,offset,presEff));
            end
            [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);
            
        end
end

% =========================================================================
%% Save the plots
% =========================================================================
if Work.savePlot
    if plotTyp==10
        if Work.(['FIT_intercept_' Work.whichO2quantity]) == 0, offset = 'nooffset';
        else, offset = 'okoffset';
        end
        if Work.presEff, presEff = 'okpreseff'; else, presEff = 'nopreseff'; end
        
        if Work.DODRIFT
            saveFile2 = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_data_corr_%d_%s_%s_okdrift_on%s_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,Work.whichDrift,offset,presEff));
        else
            saveFile2 = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_data_corr_%d_%s_%s_nodrift_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,offset,presEff));
        end
        [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile2);
    else
        if Work.(['FIT_intercept_' Work.whichO2quantity]) == 0, offset = 'nooffset';
        else, offset = 'okoffset';
        end
        if Work.presEff, presEff = 'okpreseff'; else, presEff = 'nopreseff'; end
        
        if Work.DODRIFT
            saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_corr_%d_%s_%s_okdrift_on%s_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,Work.whichDrift,offset,presEff));
        else
            saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_corr_%d_%s_%s_nodrift_%s_%s',Work.wmo,Work.whichCorr,Work.whichO2quantity,offset,presEff));
        end
        [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);
    end
    
    
    
end
end

