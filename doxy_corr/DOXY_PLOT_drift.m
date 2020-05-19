% DOXY_PLOT_drift trace the doxy sensor drift
%
% SYNTAX
% [ax, Work] = DOXY_PLOT_drift(DRIFT,Work,daydiff,unit,type, ax)
%
% DESCRIPTION
% DOXY_PLOT_drift trace the argo doxy sensor drift, compared the data to
% the WOA below 1500m, where the drift should be more visible.
%
% INPUT
%   DRIFT (structure)    Drift estimating structure, from argo data and
%                        climatology data.
%                        Example:
%                        DRIFT = 
%                                      doxy: [46x85 double]
%                                   doxyref: [46x85 double]
%                                   deltaO2: [46x85 double]
%                              deltaO2_Mean: [1x85 double]
%                            doxy_corrected: [85x120 double]
%                          fitRegression: [1x85 double]
%
%   Work (struct)        float working structure, issued and computed from
%                        the argo data and the reference data (WOA or REF).
%                        Carries plot informations among other things.
%                        Example:
%                         Work = 
%                                  readme: [1x1 struct]
%                                    unit: 'mumol/kg'
%                              R2threshold: 0.8000
%                                     wmo: 1901205
%                                  sensor: 'Optode'
%                                makePlot: 1
%                                savePlot: 1
%                                 dirPlot: '/home1/homedir4/perso/ebrion/NAOS/NAOS_2015(2)/plots/REF'
%                           pres_adjusted: [1x1 struct]
%                           temp_adjusted: [1x1 struct]
%                           psal_adjusted: [1x1 struct]
%                                DOXY_WOA: [1x1 struct]
%
%   daydiff (double)      argo time array (directly from data). Julian day.
%
%   unit (string)        doxy unit. Example : 'mumol/L'
%
%   type (double)        kind of plot. To use from 1 to 3?
%
%   ax (handle)          axes handles where to plot.
%
% OUTPUT
%   Work (struct)        float working structure updated
%
%   ax (handle)          axes handles where to plot.    
%
% CALL :
%
% SEE ALSO
%   DOXY_drift, DOXY_drift_isopycn, DOXY_corr_main

% HISTORY
%   $created: 10/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%
%             26/10/2018      Marine GALLIAN, Altran Ouest :
%                             Addition of the possibility to apply a 
%                             constant drift from a certain day
%             27/03/2019      Marine GALLIAN, Altran Ouest 
%                             Improve plotting for drift computed on NCEP
%                             data
%             10.04.2020      T.Reynaud temporal ==> time


function [ax,Work] = DOXY_PLOT_drift(DRIFT,Work,daydiff,unit,type,ax)


% =========================================================================  
%% Initialisation
% =========================================================================  

switch type
    case 1
        hFig = figure('Name',sprintf('TIME DRIFT - %d',Work.wmo),'NumberTitle','off',...
                      'unit','normalized','OuterPosition',[0.67 0.50 0.33 0.50]);
        
        set(gca,'fontweight','bold');
        
        % =========================================================================  
        %% Plot interpolated doxy from argo and doxy from WOA respect to the number
        %% of day after deployement
        % =========================================================================  
        subplot(2,1,1);
        hold on;
        if iscell(DRIFT.doxy)
            doNotPlot = all(cellfun(@all,cellfun(@isnan,DRIFT.doxy,'UniformOutput',false)));
        else
            doNotPlot = all(all(isnan(DRIFT.doxy)));
        end
        
        %If do not plot
        if doNotPlot
            ax = subplot(2,1,2);
            text(0.1,0.2,sprintf('No Drift Computation: no deep data or not enough cycle'));
            return
        end
        
        % When NCEP drift, we plot the average also 
        if strcmp(Work.whichDrift,'NCEP')
            deb=0;
            nbr_meas=Work.nbrMeas.inAir;
            daydiff_tmp=zeros(size(nbr_meas));
            DRIFT_doxy_tmp=zeros(size(nbr_meas));
            for i = 1 : length(nbr_meas)  % MG
                daydiff_tmp(i)=nanmean(daydiff(deb+1:(deb+nbr_meas(i))));
                DRIFT_doxy_tmp(i)=nanmean(DRIFT.doxy(deb+1:deb+nbr_meas(i)));
                deb = deb + nbr_meas(i);
            end
            plot(daydiff_tmp, DRIFT_doxy_tmp,'o-r')                
            plot(daydiff,DRIFT.doxy,'+k');
            plot(daydiff,DRIFT.doxyref,'+-b');
        else
            if iscell(DRIFT.doxy)
                itoPlot = find(~cellfun(@isempty,DRIFT.doxy));
                for i = itoPlot
                    plot( {i},DRIFT.doxy{i},'+-k');
                    plot(daydiff{i},DRIFT.doxyref{i},'+-b');
                end               
            else
                plot(daydiff,DRIFT.doxy,'+-k');
                plot(daydiff,DRIFT.doxyref,'+-b');
            end
        end
        xlim([min(daydiff) max(daydiff)])
        xlabel('Number of days','fontweight','bold','fontsize',Work.fontsize);        
        hold off;

        
        if ~strcmp(Work.whichCorr,'INAIR')
            title(sprintf('%d - [O2] at depth below %d m (Argo in black, WOA in blue)',Work.wmo,Work.min_drift_depth));
            ylabel(sprintf('[O2] en %s',unit),'fontweight','bold','fontsize',Work.fontsize);
        else
            title(sprintf('%d - PO2 (Argo inflated in black, average per cycle in red and NCEP in blue)',Work.wmo));
            ylabel(sprintf('PO2 en %s',unit),'fontweight','bold','fontsize',Work.fontsize);
        end
        grid on;

        % =========================================================================  
        %% Plot the difference between argo doxy and woa doxy, respect to the number
        %% of day after deployement
        % =========================================================================  
        ax = subplot(2,1,2);
        hold on;
        
        % When NCEP drift, we plot the average also
        if strcmp(Work.whichDrift,'NCEP')
            
            deb=0;
            deltaO2_tmp=zeros(size(nbr_meas));
            for i = 1 : length(nbr_meas)  % MG
                deltaO2_tmp(i)=nanmean(DRIFT.deltaO2(deb+1:(deb+nbr_meas(i))));
                deb = deb + nbr_meas(i);
            end
            plot(daydiff_tmp,deltaO2_tmp ,'o-c','LineWidth',2)                
            plot(daydiff,DRIFT.deltaO2,'+g');
            
        else
            if iscell(DRIFT.doxy)
                itoPlot = find(~cellfun(@isempty,DRIFT.deltaO2));
                for i = itoPlot
                    plot(daydiff{i},DRIFT.deltaO2{i},'+-g');
                end
            else
                plot(daydiff,DRIFT.deltaO2,'+-g');
            end
        end
        xlim([min(daydiff) max(daydiff)])
        xlabel('Number of days','fontweight','bold','fontsize',Work.fontsize);
        grid on;
        if ~strcmp(Work.whichCorr,'INAIR')
            title(sprintf('%d - delta 02 = [O2]argo - [O2]woa at depth below %d m',Work.wmo,Work.min_drift_depth));
            ylabel('delta O2','fontweight','bold','fontsize',Work.fontsize);
        else
            title(sprintf('%d - delta P02 = PO2argo - PO2ncep at surface in green and average per cycle in cyan',Work.wmo));            
            ylabel('delta PO2','fontweight','bold','fontsize',Work.fontsize);
        end
        
        
    % =========================================================================  
    %% Plot the linear regression over the difference ago doxy - WOA doxy
    % =========================================================================  
    case 2
        hFig = gcf;
        hold(ax,'on');
	    plot(ax,double(daydiff),DRIFT.poliv_c,':k');        

    case 3
        hFig = gcf;
        hold(ax,'on');
        plot(ax,double(daydiff),DRIFT.fitRegression,'r');
end
set(gca,'fontweight','bold')
hold off
drawnow;

    % =========================================================================  
    %% Add a constant drift depending on the choice made
    % =========================================================================  
switch type
    case 3
        a=questdlg('Do you want to apply a constant drift from a certain day ?',sprintf('CONSTANT DRIFT - %d',Work.wmo),'Yes','No','No');
        
        while strcmp(a,'Yes') && ~isfield(Work,'ind_drift_stop')
            jour_der = cell2mat(inputdlg(sprintf('Enter the day : \n \n *It must be between 1 and %d .',round(max(daydiff))-1),sprintf('CONSTANT DRIFT - %d',Work.wmo)));
            if ~isempty(jour_der)
                for i=1:length(daydiff)
                    if daydiff(i)>str2num(jour_der) && str2num(jour_der)>-1
                        Work.ind_drift_stop=[i,str2num(jour_der)];
                        break
                    end
                end
                if ~isfield(Work,'ind_drift_stop')
                    continue
                end

                line([str2num(jour_der),str2num(jour_der)],[ax.YLim(1),ax.YLim(2)])
                text(str2num(jour_der)+1,mean(ax.YLim),'Constant drift applied from this point','FontSize',12)
            else
                fprintf('INFO >>>>>   No constant drift applied \n    ')
                break
            end
                
        end
end

% =========================================================================  
%% Save the graph
% =========================================================================  
if Work.savePlot == 1
    if Work.presEff, presEffStr = 'okpreseff'; else, presEffStr = 'nopreseff'; end    
    if strcmp(Work.whichCorr,'REF') || strcmp(Work.whichCorr,'WOA')
        drift_type='onWOA'; 
    elseif  strcmp(Work.whichCorr,'INAIR')
        drift_type='onNCEP'; 
    end
    saveFile = fullfile(Work.dirPlot,sprintf('DOXY_drift_%s_%s_%d',drift_type,presEffStr,Work.wmo));
    DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);        
end

