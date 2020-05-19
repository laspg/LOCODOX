% DOXY_PLOT_PTS_QC plot QCs of the Pression/Temperature/Salinity data
%
% SYNTAX
% [hFig] = DOXY_PLOT_PTS_QC(varargin)
%
% DESCRIPTION
% DOXY_PLOT_PTS_QC plot QCs of the Pression (subplot1) 
%                                  Temperature (subplot2) 
%                                  Salinity (subplot3)
%                  as a function of cycle cumbers
%
% INPUT
% the seven first arguments are the following:
%         hFig = varargin{1};
%         cycn = varargin{2}; 
%         pres = varargin{3};
%         tab_pres_QC = str2numQC(varargin{4});
%         tab_temp_QC = str2numQC(varargin{5});
%         tab_psal_QC = str2numQC(varargin{6});
%         Work        = varargin{7};
% 
%   hFig (handle)        figure handles of the plot.
% 
%   cycn                 cycle number (vector length: nc)
% 
%   pres                 pressure (vector length: nz)
%
%   tab_pres_QC (string) tables of QC for pressure, temperature     
%   tab_temp_QC (string) and salinity data respectively
%   tab_psal_QC (string) dimension: nc x nz (char)
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
% OUTPUT
%   hFig (handle)        figure handles of the plot
%
% CALL
%   str2numQC
%
% SEE ALSO
%

% HISTORY
%   $created: 26/01/2017 $author: Anne Piron, Altran Ouest 
%   $Revision: version $Date: $author:
%

function [hFig] = DOXY_PLOT_PTS_QC(varargin)

    hFig = varargin{1};
    cycn = varargin{2};
    pres = varargin{3};
    if length(varargin) == 7
        tab_pres_QC = str2numQC(varargin{4});
        tab_temp_QC = str2numQC(varargin{5});
        tab_psal_QC = str2numQC(varargin{6});
        Work        = varargin{7};
    end

    QC = [0:1:9];
    mkz_std = 3;
    w(:,:)=[0.8 , 0.8 , 0.8 ;... % Gray      (QC = 0: no QC was performed)
            0.3 , 1   , 0.3 ;... % Green     (QC = 1: Good Data)
            1   , 1   , 0.3 ;... % Yellow    (QC = 2: Probably Good Data)
            1   , 0.7 , 0.3 ;... % Orange    (QC = 3: Bad Data, Potentially Correctable)
            1   , 0.3 , 0.3 ;... % Red       (QC = 4: Bad Data)
            0.8 , 0   , 0.8 ;... % Magenta   (QC = 5: Value Changed)
            0   , 0   , 0   ;... % Black     (QC = 6: not used)
            0   , 0   , 0   ;... % Black     (QC = 7: not used)
            0   , 0.8 , 0.8 ;... % Cyan      (QC = 8: Interpolated Value)
           0.3  , 0.3 , 1   ];   % Blue      (QC = 9: Missing Value)
    param = {'pressure';'temperature';'salinity'};
   
    % =========================================================================
    %% Plot figure
    % =========================================================================

    box on

    % Loop on param P/T/S (NB: fourth panel used to plot the colorbar only!)
    for ip = 1:3 % 
        subplot(4,1,ip)
        hold on
        if ip == 1 % PRES
            tab_QC = tab_pres_QC;
        elseif ip == 2 % TEMP
            tab_QC = tab_temp_QC;
        elseif ip >= 3 % PSAL
            tab_QC = tab_psal_QC;            
        end
        cmpt = 0;
        for ic = 1:size(tab_QC,1)
            for iq = 1:length(QC)
                [i] = find(squeeze(tab_QC(ic,:)) == QC(iq));
                if QC(iq) > 2
                    mkz = mkz_std + 2;
                else
                    mkz = mkz_std;
                end
                if ~isempty(i)
                    cmpt = cmpt + 1;
                    veccyc = NaN(length(i),1); veccyc(:) = cycn(ic);
                    hp(cmpt) = plot(veccyc,-squeeze(pres(ic,i)),'s','markerfacecolor',w(QC(iq)+1,:),...
                                                                    'markeredgecolor',w(QC(iq)+1,:),...
                                                                    'markersize',mkz);            
                end
            end
        end
        xlim([min(cycn)-1 max(cycn)+1])
        ylim([-(max(pres(:))+10) -(min(pres(:))-1)])
        if ip <= 3
            title(sprintf('%d - %s data quality',Work.wmo,param{ip}),'interpreter', 'none','fontweight','normal','fontsize',Work.fontsize)
        end
        colormap(w)
        
        if ip == 3
            xlabel('Cycle Number','fontsize',Work.fontsize)
        end
        ylabel('Depth (m)','fontsize',Work.fontsize)
    end
    
    %Add colorbar %25/06/19 marine
    subplot(4,1,4)
    hold on
    xlim([min(cycn)-1 max(cycn)+1])
    ylim([-(max(pres(:))+10) -(min(pres(:))-1)])
    cl = colorbar;
    set(cl,'Ticks',[0+0.1/2:1/10:(size(w,1)-1)/10+0.1/2],'TickLabels',...
        ['0';...
        '1';...
        '2';...
        '3';...
        '4';...
        '5';...
        '-';...
        '-';...
        '8';...
        '9'],'location','north','fontsize',Work.fontsize)
    cl.Label.String = 'QC';
    set(gca,'Visible','off')
               
    % =========================================================================
    %% Save the plot
    % =========================================================================
    if Work.savePlot
        saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_PTS_QC_%d_%s',Work.wmo,Work.whichCorr));
        [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);

    end
end


