% DOXY_PLOT_QC plot profiles and their QC.
%
% SYNTAX
% [] = DOXY_PLOT_QC(varargin)
% [] = DOXY_PLOT_QC(hfig,plotTyp,Work,argo,argoWork)
%
% DESCRIPTION
% DOXY_PLOT_QC plot O2 data for primary, secondary or near surface sampling.
%
% INPUT
%   hFig (handle)        Figure handles of the plot.
%
%   plotTyp (double)     indicate the type of plot you want to map
%                        1: Plot the Argo Pressure VS the Argo Density
%                        2: plot the Argo doxy VS the choosen pressure,  
%                           with QC flag colored
%                        3: plot the Argo doxy VS the choosen pressure, and
%                           color the remaining points
%
%   argowork (structure) Float working structure, issued and computed
%                           from argo float data.
%                           Example:
%                           argoWork = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%
%   Work (structure)     Doxy correction working structure, issued and
%                           computed from teh float data. Modified and
%                           completed  
%
%   argo (structure)     Argo float structure (get directly from data).
%                        Example :
%                        argo =
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
%
% OUTPUT
%
% CALL :
% DOXY_PLOT_settingsToPrint
%
% SEE ALSO
%   DOXY_argo_prepare, DOXY_argo_prepare_main
%
% HISTORY
%   $created: 09/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%              v3.0  20/02/2019  Marine GALLIAN, Altran Ouest
%                                Delete the hook removal
%                                Change figure showing data and their QC



function [] = DOXY_PLOT_QC(varargin)

% =========================================================================
%% Initialisation
% =========================================================================
hFig = varargin{1};
plotTyp = varargin{2};
if length(varargin) == 5
    Work = varargin{3};
    argo = varargin{4};
    argoWork = varargin{5};
end

try unit = Work.unit;
   doxyLabel = sprintf('[O2] (%s)',unit);
catch
    try unit = argoWork.doxy_adjusted.units;
        doxyLabel = sprintf('[O2] (%s)',unit);        
    end    
end
presLabel = 'Pres (dbar)';

set(0,'currentFigure',hFig);

% =========================================================================
%% Plots
% =========================================================================
switch plotTyp

    case 1
        % -----------------------------------------------------------------
        % Plot the Argo Pressure VS the Argo Density
        % -----------------------------------------------------------------
        subplot(1,3,1);
        plot(Work.DENSITY-1000,argoWork.pres_adjusted.data,'.k')
        set(gca,'YDir','reverse');
        set(gca,'fontweight','bold');
        xlabel('argo density','fontweight','bold','fontsize',Work.fontsize);
        ylabel(presLabel,'fontweight','bold','fontsize',Work.fontsize);
        grid on;
        hold off;
        drawnow;
    
    case 2
        % -----------------------------------------------------------------
        % QC flag : plot the Argo doxy VS the choosen pressure, with QC flag colored
        % -----------------------------------------------------------------
        subplot(1,3,2);
        hold on;
        QC = [1 2 3 4 5 8 9];
        GRBColor = [[0 1 0];[0 1 1];[1 0.5 0];[1 0 0];[1 0 1]; [1 1 0];[1 1 1]];        
        okdim = strcmp(argoWork.pres_adjusted.dim,'N_PROF');        
        for i = 1 : size(argoWork.pres_adjusted.data, find(okdim))           
            for k = 1:length(QC)
                isok = strfind(argo.doxy_qc.data(i,:),num2str(QC(k)));
                plot(argo.doxy.data(i,isok),argoWork.pres_adjusted.data(i,isok),'.','color',GRBColor(k,:));
            end
        end
        set(gca,'YDir','reverse');
        title('[O2] data and their QC (1=G,2=C,3=O,4=R,5=M,8=Y,9=B)','fontweight','bold','fontsize',Work.fontsize);
        xlabel(doxyLabel,'fontweight','bold','fontsize',Work.fontsize); 
        ylabel(presLabel,'fontweight','bold','fontsize',Work.fontsize);
        set(gca,'fontweight','bold');
        grid on;
        hold off;
        drawnow;
        
    case 3
        % -----------------------------------------------------------------
        % Final Check : plot the Argo doxy VS the choosen pressure, and
        % color the remaining points.
        % -----------------------------------------------------------------
        subplot(1,3,3);
        hold on;
        
        %%%%P%%% represent only data with QC <= 3 in order to avoid
        %%%%representing DOXY concentration of margnitud ~ 10^4 mumol/kg
        QC = [1 2 3];   
        okdim = strcmp(argoWork.pres_adjusted.dim,'N_PROF');        
        for i = 1 : size(argoWork.pres_adjusted.data, find(okdim))                    
            for k = 1:length(QC)
                isok = strfind(argo.doxy_qc.data(i,:),num2str(QC(k)));
                plot(argoWork.doxy_adjusted.data(i,isok),argoWork.pres_adjusted.data(i,isok),'c.');      
            end
        end
        
        set(gca,'YDir','reverse');
        drawnow;
        title('Final [O2] data for correction in cyan','fontweight','bold','fontsize',Work.fontsize);
        xlabel(doxyLabel,'fontweight','bold','fontsize',Work.fontsize);        
        set(gca,'fontweight','bold');
        grid on;    
        hold off;       

        
end

% =========================================================================
%% Save the plot
% =========================================================================
if Work.savePlot
    save_name= hFig.Name(~isspace(hFig.Name));
    isSpecial = ~(save_name ~= ':' & save_name ~= '-');
    save_name(isSpecial)='_';
    saveFile = fullfile(Work.dirPlot,sprintf('%s_%s',save_name, Work.whichCorr));
    [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);
    pause(1) %avoid plotting google map in the last subplot of data presentation
end

end
