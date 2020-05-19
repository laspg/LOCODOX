% DOXY_PLOT_interpolation trace the data interpolated to control the job.
%
% SYNTAX
% DOXY_PLOT_interpolation(Work, argoWork, WOAorREF, icycle, refcyc)
%
% DESCRIPTION
% DOXY_PLOT_interpolation trace the control plot for DOXY interpolation.
% WOA is interpolated over argo data.
%
% INPUT
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the reference data (WOA, REF)
%                        Example:
%                        Work = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%   
%   argowork (structure)  float working structure issued from the argo data
%                         Example:
%                         argoWork = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                                   doxySat: [1x1 struct]
%                                  doxyPSat: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%
%
%   WOAorREF (double)      (1/2). The Reference data are WOA if WOAorREF=1,
%                          and REF if WOAorREF=2.
%
%   icycle (double)        OPTIONNAL. Give the index of the argo reference
%                          cycle for the REF correction type. The
%                          reference cyle is the corresponding argo float
%                          cycle referenced to the insitu Reference CTD.
%
%   refcyc (double)        OPTIONNAL. Give the reference cycle of the argo
%                          float for the REF correction type. the
%                          reference cyle is the corresponding argo float
%                          cycle referenced to the insitu Reference CTD.
% OUTPUT
%
% CALL :
%
% SEE ALSO
%   DOXY_corr_main, DOXY_interp_WOA_main, DOXY_interp_WOA_2_argo

% HISTORY
%   $created: 09/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%              v2   26/01/2017   Anne Piron, Altran Ouest
%                                argo - ameliorations 2017 phase 1
%                                (removed subplot)
%              v2.1 06/09/2017   Emilie Brion, Altran Ouest
%                                add argument refcyc
%              v2.2 09/07/2018   Emilie Brion, Marine Gallian, Altran Ouest
%                                Better management of grouped handles.
%                                Float data profiles colored by time
%              v2.3 06/03/2020   Thierry Reynaud
%                                Bug corrected for the output PNG file

function [] = DOXY_PLOT_interpolation(Work, argoWork, WOAorREF, icycle, refcyc)

% =====================================================================
%% Initialisation
% =====================================================================
if WOAorREF == 1
    cmpl = 'WOA';
elseif WOAorREF == 2
    cmpl = 'REF';
end

if nargin <= 4
    okdim = strcmp(argoWork.pres_adjusted.dim,'N_PROF');
    dim = size(argoWork.pres_adjusted.data);
    icycle = 1:dim(okdim);
end

% =====================================================================
%% PLOT
% =====================================================================
if strcmp(cmpl,'WOA')
    outerpos = [0.52 0.15 0.24 0.50];
elseif strcmp(cmpl,'REF')
    outerpos = [0.76 0.15 0.24 0.51];
end

% Create the plot
hFig = figure('unit','normalized','Outerposition',outerpos,...
    'Name', sprintf('%s DATA INTERPOLATED on Argo VSS - %d',cmpl,Work.wmo),'NumberTitle','off');
set(0,'CurrentFigure',hFig);

hold on;

%If it is not a correction made with REF
if WOAorREF ~= 2
    % a = plot(argoWork.doxy_adjusted.data, argoWork.pres_adjusted.data,'*b');
    % b = plot( Work.(['DOXY_' cmpl '_interpP']), argoWork.pres_adjusted.data(icycle,:),'*k');
    hgWoa = hggroup;
    plot( Work.(['DOXY_' cmpl '_interpP']), argoWork.pres_adjusted.data(icycle,:),'*k','Parent',hgWoa);
    nbr_prof = size(argoWork.doxy_adjusted.data,1);
    map = colormap(jet(nbr_prof));
    hgFloat = hggroup;
    for i_mar = 1:nbr_prof
        plot(argoWork.doxy_adjusted.data(i_mar,:), argoWork.pres_adjusted.data(i_mar,:),'Color',map(i_mar,:),'Parent',hgFloat);
    end
    
%If correction made with ref
elseif WOAorREF == 2
    nbr_prof = size(argoWork.doxy_adjusted.data,1);
    map = colormap(jet(nbr_prof));
    hgFloat = hggroup;
    for i_mar = 1:nbr_prof
        plot(argoWork.doxy_adjusted.data(i_mar,:), argoWork.pres_adjusted.data(i_mar,:),'Color',map(i_mar,:),'Parent',hgFloat);
    end
    
    hgWoa = hggroup;
    plot( Work.(['DOXY_' cmpl '_interpP']), argoWork.pres_adjusted.data(icycle,:),'k','LineWidth',3,'Parent',hgWoa);
    
    hgRef = hggroup;
    plot(argoWork.doxy_adjusted.data(icycle,:), argoWork.pres_adjusted.data(icycle,:),'r','LineWidth',3,'Parent',hgRef);    
end

grid on;

legend;

set(get(get(hgFloat,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
set(get(get(hgWoa,'Annotation'),'LegendInformation'),...hD
    'IconDisplayStyle','on');

if WOAorREF == 2
    set(get(get(hgRef,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');    
    legend('O2 float',['O2 ' cmpl],sprintf('O2 float, cycle %d',refcyc),'Location','Southwest');    
else
    legend(['O2 ' cmpl ' (in black)'],'O2 float (in colors)','Location','Southwest');
end

set(gca,'YDir','reverse');
set(gca,'fontweight','bold');
xlabel('[O2] (umol/kg)','fontweight','bold','fontsize',Work.fontsize);
ylabel('Pres (dbar)','fontweight','bold','fontsize',Work.fontsize);
title(sprintf('RAW DATA \n  %d', Work.wmo),...
    'fontweight','bold','fontsize',Work.fontsize);
ylim([0,nanmax(nanmax(argoWork.pres_adjusted.data))]);
hold off;
drawnow;

%% Save the plot
% =====================================================================
if Work.savePlot == 1
    if Work.presEff, presEffStr = 'okpreseff'; else, presEffStr = 'nopreseff'; end
    %saveFile =
    %fullfile(Work.dirPlot,sprintf('DOXY_PLOT_interpolation_%d_on%s_%s',Work.wmo,cmpl,presEffStr));
    saveFile = fullfile(Work.dirPlot,sprintf('DOXY_PLOT_interpolation_%d_on%s_%s',Work.wmo,Work.whichCorr,presEffStr));% TR 06/03/2020
    DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);
end