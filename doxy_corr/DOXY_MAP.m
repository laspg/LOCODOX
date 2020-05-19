% DOXY_MAP trace the float trajectory.
%
% SYNTAX
% [] = DOXY_MAP(argo,Work,plotTyp,REF,iok)
%
% DESCRIPTION
% DOXY_MAP trace the argo trajectory control plot for DOXY correction.
%
% INPUT
%   argo (structure)     Argo float structure (get directly from data).
%                        Example :
%                        argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%
%                         argo.pres=
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
%   Work (struct)        float working structure, issued and computed from
%                        the argo data and the reference data (WOA or REF).
%                        Carries plot informations among other things.
%                        Example:
%                         Work = 
%                                  readme: [1x1 struct]
%                                    unit: 'mumol/kg'
%                              R2threshold: 0.8000
%                                     wmo: 1901205
%                                makePlot: 1
%                                savePlot: 1
%
%
%   plotTyp (double)     indicate the type of plot you want to map
%                        1: plot all float data, coloring depending on time
%                        2: plot all float data, coloring depending on time
%                        + plot a selection of reference data (REF and iok)
%                        3: plot all float data, coloring depending on time
%                        + plot all the reference data (REF)
%
%   REF (struct)         Structure of the REFERENCE data (from
%                        gen_bdd_O2ref.m)
%                         REF = 
%                                          id: {361x1 cell}
%                                       presi: [1x401 double]
%                                         lon: [361x1 double]
%                                         lat: [361x1 double]
%                                         sta: [361x1 single]
%                                        juld: [361x1 single]
%                                        temp: [361x401 single]
%                                      theta0: [361x401 single]
%                                        psal: [361x401 single]
%                                        sig0: [361x401 single]
%                                        doxy: [361x401 single]
%                                        pres: [361x401 double]
%                                         sat: [361x401 single]
%                                        psat: [361x401 single]
%                                     doxy_CV: [361x1 single]
%                                     density: [361x1 single]
%
%   iok (boolean)        indices in REF of selected ref data you want to
%                        plot. Corresponding to REF.id(iok).
%
% OUTPUT
%
% CALL :
% plot_google_map.m, DOXY_PLOT_settingsToPrint
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 09/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:

function [] = DOXY_MAP(argo,Work,plotTyp,REF,iok)

% =========================================================================
%% Initialisation
% =========================================================================
hFig = figure('unit','normalized','OuterPosition',[0.35 0.50 0.35 0.50],...
    'Name',sprintf('ARGO GEO - %d',Work.wmo),'NumberTitle','off');

cycNum = argo.cycle_number.data;
Cmap = colormap(jet(length(cycNum)));

% Geographic extension
lonLim = [-180 180];
latLim = [-85 85];

% arguments
if nargin <= 2
    plotTyp = 1;
end
    
% =========================================================================
%% Plot the argo position on the world map
% =========================================================================
shg
subplot(2,1,1)

hold on;
%  Reads climatology coordinates and plot
for i=1:length(cycNum)
    plot(argo.longitude.data(i),argo.latitude.data(i),'*','color',Cmap(i,:,:),'MarkerSize',6);
end

% if plotTyp = 2, plot the reference data
if plotTyp == 2
    plot(REF.lon(iok),REF.lat(iok),'yo','Markersize',6);
elseif plotTyp == 3
    plot(REF.lon(1:length(REF.lon)),REF.lat(1:length(REF.lon)),'yo','Markersize',6);
    plot(REF.lon(iok),REF.lat(iok),'rs','Markersize',10);
end

% Resize
xlim(lonLim)
ylim(latLim)

% formatting
title(sprintf('Trajectory of the float %d', Work.wmo),...
    'fontweight','bold','FontSize',12)
xlabel('Longitude','fontweight','bold','fontsize',14);
ylabel('Latitude','fontweight','bold','fontsize',14);

plot_google_map('MapType','satellite','APIkey','AIzaSyCz69KOkqt4HoRKE0J6TvclvmWYXfoSH-s');
ylim(latLim)


% =========================================================================
%% Plot the argo position, ZOOM
% =========================================================================
subplot(2,1,2)
% Resize
xlim([min(argo.longitude.data)-4 max(argo.longitude.data)+4])
ylim([min(argo.latitude.data)-4 max(argo.latitude.data)+4])

hold on;
for i=1:length(cycNum)
    if i == 1
        plot(argo.longitude.data(i),argo.latitude.data(i),'r^','MarkerSize',10);
    elseif i == length(cycNum)
        plot(argo.longitude.data(i),argo.latitude.data(i),'rv','color',Cmap(i,:,:),'MarkerSize',10);
    else
        plot(argo.longitude.data(i),argo.latitude.data(i),'*','color',Cmap(i,:,:),'MarkerSize',6);
    end
    text(argo.longitude.data(i)+0.05,argo.latitude.data(i),num2str(i),'color','w','fontsize',8);
end

% if plotTyp = 2, plot the reference data
if plotTyp == 2
    plot(REF.lon(iok),REF.lat(iok),'yo','Markersize',6);
elseif plotTyp == 3
    plot(REF.lon(1:length(REF.lon)),REF.lat(1:length(REF.lon)),'yo','Markersize',6);
    plot(REF.lon(iok),REF.lat(iok),'rs','Markersize',10);
end
hold off

% formatting
title(sprintf('The same as above but zooming on the region where the float derived'),...
    'fontweight','bold','FontSize',12)
xlabel('Longitude','fontweight','bold','fontsize',14);
ylabel('Latitude','fontweight','bold','fontsize',14);
grid
plot_google_map('MapType','satellite','APIkey','AIzaSyCz69KOkqt4HoRKE0J6TvclvmWYXfoSH-s');

drawnow;

% =========================================================================
%% Save the figure
% =========================================================================
if Work.savePlot
    saveFile = fullfile(Work.dirPlot,sprintf('DOXY_MAP_%d_%s', Work.wmo, Work.whichCorr));
    [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);
end

end