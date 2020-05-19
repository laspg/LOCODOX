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
%   $created: 09/02/2020 $author: Thierry Reynaud
%   $Revision: version $Date: $author:
%   V2 24.04.2020 ==> Figure Title modified

function [] = DOXY_MAP2(CONFIG,argo,Work,plotTyp,REF,iok)
%function [] = DOXY_MAP2(CONFIG)
%load debug_map.mat;

% =========================================================================
%% Initialisation
% =========================================================================
addpath(CONFIG.M_MAPDir);
land_color = [238, 232, 168]/255;

hFig = figure('unit','normalized','OuterPosition',[0.35 0.50 0.35 0.50],...
    'Name',sprintf('ARGO GEO - %d',Work.wmo),'NumberTitle','off');

cycNum = argo.cycle_number.data;
Cmap = colormap(jet(length(cycNum)));

shg
% =========================================================================
%% Plot the argo position, ZOOM
% =========================================================================
% Resize
xlim=([floor(min(argo.longitude.data)-1) ceil(max(argo.longitude.data)+1)]);
ylim=([floor(min(argo.latitude.data)-1) ceil(max(argo.latitude.data)+1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading bottom topography:

if CONFIG.M_MAP_PLOT_BATHY
    [bathy2,lat_bathy,lon_bathy]=read_bathy2006([ylim(1) ylim(2) xlim(1) xlim(2)],CONFIG);
    
    lon_pl1=lon_bathy(1:end);
    lat_pl1=lat_bathy(1:end);
    i1=find(lon_pl1 >= xlim(1) & lon_pl1 <= xlim(2));
    j1=find(lat_pl1 >= ylim(1) & lat_pl1 <= ylim(2));
    dx2=abs(lon_pl1(2)-lon_pl1(1));
    dy2=abs(lat_pl1(2)-lat_pl1(1));
    
    lim=[];
    lim(:,1)=[-7000:500:-1000]';
    lim(:,2)=[-6000:500: 0000]';
    
    lim_cb=lim-0.001;
end

m_proj('miller','long',[xlim(1) xlim(2)], ...
    'lati',[ylim(1) ylim(2)]);
hold on;

if CONFIG.M_MAP_PLOT_BATHY
    
    [xi,yi]=meshgrid(lon_pl1(i1),lat_pl1(j1));
    [xpl_pcol,ypl_pcol]=m_ll2xy(xi-dx2/2,yi-dy2/2);
    [xpl1,ypl1]=m_ll2xy(xi,yi);
    
    bathy2_pl=nan*ones(size(bathy2(j1,i1)));
    for ic=1:1:size(lim,1)
        idx=find( bathy2(j1,i1)>= lim(ic,1) & bathy2(j1,i1) < lim(ic,2));
        bathy2_pl(idx)=ic;
    end
    
    [h1]=pcolor(xpl_pcol,ypl_pcol,bathy2_pl);
    caxis([1 size(lim,1)+1]);
    shading('interp');
    gray_map=flipud(gray(size(lim,1)));gray_map(find(gray_map(:,:) == 0))=0.05;
    c=colormap(gray_map);
    brighten(0.5);
end



m_gshhs_i('patch',land_color,'edgecolor','none');hold on;

lon_dd=argo.longitude.data;
lat_dd=argo.latitude.data;

iwrap = find(lon_dd>180);
if ~isempty(iwrap)
     lon_dd(iwrap) = lon_dd(iwrap) - 360;
end

iwrap = find(lon_dd<-180);
if ~isempty(iwrap)
     lon_dd(iwrap) = lon_dd(iwrap) + 360;
end

for i=1:length(cycNum)
    if i == 1
        m_plot(lon_dd(i),lat_dd(i),'r^','MarkerSize',10);
    elseif i == length(cycNum)
        m_plot(lon_dd(i),lat_dd(i),'rv','color',Cmap(i,:,:),'MarkerSize',10);
    else
        m_plot(lon_dd(i),lat_dd(i),'*','color',Cmap(i,:,:),'MarkerSize',6);
    end
    m_text(lon_dd(i)+0.05,lat_dd(i),num2str(i),'color','w','fontsize',8);
end

% if plotTyp = 2, plot the reference data
if plotTyp == 2
    lon_dd2=REF.lon(iok);
    lat_dd2=REF.lat(iok);
    
    iwrap = find(lon_dd2>180);
    if ~isempty(iwrap)
        lon_dd2(iwrap) = lon_dd2(iwrap) - 360;
    end
    
    iwrap = find(lon_dd2<-180);
    if ~isempty(iwrap)
        lon_dd2(iwrap) = lon_dd2(iwrap) + 360;
    end   
    
    m_plot(lon_dd2,lat_dd2,'yo','Markersize',6);
elseif plotTyp == 3

    lon_dd2=REF.lon(1:length(REF.lon));
    lat_dd2=REF.lat(1:length(REF.lon));
    
    iwrap = find(lon_dd2>180);
    if ~isempty(iwrap)
        lon_dd2(iwrap) = lon_dd2(iwrap) - 360;
    end
    
    iwrap = find(lon_dd2<-180);
    if ~isempty(iwrap)
        lon_dd2(iwrap) = lon_dd2(iwrap) + 360;
    end   
    
    lon_dd3=REF.lon(iok);
    lat_dd3=REF.lat(iok);
    
    iwrap = find(lon_dd3>180);
    if ~isempty(iwrap)
        lon_dd3(iwrap) = lon_dd3(iwrap) - 360;
    end
    
    iwrap = find(lon_dd3<-180);
    if ~isempty(iwrap)
        lon_dd3(iwrap) = lon_dd3(iwrap) + 360;
    end
       
    m_plot(lon_dd2,lat_dd2,'yo','Markersize',6);
    m_plot(lon_dd3,lat_dd3,'rs','Markersize',10);
end
hold off

% formatting
title(sprintf('Region where the float derived'),...
    'fontweight','bold','FontSize',12)
xlabel('Longitude','fontweight','bold','fontsize',CONFIG.fontsize);
ylabel('Latitude','fontweight','bold','fontsize',CONFIG.fontsize);
%m_grid('box','fancy','color',[0 0 0],'linestyle','-.');
if diff(xlim) <= 15
    dx=1;
elseif diff(xlim) > 15 & diff(xlim) <= 30
    dx= 2;
else
    dx=10;
end

if diff(ylim) <= 15
    dy=1;
elseif diff(ylim) > 15 & diff(ylim) <= 30
    dy= 2;
else
    dy=10;
end

m_grid('xtick',[xlim(1):dx:xlim(2)],'ytick',[ylim(1):dy:ylim(2)],'color',[0 0 0],'linestyle','-.')

drawnow;

% =========================================================================
%% Save the figure
% =========================================================================
if Work.savePlot
    saveFile = fullfile(Work.dirPlot,sprintf('DOXY_MAP_%d_%s', Work.wmo, Work.whichCorr));
    [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);
end

end