% DOXY_NCEP_colocalize colocalizes the NCEP data over the argoTrajWork time/lat/lon
% regridded in a 6h daily grid.
%
% SYNTAX
% [NCEP] = DOXY_NCEP_colocalize(NCEP,argoTrajWork,CONFIG)
%
% DESCRIPTION
% DOXY_NCEP_colocalize colocalizes the NCEP data over the argoTrajWork time/lat/lon
% regridded in a 6h daily grid.
%
% INPUT
%   NCEP (structure)     filled with the variable and its attributes found
%                        in the NCEP climatology files
%                        Example: NCEP
%                              slp: [1x1 struct]
%                              air: [1x1 struct]
%                             rhum: [1x1 struct]
%
%                        with NCEP.slp = 
%                         dimorder: 'C'
%                              lat: [1x1 struct]
%                              lon: [1x1 struct]
%                              slp: [1x1 struct]
%                           recdim: 'time'
%                              obj: 'ObsInSitu'
%                        fillisnan: 1
%                             juld: [1x1 struct]
%                             time: 23884
%                     firstdimname: 'time'
%
%                       and NCEP.slp.slp =
%                            name: 'slp'
%                             dim: {'time'  'lat'  'lon'}
%                            data: [23884x73x144 single]
%                       long_name: '4xDaily Sea Level Pressure'
%                           units: 'mBar'
%                       precision: 0
%                                   ..........
%
%  argoTrajWork (structure)     float data structure (argoTrajWork) choose for the Doxy
%                               correction.
%                       Example: argoTrajWork
%                          float_serial_no: [1x1 struct]
%                            wmo_inst_type: [1x1 struct]
%                                     juld: [1x1 struct]
%                                 latitude: [1x1 struct]
%                                longitude: [1x1 struct]
%                                   ..........
%
% CONFIG (struct)       Configuration structure with data path,
%                       operator choices, ...
%                             CONFIG = 
%                        DataDir: '/home/oo26/coriolis/co05/co0508/dac/coriolis/'
%                        NCEPDataDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/data_input/'
%                        saveDataDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/results/data/'
%                        savePlotDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/results/plots/'
%                                   ..........
%                             
% OUTPUT
%   NCEP (structure)     the field coloc is filled with colocalized
%                        variables over the regridded argoTrajWork data (6h daily)
%                        Example: NCEP.coloc
%                                     juld: [1522x1 double]
%                                      slp: [1522x1 double]
%                                      air: [1522x1 double]
%                                     rhum: [1522x1 double]
%
% CALL :
%
% SEE ALSO
%   DOXY_corr_main, DOXY_NCEP_read

% HISTORY
%   $created: 09/05/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       22/03/2019      Marine GALLIAN, Altran Ouest
%                       Update function in order to manage new MC
%                       measurement code for in air data

function [NCEP] = DOXY_NCEP_colocalize(NCEP,argoTrajWork,CONFIG)

% Take inAir argoTrajWork time, lat and lon --MG
 for i = 1:length(argoTrajWork.ppox_doxy_adjusted.data)
     isInAir = ismember(argoTrajWork.measurement_code.data{i},CONFIG.inAirMC);
     if ~isempty(isInAir)
         argoTrajWork.juld.data{i} = argoTrajWork.juld.data{i}(isInAir);
         argoTrajWork.latitude.data{i} = argoTrajWork.latitude.data{i}(isInAir)  ;      
         argoTrajWork.longitude.data{i} = argoTrajWork.longitude.data{i}(isInAir);
     end
 end

% =========================================================================
% Interpolate argoTrajWork time, lat and lon over a 6h daily grid
% =========================================================================
if iscell(argoTrajWork.juld.data)
    juld = cellfun(@(x) x+datenum(strrep(argoTrajWork.juld.units,...
            'days since',''),'yyyy-mm-dd HH:MM:SS'),...
            argoTrajWork.juld.data,'UniformOutput',false);
    juld = cell2mat(juld);
else
    juld = argoTrajWork.juld.data + ...
        datenum(strrep(argoTrajWork.juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');
end


% Argo Traj Latitude and longitude
% The position of the PPOX is the first localization of Argo, for the same
% cycle_number
[myunik, idxunik] = unique(argoTrajWork.cycle_number.data);
idxunik = idxunik(isfinite(myunik));
idxunik(end+1) = length(argoTrajWork.cycle_number.data);

if iscell(argoTrajWork.latitude.data)
    lat = cell2mat(argoTrajWork.latitude.data);
    lon = cell2mat(argoTrajWork.longitude.data);
else
    lat = NaN(1,length(idxunik)-1);
    lon = NaN(1,length(idxunik)-1);
    for i = 1:length(idxunik)-1
        k = idxunik(i);
        isok = find(~isnan(argoTrajWork.latitude.data(k:idxunik(i+1))),1,'first');
        if ~isempty(isok)
            lat(i) = argoTrajWork.latitude.data(k+isok-1);
            lon(i) = argoTrajWork.longitude.data(k+isok-1);
        end
    end
end

% Define the 6h daily grid (4 times a day)
    % Force start/end date to 4-times daily grid (6h daily)
startdate = floor(min(juld)*4)/4; %--MG
enddate = ceil(max(juld)*4)/4; %--MG
    % new time grid : times at which to extract data
argo6h.juld = (startdate:1/4:enddate)';

% keep nan values in memory
%idxnan = isnan(juld); % Commented by T. Reynaud 06.11.2020
idxnan = isnan(juld.*lat);% Filter date without positions

% Interpolation over the argoTrajWork track, taking into account the new time grid
    % If any non unique values (example : repeated AET), add infinitesimal data
    % to array and then interpolate.
varToManage = {juld,lat,lon};
varNameToManage = {'juld','lat','lon'};
for v = 1:length(varToManage)
    y = varToManage{v};
    uniqueTab = unique(y);
    if length(uniqueTab) ~= length(y)
        eval([varNameToManage{v} '= y + linspace(0, 1, length(y))*1E-5;']);
    end
end

    % interpolate lat and lon over new grid time (extrapolated over NaN values)
argo6h.lat = interp1(juld(~idxnan),lat(~idxnan),argo6h.juld,'pchip',NaN);
argo6h.lat(isnan(argo6h.lat)) = interp1(juld(~idxnan),lat(~idxnan),argo6h.juld(isnan(argo6h.lat)),'nearest','extrap');
argo6h.lon = interp1(juld(~idxnan),lon(~idxnan),argo6h.juld,'pchip',NaN);
argo6h.lon(isnan(argo6h.lon)) = interp1(juld(~idxnan),lon(~idxnan),argo6h.juld(isnan(argo6h.lon)),'nearest','extrap');

% =========================================================================
% Find NCEP data near to argoTrajWork regridded data (time, lat and lon over a 6h
% daily grid)
% =========================================================================
varNames = fieldnames(NCEP.init);
varNames = varNames(~ismember(varNames,'nodata'));
ind_h=0;
for i = 1:length(varNames)
    nbdata = length(argo6h.juld);
    easyVar = varNames{i};
    dataNcepOnArgo.(easyVar).data = ones(nbdata,1)*NaN;
    dataNcepOnArgo.(easyVar).juld = argo6h.juld;
    dataNcepOnArgo.(easyVar).lat = argo6h.lat;
    dataNcepOnArgo.(easyVar).lon = argo6h.lon;
    
    % get NCEP data
    data = NCEP.init.(easyVar).(easyVar).data;
    lat = NCEP.init.(easyVar).lat.data;
    lon = NCEP.init.(easyVar).lon.data;
    juld = NCEP.init.(easyVar).juld.data + ...
        datenum(strrep(NCEP.init.(easyVar).juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');

    % look for the nearest NCEP juld from the argo6h juld
    isok = dsearchn(juld,argo6h.juld);

    %% FROM HENRY BITIIG: NCEPgetdata_netcdf4
    if length(isok) ~= length(unique(isok))
        % option 1: abort completely
        % Time index assigments to gridded data do not match uniquely!!
        % Check time resolution / multi-year assignment !
        % option 2: remove entries out of range and adjust time frame for
        % simulation
        if any(juld(isok) - argo6h.juld <= 1/4) && ind_h==0
            ind_h=1;
            idx=find(argo6h.juld>juld(end));% Added By TR 10.04.2020
            if ~isempty(idx)
                h = warndlg({sprintf('Float Time series longer than NCEP gridded data !!'),' '...
                strcat('INFO :  NCEP data are loaded up to: ',datestr(juld(end),31),'Float lasts: ',datestr(argo6h.juld(end),31)),'Update NCEP files!'});
            uiwait(h);
            else
            h = warndlg({sprintf('Time index assigments to gridded data do not match uniquely!!'),' '...
                strcat('INFO :  NCEP data are loaded for year ',num2str(CONFIG.ncepYears(1),'%4.4i'),'  to ',num2str(CONFIG.ncepYears(end),'%4.4i'),' ; see the configuration file. '),'Update NCEP files!'});
            uiwait(h);
            end
        elseif ind_h==0
            ind_h=1;
            h = warndlg({sprintf('Time index assigments to gridded data do not match uniquely!!'),' '...
                strcat('INFO :  NCEP data are loaded for year ',num2str(CONFIG.ncepYears(1),'%4.4i'),'  to ',num2str(CONFIG.ncepYears(end),'%4.4i'),' ; see the configuration file. '), 'Check time resolution / multi-year assignment or update NCEP files!'});
            uiwait(h);
        end
        % find indices to be removed at the end based on time difference
        % In particular, if the NCEP data doesn't cover the argo data time span in whole
        % or in part, operator is warned up
        rmind = find(abs(juld(isok) - argo6h.juld) >= 1/4/2); % >= 3 hours (half NCEP resolution)
        dataNcepOnArgo.(easyVar).data(rmind) = [];
        dataNcepOnArgo.(easyVar).juld(rmind)= [];
        dataNcepOnArgo.(easyVar).lat(rmind) = [];
        dataNcepOnArgo.(easyVar).lon(rmind) = [];
        nbdata = length(dataNcepOnArgo.(easyVar).juld);
    end
    
    %% Look for the nearest latitude and longitude indices
    ilat = NaN(nbdata,2);
    ilon = NaN(nbdata,2);
    for j = 1:nbdata
        isoklat = [find(lat >= dataNcepOnArgo.(easyVar).lat(j),1,'last') find(lat < dataNcepOnArgo.(easyVar).lat(j),1,'first')];
        if ~isempty(isoklat)
            % catch events when no proper index is obtained and NaN remains...
            ilat(j,:) = isoklat;
            ilon(j,:) = [find(lon - 360 < dataNcepOnArgo.(easyVar).lon(j),1,'last') ...
                         find(lon - 360 >= dataNcepOnArgo.(easyVar).lon(j),1,'first')];
        end
    end
    
    % Replace NaNs with 1 and keep index for later NaN replacement
    idxnan = isnan(ilat(:,1)) & isnan(ilat(:,2));
    ilat(idxnan,:) = 1;
    ilon(idxnan,:) = 1;
    
    % Simple 2D bilinear interpolation over lat/lon 
    % Henry Bittg comments :
    % the following lines attribute a linear weight in latitude and
    % longitude to the four corners surrounding the current interpolation
    % point. I made this to avoid a more time-consuming 2d bilinear
    % interpolation, which wasn???t readymade in Matlab available.
    wlat = 1 - abs(lat(ilat) - repmat(dataNcepOnArgo.(easyVar).lat,1,2))./...
        repmat(sum(abs(lat(ilat) - repmat(dataNcepOnArgo.(easyVar).lat,1,2)),2),1,2);
    wlon = 1 - abs(lon(ilon) - 360 - repmat(dataNcepOnArgo.(easyVar).lon,1,2))./...
        repmat(sum(abs(lon(ilon) - 360 - repmat(dataNcepOnArgo.(easyVar).lon,1,2)),2),1,2);
    
    % now extract scaled values
    for j = 1:nbdata
        dataNcepOnArgo.(easyVar).data(j) = sum(sum(double(squeeze(data(isok(j),...
            ilat(j,:),ilon(j,:)))).*(wlat(j,:)'*wlon(j,:))))/sum(sum(wlat(j,:)'*wlon(j,:)));
    end
    % replace NaN again (double-proof; wlat/wlon should be NaN already)
    dataNcepOnArgo.(easyVar).data(idxnan) = NaN;
    clear juld isok ilat ilon wlat wlon output
end

% =========================================================================
% Interpolate data to argoTrajWork profile date
% Finalize data
% =========================================================================
NCEP.coloc.juld = cell2mat(argoTrajWork.juld.data);
for i = 1:length(varNames)
    easyVar = varNames{i};
    NCEP.coloc.(easyVar) = interp1(dataNcepOnArgo.(easyVar).juld,dataNcepOnArgo.(easyVar).data,...
        NCEP.coloc.juld + datenum(strrep(argoTrajWork.juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS'));
end