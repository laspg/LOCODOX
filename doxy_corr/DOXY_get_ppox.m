% DOXY_get_ppox gets the PPOX and useful information associated to.
%
% SYNTAX
% [argoTrajWork] = DOXY_get_ppox(argoTrajWork,Dim,CONFIG)
%
% DESCRIPTION
% DOXY_get_ppox gets the PPOX and useful information associated to.
%       - Get Ascent End Time
%       - Get the argos localisation
%       - Check first available argos localisation
%       - Get the corresponding PPOX_DOXY
%
% NOTE
%  * Old TRAJ version: PPOX has to be read in the profile DOXY of BTraj.nc,
% as the last profile point, dated with AET.
% Profile measurement_code = 590
%  * Intermediate TRAJ version: directly read PPOX_DOXY, no matter the
% measurement_code
%  * Future TRAJ Version: directly read PPOX_DOXY, no matter the measurement_code.
% The PPOX_DOXY could be directly dated and localize, depending on the
% float version
%
% INPUT
%     argoTrajWork (structure)    float trajectory intermediate working
%                                 data
%
%     Dim (structure)             float trajectory intermediate working
%                                 dimensions
%
%     CONFIG (struct)       Configuration structure with data path,
%                           operator choices, ...
%                             CONFIG = 
%                        DataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/data_coriolis/INAIR/'
%                    NCEPDataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04\Data_NCEP'
%                    saveDataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/data/'
%                    savePlotDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/plots/'
%                              ...
%                              ...
%                    inAirFormat: 0
%                        inAirMC: 1090
%                              ...
%
% OUTPUT
%     argoTrajWork (structure)    float trajectory intermediate working
%                                 data completed
%
% CALL : 
%   NCR_file
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 10/05/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:

% function [argoTrajWork] = DOXY_get_ppox(argoTrajWork,Dim,CONFIG,argo1Struct,argo2Struct,argo3Struct, argo4Struct)
function [argoTrajWork] = DOXY_get_ppox(argoTrajWork,Dim,CONFIG)

% =========================================================================
%% Find the index of the cycle
% =========================================================================
[cycNum, idxunik] = unique(argoTrajWork.cycle_number.data);
idxunik = idxunik(isfinite(cycNum));
idxunik(end+1) = length(argoTrajWork.cycle_number.data)+1;

% =========================================================================
%% Get Ascent End Time (AET, Measurement_code = 600)
% Ascent End Time is the time associated to the PPOX measurement.
% for cycle#0, AET doesn't exist. Get the TSD or the FMT or the first
% localisation argos, depends on avalaibility
% =========================================================================
aet = NaN(1,length(idxunik)-1);
MCutil = [600 700 701 702 703];
imc = 1;
for i=1:length(idxunik)-1
    k = idxunik(i);
    while imc <= length(MCutil)
        isok = find(argoTrajWork.measurement_code.data(k:idxunik(i+1)-1) == MCutil(imc),1,'first');
        if ~isempty(isok) && ~isnan(argoTrajWork.juld.data(k+isok-1))
            aet(i) = argoTrajWork.juld.data(k+isok-1);
            break
        else
            imc = imc+1;
        end
    end
end

% =========================================================================
%% Get the argos localisation (position and time)
% The corresponding measurement code goes from 700 to 704. keep only the
% one with both lat/lon and juld
% =========================================================================
argosLoc.lat = NaN(length(idxunik) -1,100);
argosLoc.lon = NaN(length(idxunik) -1,100);
argosLoc.juld = NaN(length(idxunik) -1,100);
MCutil = [700 702 703 704];
for i = 1:length(idxunik) - 1
    k = idxunik(i);
    imc = 1;
    imc_pos = false;
    imc_juld = false;
    while imc <= length(MCutil)
        isok = find(argoTrajWork.measurement_code.data(k:idxunik(i+1)-1) == MCutil(imc));
        if ~isempty(isok)
            if any(~isnan(argoTrajWork.latitude.data(k+isok-1))) & ~imc_pos
                imc_pos = true;
            end
            if any(~isnan(argoTrajWork.juld.data(k+isok-1))) & ~imc_juld
                imc_juld = true;
            end
            if imc_pos
                argosLoc.lat(i,1:length(isok)) = argoTrajWork.latitude.data(k+isok-1);
                argosLoc.lon(i,1:length(isok)) = argoTrajWork.longitude.data(k+isok-1);
                imc_pos = false;
            end
            if imc_juld
                argosLoc.juld(i,1:length(isok)) = argoTrajWork.juld.data(k+isok-1);
                imc_juld = false;
            end
        end
        imc = imc+1;
    end
end

% -----------------------------------------------------------------
% Check first available argos localisation (non NaN) (FLA)
% Note: The correct position_qc (qc <= 3) has been selected when
% argoTrajWork have been created.
% -----------------------------------------------------------------
fla.juld = NaN(1,length(idxunik)-1);
fla.lat = NaN(1,length(idxunik)-1);
fla.lon = NaN(1,length(idxunik)-1);
for i = 1:length(idxunik)-1
    isok = find(isfinite(argosLoc.juld(i,:)),1,'first');
    if ~isempty(isok)
        fla.juld(i) = argosLoc.juld(i,isok);
    end
    isok = find(isfinite(argosLoc.lat(i,:)),1,'first');
    if ~isempty(isok)
        fla.lat(i) = argosLoc.lat(i,isok);
        fla.lon(i) = argosLoc.lon(i,isok);
    end
end
% =========================================================================
%% Spread the data by cycle
% =========================================================================
% PPOX is carried by MC=1090.
% Depending on floats, the PPOX in the traj could be :
%   * the last unpumped measurement reported in the trajectory.
%   * the true in-air measurements
% As the inAir measuremnt is not 


% -----------------------------------------------------------------
% Get the juld, lat and lon for surface measurement
% -----------------------------------------------------------------
argoInAirWork.juld = argoTrajWork.juld;
argoInAirWork.juld.dim = {Dim.n_cycle.name};
argoInAirWork.juld.data = cell(size(cycNum));
argoInAirWork.latitude = argoTrajWork.latitude;
argoInAirWork.latitude.data = cell(size(cycNum));
argoInAirWork.latitude.dim = {Dim.n_cycle.name};
argoInAirWork.longitude = argoTrajWork.longitude;
argoInAirWork.longitude.dim = {Dim.n_cycle.name};
argoInAirWork.longitude.data = cell(size(cycNum));
argoInAirWork.cycle_number = argoTrajWork.cycle_number;
argoInAirWork.cycle_number.dim = {Dim.n_cycle.name};
argoInAirWork.cycle_number.data = cycNum;

argoInAirWork.ppox_doxy_adjusted = argoTrajWork.ppox_doxy_adjusted;
argoInAirWork.ppox_doxy_adjusted.dim = {Dim.n_cycle.name};
argoInAirWork.ppox_doxy_adjusted.data = cell(size(aet));
argoInAirWork.psat.dim = {Dim.n_cycle.name};
argoInAirWork.psat.data = cell(size(cycNum));
argoInAirWork.psal_adjusted = argoTrajWork.psal_adjusted;
argoInAirWork.psal_adjusted.dim = {Dim.n_cycle.name};
argoInAirWork.psal_adjusted.data = cell(size(cycNum));
argoInAirWork.temp_adjusted = argoTrajWork.temp_adjusted;
argoInAirWork.temp_adjusted.dim = {Dim.n_cycle.name};
argoInAirWork.temp_adjusted.data = cell(size(cycNum));
argoInAirWork.pres_adjusted = argoTrajWork.pres_adjusted;
argoInAirWork.pres_adjusted.dim = {Dim.n_cycle.name};
argoInAirWork.pres_adjusted.data = cell(size(cycNum));
argoInAirWork.juld_ppox_doxy_adjusted = argoTrajWork.juld;
argoInAirWork.juld_ppox_doxy_adjusted.dim = {Dim.n_cycle.name};
argoInAirWork.juld_ppox_doxy_adjusted.data = cell(size(cycNum));
argoInAirWork.measurement_code = argoTrajWork.measurement_code;
argoInAirWork.measurement_code.dim = {Dim.n_cycle.name};
argoInAirWork.measurement_code.data = cell(size(cycNum));

argoInAirWork.trajVar = argoTrajWork.trajVar;
% =========================================================================
%% Get the corresponding PPOX
%  * Old TRAJ version: PPOX has to be read in the profile DOXY of BTraj.nc,
% as the last profile point, dated with AET.
% Profile measurement_code = 590
%  * Intermediate version: directly read PPOX_DOXY, no matter the
% measurement_code
%  * Future Version: directly read PPOX_DOXY, no matter the measurement_code.
% The PPOX_DOXY could be directly dated and localize, depending on the
% float version
% =========================================================================
for i = 1:length(idxunik)-1
    k = idxunik(i);
    cyc = k:idxunik(i+1)-1;
    %isok = argoTrajWork.measurement_code.data(cyc) == CONFIG.inAirMC;
    isok = ismember(argoTrajWork.measurement_code.data(cyc),[CONFIG.inWaterMC CONFIG.inAirMC]);
    if any(isok)
        argoInAirWork.ppox_doxy_adjusted.data{i} = argoTrajWork.ppox_doxy_adjusted.data(cyc(isok));
        argoInAirWork.psal_adjusted.data{i} = argoTrajWork.psal_adjusted.data(cyc(isok));
        argoInAirWork.temp_adjusted.data{i} = argoTrajWork.temp_adjusted.data(cyc(isok));
        argoInAirWork.pres_adjusted.data{i} = argoTrajWork.pres_adjusted.data(cyc(isok));
        argoInAirWork.juld_ppox_doxy_adjusted.data{i} = argoTrajWork.juld.data(cyc(isok));
        % Keep juld for InAir measurement if the number of juld is the same that
        % the number of inAir measurements (juld carried by MC=1090, i.e.
        % dated inair measurement). Otherwise, repeat AET.
        argoInAirWork.measurement_code.data{i} = argoTrajWork.measurement_code.data(cyc(isok));
        if any(isnan(argoTrajWork.juld.data(cyc(isok))))
            % set AET for InAir measurement date
            argoInAirWork.juld.data{i} = repmat(aet(i),1,sum(isok));
        else
            argoInAirWork.juld.data{i} = argoTrajWork.juld.data(cyc(isok));            
        end
        % Position : as it is now, the first available position to locate
        % the InAir measurement is the first Argo localisation. The program
        % is made to take into account a possible evolution of Argo, if the
        % InAir measurement could be locate by itself with specific argos
        % localisation.
        if any(isnan(argoTrajWork.latitude.data(cyc(isok))))
            argoInAirWork.latitude.data{i} = repmat(fla.lat(i),1,sum(isok));
            argoInAirWork.longitude.data{i} = repmat(fla.lon(i),1,sum(isok));            
        end             
    end
end

argoTrajWork = argoInAirWork;
% the version of floats  for wwhich the last unpumped measurement are used
% as "in-air" measurement (MC = 590)
% for i = 1:length(idxunik)-1
%     k = idxunik(i);
%     isok = find(argoTrajWork.measurement_code.data(k:idxunik(i+1)-1) == 590,1,'last');
%     if ~isempty(isok)
%         argoInAirWork.doxy_adjusted.data(i) = argoTrajWork.doxy_adjusted.data(k+isok-1);
%         argoInAirWork.psal_adjusted.data(i) = argoTrajWork.psal_adjusted.data(k+isok-1);
%         argoInAirWork.temp_adjusted.data(i) = argoTrajWork.temp_adjusted.data(k+isok-1);
%         argoInAirWork.pres_adjusted.data(i) = argoTrajWork.pres_adjusted.data(k+isok-1);
%     end
% end

% the version of floats for which the in-air measurement are coded with a
% "surface" measurement code (MC reserved for specific timing events in
% air)
% for i = 1:length(idxunik)-1
%     k = idxunik(i);
%     cyc = k:idxunik(i+1)-1;
%     isok = argoTrajWork.measurement_code.data(cyc) == CONFIG.inAirMC;
%     if any(isok)
%         argoInAirWork.ppox_doxy_adjusted.data(i,1:sum(isok)) = argoTrajWork.ppox_doxy_adjusted.data(cyc(isok));
%         argoInAirWork.psal_adjusted.data(i,1:sum(isok)) = argoTrajWork.psal_adjusted.data(cyc(isok));
%         argoInAirWork.temp_adjusted.data(i,1:sum(isok)) = argoTrajWork.temp_adjusted.data(cyc(isok));
%         argoInAirWork.pres_adjusted.data(i,1:sum(isok)) = argoTrajWork.pres_adjusted.data(cyc(isok));
%     end
% end


% % =========================================================================
% % Be aware of the time and spatial lag between AET and FLA.
% % =========================================================================
% % Estimate the argo surface drift speed
% % time and distance between all the argos Localisation
% argosLoc.diff.loc = NaN(length(idxunik)-1,size(argosLoc.lat,2)-1);
% argosLoc.diff.juld = NaN(length(idxunik)-1,size(argosLoc.lat,2)-1);
% surfSpeed = NaN(length(idxunik)-1,size(argosLoc.lat,2)-1);
% for i = 1:length(idxunik)-1
%     if sum(isfinite(argosLoc.lat(i,:))) > 1
%         argosLoc.diff.loc(i,:) = dist(argosLoc.lat(i,:),argosLoc.lon(i,:));
%         argosLoc.diff.juld(i,:) = diff(arargoInAirWorkargoInAirWorkgosLoc.juld(i,:));
%         surfSpeed(i,:) = argosLoc.diff.loc(i,:)./(argosLoc.diff.juld(i,:)*24*60*60);
%     else
%         fprintf('No more than one surface position for the cycle %d => no drift speed calculation\n', ...
%           argoTrajWork.cycle_number.data(i));  
%     end
% end
% surfSpeed(isinf(surfSpeed)) = NaN;
% [a,b] = find(isfinite(surfSpeed),1,'Last');
% surfSpeed = surfSpeed(:,1:b);
% 
% 
% % Eliminate data too much different
% % for i = 1:length(idxunik)-1
% %     speed = toto(i,:);   
% %     tmpStd = mynanstd(speed); avg = nanmean(speed); speed(abs(speed) > avg + 2.8*tmpStd) = NaN;    
% %     tmpStd = mynanstd(speed); avg = nanmean(speed); speed(abs(speed) > avg + 2.8*tmpStd) = NaN;    
% %     tmpStd = mynanstd(speed); avg = nanmean(speed); speed(abs(speed) > avg + 2.8*tmpStd) = NaN;    
% %     surfSpeed(i,:) = speed;
% % end
% 
% % Check for outliers data : data over 3 m/s are eliminated (cf ANDRO
% % generation)
% surfSpeed(surfSpeed > 3) = NaN;
% 
% % Time interval between AET and first available argos localisation
% deltaT = fla.juld - aet;
% 
% % mean surface speed
% meanSurfSpeed = nanmean(surfSpeed,2);
% 
% 
