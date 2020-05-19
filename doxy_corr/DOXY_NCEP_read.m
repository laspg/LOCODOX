% DOXY_NCEP_read reads the NCEP climatology file and convert variables in
% useful units.
%
% SYNTAX
% [NCEP] = DOXY_NCEP_read(CONFIG, years)
%
% DESCRIPTION
% DOXY_NCEP_read reads the NCEP climatology files and stocks it in a structure.
% The parameters are converted in a useful unit. Updates the NCEP local
% data with ftp data if CONFIG.ncepDoUpdate = 1.
% Only NCEP data for a selected years are read.
%
% INPUT
%   CONFIG (structure)   Structure read from configuration file
%                        Example:CONFIG = 
%                     CONFIGURATION_FILE: '/home1/homedir4/perso/ebrion/NAOS/NAOS_2016(1)/programmes/MATLAB/ebrion/DOXY_config.txt'
%                                ncepFtp: 'ftp.cdc.noaa.gov'
%                             ncepFtpDir: 'Datasets/ncep.reanalysis/'
%                          ncepFtpSubDir: {'surface'  'surface'  'surface'}
%                              ncepFiles: {'slp'  'air.sig995'  'rhum.sig995'}
%                           ncepDoUpdate: '0'
%                                       ...
%
%   years (double)       (optionnal). Select year over the NCEP data, to
%                        read only this(these) year(s). If non given, all
%                        the available NCEP years are read.
%                        Example: years = [2002:2016].
%
% OUTPUT
%   NCEP (structure)     filled with the variable and its attributes found
%                        in the NCEP climatology files
%                        Example: NCEP = 
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
%         least_significant_digit: -1
%                         GRIB_id: 2
%                       GRIB_name: 'PRMSL'
%                        var_desc: 'Sea Level Pressure'
%                      level_desc: 'Sea Level'
%                       statistic: 'Individual Obs'
%                     parent_stat: 'Other'
%                      FillValue_: -9.9692e+36
%                    actual_range: [93190 111720]
%                     valid_range: [87000 115000]
%                         dataset: 'NCEP Reanalysis'
%                            type: 5    
%                                   ..........
%
% CALL :
%   NCR_file
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 09/05/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
% 10.04/2020 T. Reynaud Download only CONFIG.ncepGetYears



function [NCEP] = DOXY_NCEP_read(CONFIG, years)

% CONFIG.NCEPdoUpdate [1/0] : do the NCEP udpate from the NCEP ftp

if CONFIG.ncepDoUpdate
% =========================================================================
%% UPDATE : Load S data from web
% ATTENTION : seulement un update, ou carr??ment tout charger ?
% =========================================================================
    % open ftp connection
    cd(CONFIG.NCEPDataDir);% Added by Thierry Reynaud 09.04.2020
    ncepftp=ftp(CONFIG.ncepFtp);
    cd(ncepftp,CONFIG.ncepFtpDir);
    for i=1:length(CONFIG.ncepFiles)
        disp(['load current ' CONFIG.ncepFiles{i}]);
        if i == 1
            cd(ncepftp,CONFIG.ncepFtpSubDir{i});
        else
            cd(ncepftp,['../',CONFIG.ncepFtpSubDir{i}]);
        end
        
        %VT 14/01/2020 ==> Modified by T. Reynaud 10.04.2020
        for iy=1:length(CONFIG.ncepGetYears)
            
%         [CONFIG.ncepFiles{i} '.' datestr(now,'yyyy') '.nc']
%         mget(ncepftp,[CONFIG.ncepFiles{i} '.' datestr(now,'yyyy') '.nc'],CONFIG.NCEPDataDir);
          display([CONFIG.ncepFiles{i} '.' num2str(CONFIG.ncepGetYears(iy)) '.nc'])
          mget(ncepftp,[CONFIG.ncepFiles{i} '.' num2str(CONFIG.ncepGetYears(iy)) '.nc'],CONFIG.NCEPDataDir);
        end
    end
    close(ncepftp);
    disp('NCEP files updated')

end
% =========================================================================
%% Read S data
% =========================================================================
fprintf('INFO >>>>>    Reading NCEP data :  \n    ')

for v = 1:length(CONFIG.ncepFiles)
    fileList = dir(fullfile(CONFIG.NCEPDataDir,sprintf('%s*.nc',CONFIG.ncepFiles{v})));
    % Keep years
    if nargin == 2
        tmp = regexp({fileList.name}','[0-9]{4}','match');
        tmp = [tmp{:}]';
        tmp = str2num(char(tmp));
        okYear = ismember(tmp,years);
        fileList = fileList(okYear);
    end
    
    easyVar = regexp(CONFIG.ncepFiles{v},'\.','split');
    easyVar = easyVar{1};
    if isempty(fileList)
        fprintf('WARNING : there is no %s NCEP data for the defined years CONFIG.ncepYears %s\n',easyVar,sprintf('%d,',years));
        NCEP.nodata = 1;
        NCEP.(easyVar) = [];        
    else
        NCEP.nodata = 0;
        Smult = [];
        Smult = [];
        DIMmult = [];
        fprintf('\t Read %s file ... ...\n',easyVar);
        for i = 1:length(fileList)
            [S, DIMN] = NCR_file(fullfile(CONFIG.NCEPDataDir,fileList(i).name));
            % put time to juld
            S.juld = S.time;
            S.juld.data = S.time.data/24;% + datenum([1800 01 01 0 0 0]);
            S.juld.name = 'juld';
            S.juld.name = 'juld';
            S.juld.units = 'days since 1800-01-01 00:00:0.0';
            S.juld.long_name = 'julian day';
            S.juld.standard_name = 'julian day';
            S.juld.actual_range = S.juld.actual_range / 24;
            S = rmfield(S,'time');
            [Smult,DIMmult] = cat_profile_dim(Smult,S,DIMmult,DIMN,'time');
        end
        clear S;
        NCEP.(easyVar) = Smult;
        NCEP.(easyVar) = fillValue_2_nan(NCEP.(easyVar));
    end
end

% =========================================================================
%% Convert to useful units
% =========================================================================
varNames = fieldnames(NCEP);
for i = 1:length(varNames)
    if strcmp(varNames{i},'slp')
        % Pa -> hPa, mbar
        NCEP.(varNames{i}).(varNames{i}).data = NCEP.(varNames{i}).(varNames{i}).data/100;
        NCEP.(varNames{i}).(varNames{i}).units = 'mBar';
    elseif strcmp(varNames{i},'air')
        % K -> celsius degree
        NCEP.(varNames{i}).(varNames{i}).data = NCEP.(varNames{i}).(varNames{i}).data - 273.15;
        NCEP.(varNames{i}).(varNames{i}).units = 'degC';        
    end
end
cd(CONFIG.LocodoxMainDir);% Added by Thierry Reynaud 10.04.2020