% DOXY_argo_write writes a B-monoprofile updated.
%
% SYNTAX
% DOXY_argo_write(DataDir,fileDir,saveDir,prefix,...
%                 argo1Struct, argo2Struct, argo3Struct, argo3Struct, Work)
%
% DESCRIPTION
% DOXY_argo_write writes a new monoprofile from old one, with the
% DOXY_adjusted field and some metadata fields updated.
%
% INPUT
%   DataDir (string)      full path of the data.
%
%   fileDir (string)      full path of the data.
%
%   saveDir (string)      full path of the sacing directory for data.
%
%   prefix (string)       prefix of the new file. Defined in the
%                         DOXY_config.txt
%
%   argo1Struct (struct)  Three data structures organised with the main profile
%   argo2Struct (struct)  for working is in argo1Struct (vertical
%   argo3Struct (struct)  sampling shceme of interest for the kind of
%   argo4Struct (struct)  correction), and the others are in the second ones.
%                         The program DOXY_argo_prepare_main.m describes
%                         the choice of the main and other profiles.
%                         Example :
%                             argo1Struct =
%
%                                     argo: [1x1 struct]
%                                 argoWork: [1x1 struct]
%                                     Work: [1x1 struct]
%                             with:
%                             argo1Struct.argo =
%                                          pres: [1x1 struct]
%                                 pres_adjusted: [1x1 struct]
%                                          doxy: [1x1 struct]
%                                       doxy_qc: [1x1 struct]
%                                             ...
%
%                             argo1Struct.argoWork =
%                                 pres_adjusted: [1x1 struct]
%                                 temp_adjusted: [1x1 struct]
%                                 psal_adjusted: [1x1 struct]
%                                       an_dens: [1x1 struct]
%                                       density: [1x1 struct]
%                                       doxy_qc: [1x1 struct]
%                                 doxy_adjusted: [1x1 struct]
%                                           sat: [1x1 struct]
%                                          psat: [1x1 struct]
%
%                              argo1Struct.Work =
%                                       readme: [1x1 struct]
%                                         unit: 'mumol/kg'
%                                   R2threshold: 0.8000
%                                          wmo: 1901205
%                                       sensor: 'Optode'
%                                     makePlot: 0
%                                     savePlot: 0
%                                    whichCorr: 'WOA'
%                                     DOXY_RAW: [166x117 single]
%                                        timar: [166x20 char]
%                                        datat: [166x1 double]
%                                      DENSITY: [166x117 single]
%                                        DEPTH: [166x117 double]
%                                      CAPTEUR: 'Optode'
%                                      DOXY_QC: [166x117 char]
%                                         DENS: [166x117 single]
%
%   Work (structure)      Float working structure, issued and computed from
%                         the argo data and the climatology data (WOA),
%                         completed with the annual drift.
%                         Example:
%                         Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                                     DRIFT: 4.8400
%   REF_ARGO (structure) List the WMO of argo float with DOXY sensor,
%                           and useful information : cycle, wmo, reference profile. 
%                           Example:
%                            REF_ARGO = 
%                                 cycle : [1x1 double]
%                                 refId : [1x1 string]
%                                 wmo : [1x1 double]
%
% OUTPUT
%
% CALL : 
%   DOXY_convert, NCR_file, NCW_file, DOXY_update_fields, check_profile_qc
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 19/01/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v2 25/05/2016   Emilie Brion, ALTRAN OUEST
%                       takes into account three Vertical Sampling Scheme
%                       (primary, secondary and nearsurface).
%       v3 26/01/2017   Emilie Brion, Anne Piron, Altran Ouest
%                       argo - ameliorations 2017 phase 1
%                       new repertories names to save final netcdf files, takes
%                       into account PSAT or DOXY correction in the
%                       directories names
%          21/03/2019   Marine GALLIAN, Altran Ouest, 
%                       new repertories names to save final netcdf files, takes
%                       into account drift computed on WOA or NCEP in the
%                       directories names and if a constant correction was
%                       applied


function [] = DOXY_argo_write(varargin)
% DataDir,fileDir,saveDir,prefix,argo1Struct, argo2Struct, argo3Struct, argo4Struct, Work,REF_ARGO)

DataDir = varargin{1};
fileDir = varargin{2};
saveDir = varargin{3};
prefix = varargin{4};
argo1Struct = varargin{5};
argo2Struct = varargin{6};
argo3Struct = varargin{7};
argo4Struct = varargin{8};
Work = varargin{9};

if length(varargin) == 10
    REF_ARGO= varargin{10};
else
    REF_ARGO=[];
end


% Get the list of Bio files
bio_files = dir(fullfile(DataDir,fileDir,'B*nc'));
bio_files = {bio_files(:).name}';
new_bio_files = strrep(bio_files,'BR',prefix);

% Loop over the files
for k = 1:length(bio_files)
    % Read Bio-file
    % ---------------------------------------------------------------------
    [monoProf,Dim,Globatt] = NCR_file(fullfile(DataDir,fileDir,bio_files{k}));
    cycle_monoProf = unique(monoProf.cycle_number.data);
    dir_monoProf = unique(monoProf.direction.data);
        
    % Updates scientific calib fields, history fields,
    % data_state_indicator, data_mode, date_update.
    
    [monoProf, Dim] = DOXY_update_fields(Work,monoProf,Dim,REF_ARGO,k);
    monoProf = fillValue_2_nan(monoProf);
    monoProf.hereDoxyFull.name = 'hereDoxy';
    monoProf.hereDoxyFull.dim = {'N_PROF'};
    monoProf.hereDoxyFull.data = false(1,monoProf.n_prof);
    monoProf.hereDoxyFull.long_name = 'Index of the presence of the variable DOXY';
    monoProf.hereDoxyFull.standard_name = 'DOXY presence index';
    isok = false(1,monoProf.n_prof);
    for i = 1:monoProf.n_prof
        test = strcmp(cellstr(deblank(squeeze(monoProf.station_parameters.data(i,:,:)))),'DOXY');
        if any(test)
            isok(i) = true;
        end
    end
    % keep all the doxy occurence (different from the
    % argoStruct.argo.hereDoxy.data, where only the first occurence is kept)
    monoProf.hereDoxyFull.data = isok;
    
    if sum(isok) > 2
        % More than two differents profiles carrying DOXY
        % example : two secondary profiles with DOXY, and one nearSurface
        % profile.
        VSS_Name = regexp(cellstr(monoProf.vertical_sampling_scheme.data(isok,:)),'.*(?= sampling)','match');
        VSS_Name = char([VSS_Name{:}]);
        [uniqueVSS ia ib] = unique(VSS_Name,'rows','stable');
        isok_tmp = find(isok);
        if size(uniqueVSS,1) ~= size(uniqueVSS,1)
            fprintf('VSS carrying DOXY is not unique for one kind of VSS!\n');
            % keep the first encoutered
            for n = 1:size(VSS_Name,1)
                if sum(ismember(ib,ia(n))) > 1
                    % get the first occurence as reference
                    idx = find(ib == ia(n));
                    idx = idx(2:end);
                    isok(isok_tmp(idx)) = false;
                end
            end
        end
    end
    monoProf.hereDoxy.name = 'hereDoxy';
    monoProf.hereDoxy.dim = {'N_PROF'};
    monoProf.hereDoxy.data = false(1,monoProf.n_prof);
    monoProf.hereDoxy.long_name = 'Index of the presence of the variable DOXY - First occurence';
    monoProf.hereDoxy.standard_name = 'DOXY presence first occurence index';
    monoProf.hereDoxy.data = isok;
    
    % Fill the doxy_adjusted data updated
    % WARNING: float 6901627 2DO : not always the doxy and doxy2 parameter, sometimes only doxy2
    % ---------------------------------------------------------------------
    possibleVSS = {'Primary','Secondary','Near-surface','Other'};
    if any(monoProf.hereDoxy.data)
        raf = false(size(monoProf.hereDoxy.data));
        for v = 1:length(possibleVSS)
            % look for the argoNstruct corresponding to the monoprof,
            % respect to the VSS
            argoStruct = eval(sprintf('argo%dStruct',v));
            localVSS = argoStruct.VSS;
            isokCycNum = argoStruct.argo.cycle_number.data == cycle_monoProf & strcmp(cellstr(argoStruct.argo.direction.data),dir_monoProf);
            if v ~= 4
                isokVSS = ~cellfun('isempty',strfind(cellstr(monoProf.vertical_sampling_scheme.data),localVSS)) & (monoProf.hereDoxy.data)';
                raf = raf | isokVSS';
            else
                isokVSS = ~raf & monoProf.hereDoxyFull.data;
            end
            if any(isokVSS) && any(isokCycNum & argoStruct.argo.hereDoxy.data)
                monoProf.doxy_adjusted.data(isokVSS,:) = argoStruct.argoWork.doxy_adjusted.data(isokCycNum & argoStruct.argo.hereDoxy.data,1:Dim.n_levels.dimlength);
                monoProf.doxy_adjusted_qc.data(isokVSS,:) = argoStruct.argoWork.doxy_qc.data(isokCycNum & argoStruct.argo.hereDoxy.data,1:Dim.n_levels.dimlength);
                monoProf.doxy_adjusted_error.data(isokVSS,:) = argoStruct.argoWork.doxy_adjusted_error.data(isokCycNum & argoStruct.argo.hereDoxy.data,1:Dim.n_levels.dimlength);
            end
        end
    end
    monoProf = nan_2_fillValue(monoProf);
    
    % Update the PROFILE_QC from the DOXY_ADJUSTED_QC
    % ---------------------------------------------------------------------
    monoProf = check_profile_qc(monoProf,'doxy');
    
    % Fields hereDoxy and hereDoxyFull no more useful
    % ---------------------------------------------------------------------
    monoProf = rmfield(monoProf,'hereDoxy');
    monoProf = rmfield(monoProf,'hereDoxyFull');

    % Write the new NetCDF file, with corrected DOXY fields and
    % metadata udpated
    % ---------------------------------------------------------------------
    % Name of new files : BR become BD
    dirout = fullfile(saveDir,'NC',Work.dirout,num2str(Work.wmo));
    if ~exist(dirout,'dir')
        mkdir(dirout);
    end
    ficout = fullfile(dirout,new_bio_files{k});
    NCW_file(monoProf,Dim,ficout,Globatt);
end