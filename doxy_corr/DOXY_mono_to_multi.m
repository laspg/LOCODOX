% DOXY_mono_to_multi create multiprof structure from monoprof bio_files and
% core_files
%
% SYNTAX
% [Smult, DIMmult] = DOXY_mono_to_multi(DataDir,wmo,replist,VSS,LogDir)
%
% DESCRIPTION
% DOXY_mono_to_multi concatenate the mono-profile argo float data for Doxy
% correction, as a multi-profil-like profiles, from Bio files and Core
% files.
% Different steps :
% * Get and read the Bio files ("B" files)
% * Get and read the Core files ("D" or "R" files)
% * Don't take into account descendant profiles.
% * Remove negative and constant pressure
% * Compute Depth and Density
% * Apply quality flags to DOXY
% * Convert Doxy data to new unit
%
%
% INPUT
%     DataDir (string)      full path of the data.
%     
%     wmo (double)          WMO of the float
%
%     replist (string)      wmo/profile/
%
%     VSS (string)          Vertical Sampling Scheme choice
%
%     LogDir (string)       fullpath of log directory
%
% OUTPUT
%   Smult (structure)   Contains all information found in the NetCDF
%                       monoprofile files, concatenated in a "multiprofile"
%                       structure.
%                       Example :
%                           Smult.temp.name = 'TEMPERATURE' name of the variable
%                           Smult.temp.dim  = {'N_PROF','N_LEVEL'}
%                           Smult.temp.data =  n_prof x n_level value of temperature
%
%   DIMmult (structure) Contains all the dimensions information found in
%                       the NetCDF monoprofile files, concatenated in a
%                       "multiprofile" structure.
%                       Example :
%                           DIMmult.n_prof.name = 'N_PROF'     name of the dimension
%                           DIMmult.n_prof.dimlength = 88      length of the dimension
%
% CALL : 
%   NCR_file, findstr_tab, extract_profile_dim, check_FirstDimArray,
%   cat_profile_dim
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 11/01/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%              v3  26/01/2017  Anne Piron, Altran Ouest
%                              argo - ameliorations 2017 phase 1
%
function [Smult, DIMmult] = DOXY_mono_to_multi(DataDir,wmo,replist,VSS,LogDir)

% =========================================================================
%% List and prepare files
% =========================================================================
% Bio files
bio_files = dir(fullfile(DataDir,replist,'B*nc'));
bio_files = {bio_files(:).name}';

% Core files
core_files = dir(fullfile(DataDir,replist,'*nc'));
core_files = {core_files(:).name}';
iscore = ~cellfun('isempty',regexp(core_files,'^R[0-9]*')) | ~cellfun('isempty',regexp(core_files,'^D[0-9]*'));
core_files = core_files(iscore);

% Sort file names by their cycle number : bio_files

[~,num_cy_numB] = strtok(bio_files,'_');
num_cy_numB = strrep(strrep(num_cy_numB,'_',''),'.nc','');
%--EB
[~,isort] = sort(num_cy_numB);
bio_files = bio_files(isort);
%--EB
num_cy_numB = num_cy_numB(isort);
%--EB
[~,num_cy_numC] = strtok(core_files,'_');
num_cy_numC = strrep(strrep(num_cy_numC,'_',''),'.nc','');
%--EB
[~,isort] = sort(num_cy_numC);
core_files = core_files(isort);
%--EB
num_cy_numC = num_cy_numC(isort);
%--EB

% =========================================================================
%% Read Bio_files and Core_files, make the mutliprof structure
% =========================================================================
%-- EB, 12/10/2017
Smult.n_prof = 0;
% Smult = [];
%-- EB
nlevels = [];    % to keep the individual n_levels information
DIMmult = [];
% loop over Bio mono-profiles
for k = 1:length(bio_files)
    % Read Bio-file
    [SmonoB,DIMB,~] = NCR_file(fullfile(DataDir,replist,bio_files{k}));
    
    % Read the corresponding core_file
    isOkCore = strcmp(num_cy_numB(k),num_cy_numC);
    [SmonoC,DIMC,~] = NCR_file(fullfile(DataDir,replist,core_files{isOkCore}));
    % Read Parameter_data_mode : if blank found, replace by R. List file in
    % logfile
    isblank = isspace(SmonoB.parameter_data_mode.data);
    SmonoB.parameter_data_mode.data(isblank) = 'R';
    
    % log file for blank parameter 
%     if sum(sum(isblank))>0
%         %if any(isblank)
%         %--EB
%         for nbprof = 1:DIMB.n_prof.dimlength
%             if any(isblank(nbprof,:))
%                 param = find(isblank(nbprof,:));
%                 if exist('logId','var')
%                     fprintf(logId,'%s; %s; %s\n',bio_files{k},num2str(nbprof),num2str(param));
%                 else
%                     logId = fopen(fullfile(LogDir,...
%                         sprintf('DOXY_deblank_parameter_data_mode_%d.log',wmo)),'w');
%                     fprintf(logId,'File; profile index; parameter index\n');
%                     fprintf(logId,'%s; %s; %s\n',bio_files{k},num2str(nbprof),num2str(param));
%                 end
%             end
%         end
%     end
        
    SmonoB.n_prof = DIMB.n_prof.dimlength;
    SmonoC.n_prof = DIMC.n_prof.dimlength;
    
    % ---------------------------------------------------------------------
    %% Add PTS from core
    % - TEMP and PSAL from the CTD (C-file) will be the LOCODOX reference
    % to make the correction (compute the DENS)
    % - PRES:
    %   - in the ARGO format, the variable PRES is common to the C-file
    % and the B-file, i.e. the PRES of the primary profile in the C-file is
    % the same than the PRES in the primary profile in the B-file (resp.
    % for all the other VSS).
    %   - Moreover, the variable PRES is specific to a profile: for
    %   example, the PRES of the primary profile is different from the PRES
    %   in the secondary profile (both in C-file and B-file).
    %   - PRES_ADJUSTED exists only in the C-file. If exists, it should be
    %   used as the reference pressure.
    % ---------------------------------------------------------------------
    SmonoB.pres = SmonoC.pres;
    SmonoB.pres_qc = SmonoC.pres_qc;
    SmonoB.pres_adjusted = SmonoC.pres_adjusted;
    SmonoB.pres_adjusted_qc = SmonoC.pres_adjusted_qc;
    SmonoB.temp = SmonoC.temp;
    SmonoB.temp_qc = SmonoC.temp_qc;
    SmonoB.temp_adjusted = SmonoC.temp_adjusted;
    SmonoB.temp_adjusted_qc = SmonoC.temp_adjusted_qc;
    SmonoB.psal = SmonoC.psal;
    SmonoB.psal_qc = SmonoC.psal_qc;
    SmonoB.psal_adjusted = SmonoC.psal_adjusted;
    SmonoB.psal_adjusted_qc = SmonoC.psal_adjusted_qc;
    
    % Save data_mode from Core Files
    SmonoB.data_mode_CF = SmonoC.data_mode;
    
    % station parameter : keep the index where DOXY is present in the profil
    SmonoB.hereDoxy.name = 'hereDoxy';
    SmonoB.hereDoxy.dim = {'N_PROF'};
    SmonoB.hereDoxy.data = false(1,SmonoB.n_prof);
    SmonoB.hereDoxy.long_name = 'Index of the presence of the variable DOXY';
    SmonoB.hereDoxy.standard_name = 'DOXY presence index';
    isok = false(1,SmonoB.n_prof);
    for i = 1:SmonoB.n_prof
        test = strcmp(cellstr(deblank(squeeze(SmonoB.station_parameters.data(i,:,:)))),'DOXY');
        if any(test)
            isok(i) = true;
        end
    end    
    
    if sum(isok) > 2
        % More than two differents profiles carrying DOXY
        % example : two secondary profiles with DOXY, and one nearSurface
        % profile.        
        VSS_Name = regexp(cellstr(SmonoB.vertical_sampling_scheme.data(isok,:)),'.*(?= sampling)','match');
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
    
    SmonoB.hereDoxy.data = isok;
    
    % keep information in a log
    if exist('logIdVSS','var')
        fprintf(logIdVSS,'%s; %d; %s\n',bio_files{k},DIMB.n_prof.dimlength,num2str(isok));
    else
        logIdVSS = fopen(fullfile(LogDir,...
            sprintf('DOXY_VSS_choice_%d.log',wmo)),'w');
        fprintf(logIdVSS,'File; nb of VSS; VSS logical index\n');
        fprintf(logIdVSS,'%s; %d; %s\n',bio_files{k},DIMB.n_prof.dimlength,num2str(isok));
    end

    % Concatenate monoprofile by monoprofile
    %--EB, 12/10/2017
    %   if SmonoB.n_prof ~= 0 && ~isempty(Smult)
    if SmonoB.n_prof ~= 0
        %-- EB
        if Smult.n_prof ~= 0
            [Smult,DIMmult] = cat_profile_dim(Smult,SmonoB,DIMmult,DIMB,'N_PROF');
        else
            Smult = SmonoB;
            DIMmult = DIMB;
        end
        nlevels = [nlevels;repmat(DIMB.n_levels.dimlength,DIMB.n_prof.dimlength,1)];
        %-- EB, 12/10/2017
        %    elseif isempty(Smult)
        %        Smult.n_prof = 0;
        %        DIMmult = DIMB;
    elseif isempty(Smult)
        DIMmult = DIMB;
        % -- EB
    end       
end

if exist('logId','var')
    fclose(logId);
end
if exist('logIdVSS','var')
    fclose(logIdVSS);
end

Smult.nlevels.dim = {'N_PROF'};
Smult.nlevels.data = nlevels;

end
