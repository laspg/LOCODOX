%  NCW_file create a netcdf file containing variables in VARS
%
% SYNTAX
% NCW_file(VARS,DIMS,ficout,CONFIG,verbose)
%
% DESCRIPTION
% NCW_file create a NetCDF file from the variables and its
% attributes find in the VARS structure. The DIMS structures is also used.
% The global attributes is filled in the NetCDF file if the CONFIG
% structure is given in argument.
%
% INPUT
%   VARS (structure)    Contains the variables definition. Required fields
%                       for each variable : name (string), dim (cell array
%                       of string), data (array)
%                       Example :
%                           VARS.temp.name = 'TEMPERATURE' name of the variable
%                           VARS.temp.dim  = {'N_PROF','N_LEVEL'}
%                           VARS.temp.data =  n_prof x n_level value of temperature
%
%   DIMS (structure)    Contains the dimensions definition. Required fields
%                       for each dim : name (string), dimlength (scalar)
%                       Example :
%                           Dim.n_prof.name = 'N_PROF'     name of the dimension
%                           Dim.n_prof.dimlength = 88      length of the dimension
%
%   ficout (string)     output file name
%
%	CONFIG (structure) (OPTIONNAL) Contains global attribute you want to
%                      put in the netcdf file. Required fields : name
%                      (string)
%
%  verbose (scalar)    (OPTIONNAL) 0 (default),  1 print to screen.
%
% OUTPUT
%    none
%
% CALL :
%   check_FirstDimArray_is, check_LastDimArray_is, replace_nan_byfill
%
% SEE ALSO
%   read_netcdf_allthefile
%

% HISTORY
%   $created: /05/2007 $author: Cecile Cabanes, LPO, CNRS
%   $Revision: version $Date: $author:
%       v?  /10/2008    Cecile Cabanes, LPO, CNRS
%       v?  /02/2009    Cecile Cabanes, LPO, CNRS
%       v?  /04/2010    Cecile Cabanes, LPO, CNRS
%                       adapted to Matlab up to R2008b
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       adapted to be shared to the O2 community

% -----------------------------
% REMARKS
% Matlab (or fortran F) and C manage the dimensions in different ways:
%      Matlab: X(x,y,z,t), unlimitted is the latter dimension
%           C: X(t,z,y,x), unlimitted is the first dimension
% If you do a ncdump of a NetCDF file, the dimensions are sorted like with
% C. Reading a NetCDF file, the Matlab NetCDF library sorted the dimensions in the reverse
% direction. So, to write a NetCDF file, the Matlab NetCDF library need a
% variable sorted like in Fortran (unlimitted is the latter dimension).
% -----------------------------


function NCW_file(VARS,DIMS,ficout,CONFIG,verbose)


% =========================================================================
%% Initialisation
% =========================================================================

if nargin <= 3
    CONFIG = [];
    verbose = 0;
end

if nargin == 4
    if isstruct(CONFIG)
        verbose = 0;
    else
        verbose = CONFIG;
        CONFIG = [];
    end
end

% -------------------------------------------------------------------------
% Check the variables dimensions order (Matlab way or C way)
%   VARS.dimorder='C' : dimensions are sorted in the C way (i.e unlimitted first dimension)
%   VARS.dimorder='F' : dimensions are sorted in the Matlab or Fortran way (i.e unlimitted last dimension)
% -------------------------------------------------------------------------

VARS.obj = 'ObsInSitu';
% find which dimension is the unlimitted one
REC = [];
if isstruct(VARS)
    if isfield(VARS,'recdim')
        if ~isempty(VARS.recdim)
            REC = VARS.recdim;
            %REC_length = DIMS.(lower(REC)).dimlength;
        end
    end
end

% find the dimensions ordering way
if ~isfield(VARS,'dimorder') % si rien n'est precise
    if ~isempty(REC)
        VARS = check_FirstDimArray(VARS,REC);
    end
    VARS.dimorder = 'C';
end

if strcmp(VARS.dimorder,'F')
    isneedflip = 0;
    if ~isempty(REC)
        VARS = check_LastDimArray(VARS,REC); % check and permute if needed
    end
else
    isneedflip = 1;
    if ~isempty(REC)
        VARS = check_FirstDimArray(VARS,REC); % check and permute if needed
    end
end

% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------
% NaN replaced by FillValue
VARS = nan_2_fillValue(VARS);

% Change structures into cell arrays
if isstruct(VARS)
    VARS = struct2cell(VARS);
end
if isstruct(DIMS)
    DIMS = struct2cell(DIMS);
end

NBDIMS = length(DIMS);
NBVARS = length(VARS);

disp(['writing ' num2str(NBVARS), ' variable(s)  in the file....'])
disp(ficout)


% =========================================================================
%% NetCDF writing
% =========================================================================
% Open the NetCDF file
ncid = netcdf.create(ficout, 'CLOBBER');

% -------------------------------------------------------------------------
% WRITE the Global attributes:
% -------------------------------------------------------------------------
% Check in CONFIG if the fields exist
if ~isempty(CONFIG)
    if isstruct(CONFIG)
        CONFIG = struct2cell(CONFIG);
    end
    NBATT = length(CONFIG);
    for k = 1:NBATT
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),CONFIG{k}.name, CONFIG{k}.att)
        
        if verbose
            %            disp(['Global attribute : ' expre])
        end
    end
end

% look for the memory space needed to store VARS -> memo.bytes
% memo=whos('VARS');

% -------------------------------------------------------------------------
% WRITE the Dimensions
% -------------------------------------------------------------------------
dimID = NaN(1,NBDIMS);
for k = 1:NBDIMS
    if ~isfield(DIMS{k},'name')
        disp(' ')
        keyboard
        disp(['Dimension : ' DIMS{k}])
        error ('ERROR:  Should give a name  in Dim.name ')
    end
    if ~isfield(DIMS{k},'dimlength')
        disp(' ')
        disp(['Dimension : ' DIMS{k}])
        error ('ERROR: Should give a length  in Dim.dimlength ')
    end
    
    if isfield(DIMS{k},'dimlength')
        if ~isempty(REC)
            if strcmp(DIMS{k}.name,REC)   % si c'est la dimension unlimitted
                dimID(k) = netcdf.defDim(ncid,DIMS{k}.name,netcdf.getConstant('NC_UNLIMITED'));
            else
                dimID(k) = netcdf.defDim(ncid,DIMS{k}.name,DIMS{k}.dimlength);
            end
        else
            dimID(k) = netcdf.defDim(ncid,DIMS{k}.name,DIMS{k}.dimlength);
        end
    end
    
    if verbose
        if isfield(DIMS{k},'name')
            disp(' ')
            disp(['Dimension............',DIMS{k}.name])
        end
    end
end


% -------------------------------------------------------------------------
% DEFINITION of the Variables and its attributes
% -------------------------------------------------------------------------
varid = NaN(1,NBVARS);
for k = 1:NBVARS
    missval = [];
    if verbose
        if isfield(VARS{k},'name')
            disp(' ')
            disp(['............',VARS{k}.name])
        end
    end
    
    % Get the needed information for variables
    % ---------------------------------------------------------------------
    if isfield(VARS{k},'data')   % if it is a variable (.data field)        
        if ~isfield(VARS{k},'name')
            disp([ 'Variable ............:',VARS{k}])
            error ('ERROR: should give a name  in Var.name to all variables')
        end
        if ~isfield(VARS{k},'dim')
            disp([ 'Variable ............:',VARS{k}.name])
            error ('ERROR: should give a dim  in Var.dim to all variables')
        end
        
        % the dimension name of the variable
        extract_dim = VARS{k}.dim;
        
        % screen output : look for the FillValue_ attribute if non standard
        if isfield(CONFIG,'missval')
            missval = CONFIG.missval;
        end
        if isfield(VARS{k},'FillValue_')
            missval = VARS{k}.FillValue_;
        end
        if isfield(VARS{k},'fillval')
            missval = VARS{k}.fillval;
        end
        
        % check the data type
        if isfloat(VARS{k}.data) || isinteger(VARS{k}.data)            
            % check for infinite value
            if any(isinf(VARS{k}.data))
%             if isempty(find(isinf(VARS{k}.data)))==0
                VARS{k}.data(isinf(VARS{k}.data)) = NaN;
                warning(['infinite value found in ' VARS{k}.name ' :stored as missing_value in the netcdf file']);
            end
%             % Missing .missval field : created
%             if isempty(VARI{k}.data(isnan(VARI{k}.data)))==0
%                 if isempty(missval)
%                     missval=99999;
%                 end
%             end
            
            dataTYP = netcdf.getConstant('float');
            if isa(VARS{k}.data,'double'); dataTYP = netcdf.getConstant('double');end;
            if isinteger(VARS{k}.data); dataTYP = netcdf.getConstant('int');end;            
        end
        
        if ischar(VARS{k}.data)
            dataTYP = netcdf.getConstant('char');
        end
        if isfield(VARS{k},'type')
            dataTYP = VARS{k}.type;
        end
        
        % DEFINITION of the Variables
        % -----------------------------------------------------------------
        
        % Get te ID of the variable dimensions
        siz = length(extract_dim);
        clear alldimid
        % dimension permutation to have VARS.dimorder='F'            
        if isneedflip==1
            theloop = siz:-1:1; 
        else
            theloop = 1:siz;
        end        
        ij=1;
        alldimid = NaN(1,length(theloop));
        for m = theloop
            alldimid(ij) = netcdf.inqDimID(ncid,extract_dim{m});
            ij=ij+1;
        end
        
        % Define the variables
        %fprintf('nvar_def = %d\n',k);
        varid(k) = netcdf.defVar(ncid,VARS{k}.name,dataTYP,alldimid);
        
        % WRITE the variable attributes
        % -----------------------------------------------------------------
        attnames = fieldnames(VARS{k});
        nbatt = length(attnames);
        for iatt=1:nbatt
            if verbose==1
                disp('**********************')
                disp(['attributes: ' attnames{iatt}])
                disp('-----------------')
            end
            theattname = attnames{iatt};
            theattname = strrep(theattname,'FillValue_','_FillValue');
            
            % Not used in the NetCDF file : the fields 'data', 'dim', 'type'  'name' and 'ischar2num'
%             if ~(strcmp(theattname,'data')==1|strcmp(theattname,'dim')==1|strcmp(theattname,'type')==1|strcmp(theattname,'name')==1|strcmp(theattname,'ischar2num')==1)

            if ~(strcmp(theattname,'data') || strcmp(theattname,'dim') || ...
                 strcmp(theattname,'type') || strcmp(theattname,'name') || ...
                 strcmp(theattname,'ischar2num'))
                thevalatt = VARS{k}.(attnames{iatt});
                if verbose==1
                    disp(thevalatt);
                end                
                if ~isempty(thevalatt)
                    netcdf.putAtt(ncid,varid(k),theattname,thevalatt);
                end
            end
        end 
    end    
end

% -------------------------------------------------------------------------
% WRITE the variables
% -------------------------------------------------------------------------
netcdf.endDef(ncid);
for k = 1:NBVARS
    if isfield(VARS{k},'data')
        if ~isempty(VARS{k}.data)
%             disp( VARS{k}.name)
%             VARS{k}.dim
%             if strcmp(VARS{k}.name,'DATA_TYPE')
%                 keyboard
%             end
            ndim_var = length(VARS{k}.dim);
            vecdim = size(VARS{k}.data);            
            %if ndim_var>1
            if isneedflip
                VARS{k}.data = permute(VARS{k}.data,length(size(VARS{k}.data)):-1:1);
                vecdim = fliplr(vecdim);
            end
            %end
            
            if ndim_var>1
                start_var = zeros(ndim_var,1);
                count_var = vecdim;
            else
                start_var = 0;
                count_var = length(VARS{k}.data);
            end
            %fprintf('VAR numero %d\n',k);
            %fprintf('nvar = %d\n',k);
            netcdf.putVar(ncid,varid(k),start_var,count_var,VARS{k}.data);
        end
    end
end

% Close the NetCDF file
netcdf.close(ncid)