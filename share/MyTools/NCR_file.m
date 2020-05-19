%  NCR_file read a netcdf file, for all variables of a
%  selection of variables.
%
% SYNTAX
% [S,Dim,Globatt] = NCR_file(Ficname,VarIN,verbose)
%
% DESCRIPTION
% NCR_file read all variables or a selection of variables,
% its attributes, the dimension and the global attributes of a NetCDF file.
%
% INPUT
%   Ficname (string)    name of the file (example : 'filetoread.nc')
%
%   VarIN (cell)        (OPTIONNAL) name of variables to read.
%                       By Default, read all the variables.
%                       Example : VarIN = {'PRES','TEMP'};
%
%   verbose (scalar)    (OPTIONNAL) 0 (default),  1 print to screen.
%
% OUTPUT
%   S (structure)       Contains all information found in the NetCDF file
%                       for the variables (including attributes).
%                       Example :
%                           S.temp.name = 'TEMPERATURE' name of the variable
%                           S.temp.dim  = {'N_PROF','N_LEVEL'}
%                           S.temp.data =  n_prof x n_level value of temperature
%
%   DIMS (structure)    Contains all the dimensions information found in
%                       the NetCDF file.
%                       Example :
%                           Dim.n_prof.name = 'N_PROF'     name of the dimension
%                           Dim.n_prof.dimlength = 88      length of the dimension
%
%   Globatt (structure) Contains all the global attributes found in the
%                       NetCDF file
%
% CALL :
%   check_FirstDimArray, check_LastDimArray, fillValue_2_nan
%
% SEE ALSO
%   DOXY_corr_main

% -----------------------------
% REMARKS
% Matlab (or fortran F) and C manage the dimensions in different ways:
%      Matlab: X(x,y,z,t), unlimitted is the latter dimension
%           C: X(t,z,y,x), unlimitted is the first dimension
% If you do a ncdump of a NetCDF file, the dimensions are sorted like with
% C. Reading a NetCDF file, the Matlab NetCDF library sorted the dimensions in the reverse
% direction. Here, loaded variables and dimensions are resorted in the "C"
% way, to be identical to what is obtained with ncdump.
% -----------------------------

% -----------------------------
% INFO: The variable type is identified by a number, as followed.
%       NC_BYTE .... 1
%       NC_CHAR .... 2
%       NC_SHORT ... 3
%       NC_INT ..... 4
%       NC_FLOAT ... 5
%       NC_DOUBLE .. 6
% -----------------------------

function [S,Dim,Globatt] = NCR_file(Ficname,VarIN,verbose)


% =========================================================================
%% Initialisation
% =========================================================================

if nargin==1
    verbose = 0;
    VarIN =  {};
end

if nargin==2
    if iscell(VarIN)
        verbose = 0;
    else
        verbose = VarIN;
        VarIN = {};
    end
end

% Possible name for the FillValue attribute (to overcome the files that do not comply the CF conventions for this attribute).
poss_fillval_name = {'FillValue','FillValue_','_FillValue','_fillvalue','fillvalue','missing_value'};


% =========================================================================
%% NetCDF reading
% =========================================================================

% Open the NetCDF file
ncid = netcdf.open(Ficname, 'NC_NOWRITE');

% Get the file information
%   Nbdims      : number of dimensions
%   Nbfields    : number of variables
%   Nbglob      : number of global attributes
%   theRecdimID : Id of the unlimitted dimension
[Nbdims, Nbfields, Nbglob, theRecdimID] = netcdf.inq(ncid);

% Get all the Dimensions
Dim = [];
for kd=1:Nbdims
    [namedim, dimlen] = netcdf.inqDim(ncid,kd-1);
    onedim =lower(namedim);
    Dim.(onedim).name = namedim;
    Dim.(onedim).dimlength = dimlen;
end

% Get all the global attributes
Globatt = [];
if nargout==3
    for kd=1:Nbglob
        nameatt = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),kd-1);
        Globatt.(nameatt).att = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),nameatt);
        Globatt.(nameatt).name = nameatt;
    end
end

% -------------------------------------------------------------------------
% Select the variables to be read, respect to VarIN.
% -------------------------------------------------------------------------
% Read all the variable names in the NetCDF file
kread = 1:Nbfields;
allvar = cell(1,Nbfields);
for k=kread
    [varname] = netcdf.inqVar(ncid,k-1);
    allvar{k} = varname;
end

% Compare to the variables names in VarIN
if ~isempty(VarIN)
    [~,ia,~] = intersect(allvar,VarIN);
    kread = kread(ia);
end
% 
% if isempty(kread)
%     fprintf('NCR_file: Asked variables doesn''t exist\n');
%     return
% end

% -------------------------------------------------------------------------
% Set the variables dimensions order (Matlab way or C way)
%   S.dimorder='C' : dimensions are sorted in the C way (i.e unlimitted first dimension)
%   S.dimorder='F' : dimensions are sorted in the Matlab or Fortran way (i.e unlimitted last dimension)
% -------------------------------------------------------------------------

S.dimorder = 'C';  % option prefered because dimensions in S are sorted like in the ncdump screen print.
% S.dimorder='F';
% -------------------------------------------------------------------------
% READ the variables and fill the structure S for each variable
% -------------------------------------------------------------------------
count_alldim = 0;
nbmaxdim = 100;
alldim = cell(1,nbmaxdim);

% Loop over the variable to be read in the NetCDF files
for k=kread    
    % Get the variable information
    % ---------------------------------------------------------------------
    varid = netcdf.inqVarID(ncid,allvar{k});
    [VarName,xtype,dimids,natts] = netcdf.inqVar(ncid,varid);
    oneChamp=lower(VarName);
        
    % Intialise the structure S for the variable
    % ---------------------------------------------------------------------
    S.(oneChamp).name=VarName;
    S.(oneChamp).dim=[];
    S.(oneChamp).data=[];

    vecdim=NaN(1,length(dimids)+1);
    %vecdim=[];
    
    % Get the variable attributes
    % ---------------------------------------------------------------------
    clear name_att
    name_att = cell(1,natts);
    for iat=1:natts
        % Attribute name
        attid = iat-1;
        name_att{iat} = netcdf.inqAttName(ncid,varid,attid);
        
        % Attribute value and case of the _FillValue attribute
        isthefil=0;
        for ipos=1:length(poss_fillval_name)
            if strcmp(name_att{iat},poss_fillval_name{ipos});
                isthefil=1;
            end
        end
        if ~isthefil
            S.(oneChamp).(name_att{iat}) = netcdf.getAtt(ncid,varid,name_att{iat} );
        else
            S.(oneChamp).FillValue_ = netcdf.getAtt(ncid,varid,name_att{iat} );
        end
    end
    
    if ~isfield(S.(oneChamp),'long_name')
        S.(oneChamp).long_name = [];
    end
    if ~isfield(S.(oneChamp),'units')
        S.(oneChamp).units = [];
    end
    
    % Get the variable dimensions
    % ---------------------------------------------------------------------
    S.(oneChamp).dim = '';
    count_vardim = 0;
    % dimension permutation to have S.dimorder='C'
    if strcmp(S.dimorder,'C')
        theloop = length(dimids):-1:1;
    else
        theloop = 1:length(dimids);
    end
    
    for idim = theloop
        [thedimname, dimlen] = netcdf.inqDim(ncid,dimids(idim));        
        count_alldim = count_alldim+1;
        count_vardim = count_vardim+1;
        S.(oneChamp).dim{count_vardim} = thedimname;
        vecdim(count_vardim) = dimlen;
        alldim{count_alldim} = thedimname;
    end
    
    % Get the variable type
    % ---------------------------------------------------------------------
    S.(oneChamp).type = xtype;
    
    % Get the variable data and reshape along dimensions
    % ---------------------------------------------------------------------
    %      if sum(vecdim==0)==0
    if all(vecdim~=0)
        % Get the data value
        S.(oneChamp).data = netcdf.getVar(ncid,varid);

        % Permute to sort the data in the 'C' way
        if length(vecdim)>1
            if strcmp(S.dimorder,'C')
                S.(oneChamp).data = permute(S.(oneChamp).data,length(size(S.(oneChamp).data)):-1:1);
            end
        end
        
        % Reshape the S.(oneChamp).data to coincide with vecdim        
        if any(vecdim==1)        % presence of singleton dimension
            vecdim(end)= 1;      % prevent from a bug for scalar 
            S.(oneChamp).data = reshape(S.(oneChamp).data,vecdim);
        end
    end
    
    % Get the variable units
    % ---------------------------------------------------------------------
    if ~isfield(S.(oneChamp),'units')
        S.(oneChamp).units=[];
    end
    if verbose
        disp(['............',oneChamp])
    end
end

% =========================================================================
%% Finalisation
% =========================================================================
% Clean not used dimensions for variables that are not retrieved
namedim = fieldnames(Dim);
alldim = alldim(~cellfun('isempty',alldim));
rmdim = setdiff(namedim,lower(alldim));
for jdim=1:length(rmdim)
    Dim = rmfield(Dim,rmdim{jdim});
end

if ~isempty(theRecdimID)
    if theRecdimID>0
        if ~isempty(rmdim)
            recdimname = netcdf.inqDim(ncid, theRecdimID);
%             if strcmpi(recdimname,rmdim)==0
            if all(strcmpi(recdimname,rmdim))
                S.recdim = recdimname;
            end
        else
            S.recdim = netcdf.inqDim(ncid, theRecdimID);
        end
    end
end

% Additional fields
S.obj='ObsInSitu';
if isfield(VarIN,'obj')
    S.obj=[S.obj '/' VarIN.obj];
end

S.fillisnan = 0;

% Close the NetCDF file
netcdf.close(ncid)

