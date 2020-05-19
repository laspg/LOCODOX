% 
% funcidtion [lat,lon,juld,platform,varargout]=read_ncid_base_hydro(file_name,varargin)
%
% Programme de lecture des fichiers cr???er par le programme base_hydro.m de
% C. Lagadec
%
% Creation : V. Thierry, jan 2009
% 

function [station_id,lat,lon,juld,platform,varargout] = read_nc_base_hydro(file_name,varargin)


ncid = netcdf.open(file_name,'nowrite');
%varid = netcdf.inqVarID(ncid,'CRUISE_NAME');
%platform = netcdf.getVar(ncid,varid)';
varid = netcdf.inqVarID(ncid,'STATION_NUMBER');
statNum = netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'LATITUDE_BEGIN');
lat=netcdf.getVar(ncid,varid);
varid=netcdf.inqVarID(ncid,'LONGITUDE_BEGIN');
lon = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'JULD_BEGIN');
juld = netcdf.getVar(ncid,varid);
clear varid

% platform (refId) = identifiant du fichier Netcdf
%                    compris entre le dernier / et les 5 caracteres
%                    precedant le .devant nc
% ex : file = /home/lpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/cat12_PRES.nc
% on prend cat12 comme refId

ideb=findstr(file_name,'/');
ifin=findstr(file_name,'.');
platform=file_name(ideb(end)+1:ifin(end) - 6);
station_id = cell(size(statNum));
for i=1:length(statNum)
    %tmp = deblank(platform(i,:));
    %tmp(isspace(tmp)) = '';
    %station_id{i} = sprintf('%s_%d',tmp,statNum(i));
    station_id{i} = sprintf('%s_%d',platform,statNum(i));
end

if nargin > 1
  for i=2:nargin
    tempo = varargin{i-1};
    varid = netcdf.inqVarID(ncid,tempo);
    varia = netcdf.getVar(ncid,varid)';
    fillvalue=netcdf.getAtt(ncid,varid,'_FillValue');
    varia(varia==fillvalue) = NaN;
    varargout{i-1} = varia;
  end
end

netcdf.close(ncid)
  


