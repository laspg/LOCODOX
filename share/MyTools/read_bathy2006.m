
function [topo2,lat2,lon2] = read_bathy2006(region,CONFIG)

% Function read_bathy  Read bathymetry data from Sandwell Database
%
%      [topo2,lat2,lon2] = read_bathy2006(ymin ymax xmin xmax)
%
% This function reads etopo2 from Smith and Sandwell 2006


ymin=region(1);
ymax=region(2);
xmin=region(3);
xmax=region(4);

ftopo2 =  [CONFIG.BathyDataDir,'ETOPO2v2c_f4.nc'];

ncid = netcdf.open(ftopo2, 'NC_NOWRITE');

varid = netcdf.inqVarID(ncid,'x');
lon2= netcdf.getVar(ncid,varid)  ;

varid = netcdf.inqVarID(ncid,'y');
lat2= netcdf.getVar(ncid,varid)  ;

ji1 = find(lon2 > xmin, 1 ); 
ji2 = find(lon2 > xmax, 1 )-1;
if isempty(ji2) 
    ji2=size(lon2,1);
end
jj1 = find(lat2 > ymin, 1 ); 

jj2 = find(lat2 > ymax, 1 )-1;
if isempty(jj2) 
    jj2=size(lat2,1);
end

varid = netcdf.inqVarID(ncid,'z');
topo2 = netcdf.getVar(ncid,varid,[(ji1-1) (jj1-1) ],[(ji2-ji1+1) (jj2-jj1+1)]);
netcdf.close(ncid);

topo2=double(topo2');
lat2=double(lat2(jj1:1:jj2));
lon2=double(lon2(ji1:1:ji2));
