% ------------------------------------------------------------------------------
%c
%c Conv_deg       - Conversion degres decimaux --> deg-min dec
%                                                   deg-min-sec dec
%c
%c      function [NS, ylat1, ylat2, EW, xlon1, xlon2] = Conv_deg(ylat, xlon);
%c 
% ------------------------------------------------------------------------------
%  version:
%  --------
%  1.01                                               26/06/97   F.Gaillard
%
% ------------------------------------------------------------------------------

function [NS, ylat1, ylat2, EW, xlon1, xlon2] = Conv_deg(ylat, xlon);


nn = length(ylat);
NS = [];
EW = [];
for i = 1:nn
   NS = [NS; 'N'];
   EW = [EW; 'E'];
end
is_sud = find(ylat<0);
NS(is_sud) = 'S'*ones(size(is_sud));
is_ouest = find(xlon<0);
EW(is_ouest) = 'W'*ones(size(is_ouest));

lat_dr = abs(ylat);
lat_di = floor(lat_dr);
lat_mr = (lat_dr - lat_di)*60;
lat_mi = floor(lat_mr);
lat_sr = (lat_mr - lat_mi)*60;
ylat1 = [lat_di lat_mr];
ylat2 = [lat_di lat_mi lat_sr];

lon_dr = abs(xlon);
lon_di = floor(lon_dr);
lon_mr = (lon_dr - lon_di)*60;
lon_mi = floor(lon_mr);
lon_sr = (lon_mr - lon_mi)*60;
xlon1 = [lon_di lon_mr];
xlon2 = [lon_di lon_mi lon_sr];
