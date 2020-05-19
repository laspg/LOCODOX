% DOXY_WOA_read reads the WOA climatology file and convert in new unit
%
% SYNTAX
% [WOA] = DOXY_WOA_read(WOAfile, unit)
%
% DESCRIPTION
% DOXY_WOA_read reads the WOA climatology file and stocks it in a structure.
% The DOXY parameter is converted in a new unit.
%
% INPUT
%   WOAfile (string)     full name of the NetCDF WOA climatology file to be
%                        read
%
%   unit (string)        unit in which the doxy parameter has to be
%                        converted
%                        Example: 'mumol/L'
%
% OUTPUT
%   WOA (structure)      filled with the variable and its attributes found
%                        in the NetCDF WOA climatology file
%                        Example:
%                         WOA = 
% 
%                                dimorder: 'C'
%                                latitude: [1x1 struct]
%                               longitude: [1x1 struct]
%                                   depth: [1x1 struct]
%                                    time: [1x1 struct]
%                                 doxywoa: [1x1 struct]
%                                 psatwoa: [1x1 struct]
%                                psal_woa: [1x1 struct]
%                                temp_woa: [1x1 struct]
%                                   ..........
%
% CALL :
%   NCR_file, DOXY_convert, regexp, O2ctoO2s
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 26/11/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:


function [WOA] = DOXY_WOA_read(WOAfile, unit)


% =====================================================================
%% Read the WOA monthly climatology (levels : 0 to 5500m ; size = (360,180,33,12))
% =====================================================================
[WOA,~,~]=NCR_file(WOAfile);
WOA.an_dens_woa = WOA.density;
WOA.an_dens_woa.data = WOA.density.data-1000;

% =====================================================================
%% Convert in the choosen unit
% =====================================================================
WOA.doxywoa_CV = WOA.doxywoa;
WOA.doxywoa_CV.units = unit;
% replace unit in long_name
unUnit = regexp(WOA.doxywoa_CV.long_name,'([a-z]*/[a-z]*)','split');
WOA.doxywoa_CV.long_name = sprintf([unUnit{1} '%s' unUnit{2}],unit);

% convert data
WOA.doxywoa_CV.data = DOXY_convert(WOA.doxywoa.data,WOA.doxywoa.units,unit,WOA.an_dens_woa.data);

% =====================================================================
%% Compute psatwoa taking into account the pressure effect
% =====================================================================
% In WOA, the psatwoa provided has been computed using the seawater
% formulas (Weiss, R.F., 1970): sw_sat02 + (doxywoa/sat)*100. The new recommendation leads to
% re-compute PSAT using the Gordon and Garcia, 1992 formulas.
% change long_name of original psatwoa
psatLongName = WOA.psatwoa.long_name;
WOA.psatwoa.long_name = [psatLongName ', Weiss 1970 - no pressure effect'];
[WOA.psatwoa_original] = WOA.psatwoa;

% Compute with Goordon and Garcia, 1992 : no pressure effet
WOA.psatwoa.long_name = [psatLongName ', Gordon and Garcia, 1992 - no pressure effect : P = 0'];
tmpDoxy = DOXY_convert(WOA.doxywoa_CV.data,unit,'mumol/L',WOA.an_dens_woa.data);
WOA.psatwoa.data = O2ctoO2s(tmpDoxy,WOA.temp_woa.data,WOA.psal_woa.data);

% compute new psatwoa, with P = preswoa
WOA.psatwoa.long_name = [psatLongName ', Gordon and Garcia, 1992 - pressure effect : P = preswoa'];
WOA.psatwoa_preswoa = WOA.psatwoa;
tmpDoxy = DOXY_convert(WOA.doxywoa_CV.data,unit,'mumol/L',WOA.an_dens_woa.data);
WOA.psatwoa_preswoa.data = O2ctoO2s(tmpDoxy,WOA.temp_woa.data,WOA.psal_woa.data,-WOA.preswoa.data);
WOA.psatwoa_preswoa.long_name = [psatLongName ', pressure effect taking into account : P = preswoa'];
