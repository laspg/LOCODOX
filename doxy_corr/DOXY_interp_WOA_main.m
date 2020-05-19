% DOXY_interp_WOA_main interpolate WOA data over argo horizontal and
% vertical grid.
%
% SYNTAX
% [WOA, Work] = DOXY_interp_WOA_main(maskFile, varMask, argo, WOA,argoWork, Work)
%
% DESCRIPTION
% DOXY_interp_WOA_main interpolate WOA data over argo horizontal and
% vertical grid. The vertical grid could be pressure (interpP) and over
% density (interpD). The working structure is completed. Attention is given
% to near-shore data.
%
% INPUT
%   maskFile (string)   Pathname of the mask file
%
%   varMask (string)    The variable name of the mask file useful for this
%                       program.
%                       Example : varMask = 'landsea';
%
%   Work (structure)    Float working structure, issued and computed from
%                       the argo data and the climatology data (WOA)
%                           Example:
%                           Work = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                           
%                           Work.DOXY_WOA=
%                           name: 'DOXY_ADJUSTED'
%                            dim: {'N_PROF'  'N_LEVELS'}
%                           data: [85x120 double]
%                      long_name: 'DISSOLVED OXYGEN adjusted and converted'
%                  standard_name: 'moles_of_oxygen_per_unit_mass_in_sea_water'
%                     FillValue_: 99999
%                          units: 'mumol/kg'
%                      valid_min: 0
%                      valid_max: 650
%                       C_format: '%9.3f'
%                 FORTRAN_format: 'F9.3'
%                     resolution: 0.0010
%                           type: 5
%              
%   argowork (structure)  Float working structure issued from the argo data
%                         argoWork = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                                       sat: [1x1 struct]
%                                      psat: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%
%   argo (structure)      Argo float structure (get directly from data).
%                         argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%
%                               argo.pres =
%                                       name: 'DOXY'
%                                        dim: {'N_PROF'  'N_LEVELS'}
%                                       data: [85x120 single]
%                                  long_name: 'DISSOLVED OXYGEN'
%                              standard_name: 'moles_of_oxygen_per_unit_mass_in_sea_water'
%                                 FillValue_: 99999
%                                      units: 'micromole/kg'
%                                  valid_min: 0
%                                  valid_max: 650
%                                   C_format: '%9.3f'
%                             FORTRAN_format: 'F9.3'
%                                 resolution: 0.0010
%                                       type: 5
%
%   WOA (structure)      filled with the variable and its attributes found
%                        in the NetCDF WOA climatology file, and variables
%                        computed for the DOXY correction.
%                        Example:
%                         WOA = 
%                                dimorder: 'C'
%                                latitude: [1x1 struct]
%                               longitude: [1x1 struct]
%                                   depth: [1x1 struct]
%                                    time: [1x1 struct]
%                                 doxywoa: [1x1 struct]
%                                 psatwoa: [1x1 struct]
%                                psal_woa: [1x1 struct]
%                                temp_woa: [1x1 struct]
%                                 density: [1x1 struct]
%                              doxywoa_CV: [1x1 struct]
%                                   ..........
%
% OUTPUT
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA).
%                        Example:
%                        Work = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%
%   WOA (structure)      filled with the variable and its attributes found
%                        in the NetCDF WOA climatology file, and variables
%                        computed for the DOXY correction.
%                        Example:
%                         WOA = 
%                                dimorder: 'C'
%                                latitude: [1x1 struct]
%                               longitude: [1x1 struct]
%                                   depth: [1x1 struct]
%                                    time: [1x1 struct]
%                                 doxywoa: [1x1 struct]
%                                 psatwoa: [1x1 struct]
%                                psal_woa: [1x1 struct]
%                                temp_woa: [1x1 struct]
%                                 density: [1x1 struct]
%                              doxywoa_CV: [1x1 struct]
%                                   ..........
%
%
%
% CALL : 
%   DOXY_interp_WOA_2_argo, DOXY_PLOT_interpolation
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 07/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:


function [WOA, Work] = DOXY_interp_WOA_main(maskFile, varMask, argo, WOA, ...
                                            argoWork, Work)

% ========================================================================= 
%% find WOA near argo
% ========================================================================= 
mask = NCR_file(maskFile,varMask,0);% est ce qu'on est proche des côtes ?
cycNum = argo.cycle_number.data;
nearShore = NaN(size(cycNum));

%  Reads climatology coordinates and locates surrounding points
for i=1:length(cycNum)
    %latitude
    i1 = find(WOA.latitude.data <= argo.latitude.data(i), 1, 'last');
    i2 = find(WOA.latitude.data > argo.latitude.data(i), 1, 'first');
    I=[i1,i2]; 
    %longitude
    isneg = argo.longitude.data(i) < 0;
    if any(isneg)
        adjPos = 360;
    else
        adjPos = 0;        
    end
    j1 = find(WOA.longitude.data <= argo.longitude.data(i)+adjPos, 1, 'last');
    j2 = find(WOA.longitude.data > argo.longitude.data(i)+adjPos, 1, 'first');
    J=[j1,j2];

    % est ce qu'on est proche des côtes ?
    % Important pour l'interpolation
    limite = mask.(varMask).data(I,J);
    if min(limite(:)) == 1
        nearShore(i) = 1;
    else
        nearShore(i) = 0;
    end
end

% ========================================================================= 
%% Interpolation : interpolate WOA over argo
% ========================================================================= 
% The interpolation is made over pressure (interpP) 
[WOA] = DOXY_interp_WOA_2_argo(WOA,argoWork,argo,{'preswoa','pres_adjusted'}, nearShore, Work.datat);
Work.('DOXY_WOA_interpP') = WOA.('doxywoa_CV_interpP').data;
Work.('PSAT_WOA_interpP') = WOA.('psatwoa_interpP').data;
Work.('PSAL_WOA_interpP') = WOA.('psal_woa_interpP').data;
Work.('TEMP_WOA_interpP') = WOA.('temp_woa_interpP').data;


% ========================================================================= 
%% Corrections applied on QC <=2
% ========================================================================= 
isnok = str2num(argoWork.doxy_qc.data') > 1; %#ok<ST2NM>
cmpl = {'_interpP','_interpD'};
for i = 1:2
    WOA.(['doxywoa_CV' cmpl{i}]).data(isnok) = NaN;
    WOA.(['psatwoa' cmpl{i}]).data(isnok) = NaN;
end

% ========================================================================= 
%% Control plot
% ========================================================================= 
if Work.makePlot
    DOXY_PLOT_interpolation(Work, argoWork,  1)
end