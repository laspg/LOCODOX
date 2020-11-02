% DOXY_QC undertakes a quality control of doxy data, applying flags.
%
% SYNTAX
% [Swork,corrWork] = DOXY_QC(wmo,Sin,Swork,corrWork)
%
% DESCRIPTION
% DOXY_QC undertakes a small quality control of doxy data and change QC flags
% and values according to the choice made by the user in the configuration
% file.

%
% INPUT
%     wmo (double)          wmo of the float
%
%
%     Sin (structure)       float structure
%                           Example :
%                               Sin =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%
%                               Sin.pres
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
%     Swork (structure)     float working structure, issued and computed from Sin
%                           Example:
%                           Swork = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%   
%     corrWork (structure)  Doxy_corr working structure, issued and
%                           computed from Sin
%                           Example:
%                             corrWork = 
%                                   readme: [1x1 struct]
%                                 DOXY_RAW: [85x120 single]
%                                  CAPTEUR: {'Optode'}
%                                  DOXY_QC: [85x120 char]
%
% OUTPUT
%     Swork (structure)     float working structure, issued and computed
%                           from Sin. Completed with DOXY_QC field.
%
%     corrWork (structure)  Doxy_corr working structure, issued and
%                           computed from Sin. Completed with doxy_qc
%                           field.
%
% CALL :
%
% SEE ALSO
%   DOXY_corr_main, DOXY_get_PTS

% HISTORY
%   $created: 25/11/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%              v2  26/01/2017  Anne Piron, Altran Ouest
%                              argo - ameliorations 2017 phase 1
%                              (take into account user's configuration
%                              parameters for QC data)
%
%              v2.1  09/07/2018 Emilie Brion, Altran Ouest
%                               doxy_adjusted set to NaN if pres_adjusted
%                               unfinite
%                               Swork.doxy_adjusted.data(~isfinite(Swork.pres_adjusted.data)) = NaN;
%                               Correctly manage the "D" directio profile
%                               (set QC = 4 for D profile at qc ~=
%                               FillValue)
%              v3.2  18/06/2019 Marine GALLIAN, Altran Ouest
%                               Supression of  range test, spike test, 
%                               gradient test : quality control have to be
%                               done before using locodox
%              v3.4  02/11/2020 Virginie THIERRY, Ifremer
%                               Add test on QC to move all QC flag = 3 to 1
%                               before computing and applying the
%                               correction
%                               The patch is included while the bug on QC
%                               selection is not found


function [Swork,corrWork] = DOXY_QC(wmo,Sin,Swork,corrWork)


% =========================================================================
%% Initialize the QC array
% =========================================================================
for i=1:size(Sin.doxy_qc.data,1)
    for j=1:size(Sin.doxy_qc.data,2)
        if strcmp(Sin.doxy_qc.data(i,j),'3') == 1
            Sin.doxy_qc.data(i,j)='1';
        end
    end
end

tabdoxyqc = Sin.doxy_qc.data;


% =========================================================================
%% Set of QC tests
% =========================================================================
% Do not take into acount the first descendant profile to compute
% correction
% -------------------------------------------------------------------------
idMoinsUn = find( Sin.cycle_number.data==-1);
if ~isempty(idMoinsUn)
    tabdoxyqc(idMoinsUn,:) = '4';
end

% Do not take into acount descendant profile to compute
% correction
% -------------------------------------------------------------------------
isnok = find(strcmp(cellstr(Sin.direction.data),'D'));
tmp = tabdoxyqc(isnok,:);
isNotBlank = tmp ~= ' ';
tmp(isNotBlank) = '4';
tabdoxyqc(isnok,:) = tmp;

% Remove doxy when no pressure data
% -------------------------------------------------------------------------
tabdoxyqc(~isfinite(Swork.pres_adjusted.data)) = ' ';

%QC variables
% -------------------------------------------------------------------------
corrWork.DOXY_QC = tabdoxyqc;           
Swork.doxy_adjusted_qc = Sin.doxy_qc;
Swork.doxy_qc = Sin.doxy_qc;
Swork.doxy_adjusted_qc.data = tabdoxyqc;             

% Change doxy data in NaN if :
% - doxy_QC not in the limits defined in
% DOXY_config
% - pres unfinite
% -------------------------------------------------------------------------
Swork.doxy_adjusted = Sin.doxy;
Swork.doxy_adjusted. name = [Sin.doxy.name '_ADJUSTED'];
Swork.doxy_adjusted.long_name = [Sin.doxy.long_name ' adjusted and converted'];
Swork.doxy_adjusted.units = corrWork.unit;
Swork.doxy_adjusted.data(~isfinite(Swork.pres_adjusted.data)) = NaN;

%%--AP:
QCp = 1:1:9;       % all possible QCs
QCnan = NaN(1,length(QCp)-length(corrWork.QC_O));
cmpt = 0;
for iq = 1:length(QCp)
    i = find(QCp(iq) == corrWork.QC_O, 1);
    if isempty(i)
        cmpt = cmpt + 1;
        QCnan(cmpt) = QCp(iq);
    end
    clear i
end
for iq = 1:length(QCnan)
    Swork.doxy_adjusted.data(tabdoxyqc==num2str(QCnan(iq))) = NaN;
end
clear tabdoxy;
clear tabdoxyqc;