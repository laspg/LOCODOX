% DOXY_write_prepare prepare the data to be writen in monoprofile.
%
% SYNTAX
% [argo1Struct, argo2Struct,argo3Struct,argo4Struct, CONFIG] = DOXY_write_prepare(argo1Struct, argo2Struct, argo3Struct, argo4Struct, CONFIG, Work)
%
% DESCRIPTION
% DOXY_write_prepare prepare the data to be writen in monoprofile : doxy
% unit conversion to original one, fill the DOXY_adjusted_error field,
% replace NaN by FillValue.
%
% INPUT
%   argo1Struct (struct)    Three data structures organised with the main profile
%   argo2Struct (struct)    for working is in argo1Struct (vertical
%   argo3Struct (struct)    sampling shceme of interest for the kind of
%   argo4Struct (struct)    correction), and the others are in the second ones.
%                           The program DOXY_argo_prepare_main.m describes
%                           the choice of the main and other profiles.
%                           Example :
%                             argo1Struct = 
%                                     argo: [1x1 struct]
%                                      Dim: [1x1 struct]
%                                     Work: [1x1 struct]
%                                 argoWork: [1x1 struct]
%                                      VSS: 'Primary sampling'
% 
%                             with: 
%                             argo1Struct.argo =
%                                          pres: [1x1 struct]
%                                     pres_doxy: [1x1 struct]
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
%                                   R2treshold: 0.8000
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
%   CONFIG (struct)       Structure bearing the configuration information.
%                       Created from the configuration file
%
%   Work (structure)  Doxy correction working structure, issued and
%                       computed from argo float data.
%                       Example:
%                             Work = 
%                                   readme: [1x1 struct]
%                                 DOXY_RAW: [85x120 single]
%                                  CAPTEUR: {'Optode'}
%                                  DOXY_QC: [85x120 char]
% OUTPUT
%
%   argo1Struct, argo2Struct, argo3Struct, argo4Struct and CONFIG updated.
%
% CALL : 
%   DOXY_convert, NCR_file, NCW_file
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 19/01/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v2 25/05/2016   Emilie Brion, ALTRAN OUEST
%                       takes into account three Vertical Sampling Scheme
%                       (primary, secondary and nearsurface).

function [argo1Struct, argo2Struct, argo3Struct, argo4Struct, CONFIG] = ...
    DOXY_write_prepare(argo1Struct, argo2Struct, argo3Struct, argo4Struct, CONFIG, Work)

% =========================================================================
%% *Convert DOXY field to the initial unit*
% =========================================================================
for n = 1:4
    argoStruct = eval(['argo' num2str(n) 'Struct;']);    
    if argoStruct.argo.n_prof ~= 0
        %--EB
        % Apply the constant correction if exists
        if isfield(argoStruct.Work,['OFFSET_' Work.whichO2quantity])
            argoStruct.argoWork.doxy_adjusted.data = ...
                DOXY_convert(argoStruct.Work.(['DOXY_OFFSETCORR_' Work.whichO2quantity]),...
                argoStruct.argo.doxy.units,Work.unit,argoStruct.argoWork.an_dens.data);
        else
            argoStruct.argoWork.doxy_adjusted.data = ...
                DOXY_convert(argoStruct.Work.(['DOXY_LINCORR_' Work.whichO2quantity]),...
                argoStruct.argo.doxy.units,Work.unit,argoStruct.argoWork.an_dens.data);
        end
        %--EB
    end
    eval(['argo' num2str(n) 'Struct = argoStruct;']);
end
% =========================================================================
%% *Manage QC :*
%   - Replace QC=2 by 1 and 3 by 4
%   - Find doxy_adjusted with NaN and qc not ' ' or '9' : the QC is put to '4'
%   - a doxy_qc = 4 should be associated with a doxy_adjusted_qc = 4 and a
%   doxy_adjusted at FillValue
% =========================================================================
for n = 1:4
    argoStruct = eval(['argo' num2str(n) 'Struct']);
    if argoStruct.argo.n_prof ~= 0                
        % QC 2 are replaced by 1
        is2 = argoStruct.argoWork.doxy_qc.data == '2';
        argoStruct.argoWork.doxy_qc.data(is2) = '1';
        
        % QC 3 are replaced by 4
        is3 = argoStruct.argoWork.doxy_qc.data == '3';
        argoStruct.argoWork.doxy_qc.data(is3) = '1';
                
        % QC = 4 and doxy_adjusted to FillValue
        is4 = argoStruct.argoWork.doxy_qc.data == '4';
%         argoStruct.argoWork.doxy_adjusted_qc.data(is4) = NaN;
        argoStruct.argoWork.doxy_adjusted.data(is4) = NaN;
        
        % Doxy_adjusted fillValue associated with a doxy_QC different from
        % 9 or blank is put to 4
        isFill = isnan(argoStruct.argoWork.doxy_adjusted.data);
        isNotBlankOr9 = argoStruct.argoWork.doxy_qc.data ~= ' ' & ...
                        argoStruct.argoWork.doxy_qc.data ~= '9';
        isToChange = isFill & isNotBlankOr9;
        argoStruct.argoWork.doxy_qc.data(isToChange) = '4';                
        
        eval(['argo' num2str(n) 'Struct = argoStruct;']);
    end
end
% =========================================================================
%% *Put the adjusted error values*
% =========================================================================
% modif C.Lagadec 12/04/17
% CONFIG.adjusted_error_abs et CONFIG.adjusted_error_rel sont des valeurs  numériques
% j'ai donc supprimé str2double dans la multiplication
% sinon l'on obtient NaN dans doxy_adjusted_error
% sans doute pas bon dans le cas de CHANGE_VALUES, mais rien modifié
%--EB
% CONFIG.adjusted_error_abs, _rel are already double
% (feval(DOXY_config.m) => str2double/str2num is no longer required
%--EB
% =========================================================================

cv = 1;
while cv
    
    whichError = questdlg({'Two kind of error : ABSOLUTE or RELATIVE';...
        sprintf('In the configuration : ABSOLUTE = %d, RELATIVE = %d%%',...
        CONFIG.adjusted_error_abs, CONFIG.adjusted_error_rel);...
        'Which kind of error would you like to apply in doxy_adjusted_error ?'},...
        'ERROR ON DOXY ADJUSTED','ABSOLUTE','RELATIVE','CHANGE VALUES','RELATIVE'); 
   
    switch whichError
        case 'ABSOLUTE'
            cv = 0;
            for n=1:4
                argoStruct = eval(['argo' num2str(n) 'Struct']);
                if argoStruct.argo.n_prof ~= 0
                    argoStruct.argoWork.doxy_adjusted_error.data = (argoStruct.argoWork.doxy_adjusted.data).*0 + CONFIG.adjusted_error_abs;
                end
                eval(['argo' num2str(n) 'Struct = argoStruct;']);
            end
            
            % --------   SUMMARY FILE   --------------
            fprintf(Work.log_summ,'\n');
            fprintf(Work.log_summ,sprintf('The error applied on doxy adjusted field is an absolute error equal to %s umol/kg (DOXY_ADJUSTED_ERROR = %s umol/kg)\n',num2str(CONFIG.adjusted_error_abs),num2str(CONFIG.adjusted_error_abs)));
            fprintf(Work.log_summ,'\n');
            fprintf(Work.log_summ,'----------------------------------------------------\n');
            
        case 'RELATIVE'
            cv = 0;
            for n=1:4
                argoStruct = eval(['argo' num2str(n) 'Struct']);
                if argoStruct.argo.n_prof ~= 0
                    argoStruct.argoWork.doxy_adjusted_error.data = argoStruct.argoWork.doxy_adjusted.data.*CONFIG.adjusted_error_rel/100;
                end
                eval(['argo' num2str(n) 'Struct = argoStruct;']);
            end
            
            % --------   SUMMARY FILE   --------------
            fprintf(Work.log_summ,'\n');
            fprintf(Work.log_summ,sprintf('The error applied on doxy adjusted field is a relative error equal to %s percent (DOXY_ADJUSTED_ERROR = DOXY_ADJUSTED * %s/100)\n',num2str(CONFIG.adjusted_error_rel),num2str(CONFIG.adjusted_error_rel)));
            fprintf(Work.log_summ,'\n');
            fprintf(Work.log_summ,'----------------------------------------------------\n');
                       
        case 'CHANGE VALUES'
            answer = inputdlg({'change the raw error:';...
                               'change the relative error:'},...
                               'Change adjusted error',1,...
                               {num2str(CONFIG.adjusted_error_abs),num2str(CONFIG.adjusted_error_rel)});
            CONFIG.adjusted_error_abs = str2num(answer{1}); %#ok<ST2NM>
            CONFIG.adjusted_error_rel = str2num(answer{2}); %#ok<ST2NM>
    end
end


% =========================================================================
%% * Put NaN to FillValue*
% =========================================================================
argo1Struct.argo = nan_2_fillValue(argo1Struct.argo);
argo2Struct.argo = nan_2_fillValue(argo2Struct.argo);
argo3Struct.argo = nan_2_fillValue(argo3Struct.argo);
