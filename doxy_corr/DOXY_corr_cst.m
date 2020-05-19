% DOXY_corr_cst computes the constant correction.
%
% SYNTAX
% [Work] = DOXY_corr_cst(CORR,Work,argoWork,dec)
%
% DESCRIPTION
% DOXY_corr_cst computes the constant correction to the PSAT or DOXY data.
% The offset factor is computed in DOXY_corr_apply_main.m. 
%
% INPUT
%   CORR (struct)        Correction structure, with main information and
%                        data from the main profile (choosen vertical
%                        schemes).
%                        Example :
%                             CORR = 
%                          presOrDens: 'Pressure'
%                               level: [164x117 single]
%                                cmpl: '_interpP'
%                           whichCorr: 'WOA'
%                              refCyc: [1x164 double]
%                     whichO2quantity: 'PSAT'
%                            strField: 'psat'
%                                psat: [164x117 single]
%                             psatref: [164x117 double]
%                             psat_ok: [164x117 double]
%                                   ...
%
%   Work (struct)        Working argo data structure
%                        Example:
%                         Work = 
% 
%                                         readme: [1x1 struct]
%                                     R2treshold: '0.80'
%                             adjusted_error_abs: '0'
%                             adjusted_error_rel: '5'
%                                           unit: 'mumol/kg'
%                                        dirPlot: '/home1/homedir4/perso/ebrion/NAOS/NAOS_2016(1)/plots2/REF'
%                                       makePlot: 1
%                                       savePlot: 1
%                                       fontsize: 8
%                                        history: [1x1 struct]
%                                            wmo: 5902302
%                                         sensor: 'Optode'
%                                      whichCorr: 'REF'
%                                       DOXY_RAW: [175x110 single]
%                                          timar: [175x1 datetime]
%                                          datat: [175x1 double]
%                                        DENSITY: [175x110 single]
%                                          DEPTH: [175x110 double]
%                                        DOXY_QC: [175x110 char]
%                                           DENS: [175x110 single]
%                                           PSAT: [175x110 single]
%                                           DOXY: [175x110 single]
%                                           PSAL: [175x110 single]
%                                           TEMP: [175x110 single]
%                                           PRES: [175x110 single]
%                                     PresOrDens: 'Pressure'
%                               DOXY_WOA_interpP: [175x110 double]
%                               PSAT_WOA_interpP: [175x110 double]
%                                               ...
%
%   argoWork (struct)    Intermediate working argo data structure
%
%   dec (double)         Offset to apply. Computed as the mean of the
%                        difference between argo and the reference data
%
% OUTPUT
%   Work (struct)        Working argo data structure completed
%
% CALL : 
%   DOXY_convert, O2stoO2c, O2ctoO2s
%
% SEE ALSO
%   DOXY_corr_main, DOXY_corr_compute, DOXY_corr_apply_main
%
% HISTORY
%   $created: 25/05/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%             02/04/2019    Gallian Marine, Altran Ouest
%                           Error of correction : the offset wasn't apply
%                           to PSAT for the case "DOXY" and to DOXY for the
%                           case "PSAT"

function [Work] = DOXY_corr_cst(CORR,Work,argoWork,dec)

% -------------------------------------------------------------------------
% Offset correction
% -------------------------------------------------------------------------
if CORR.presEff == 1, P = argoWork.pres_adjusted.data; else, P = 0; end

switch CORR.whichO2quantity
    case 'DOXY'
        doxy_cte = CORR.doxy + dec;
        if ~isfield(CORR,'no_log')
            fprintf('\t offset of %2.3f \n', dec);
        end
        Work.(['OFFSET_' CORR.whichO2quantity]) = dec;
        Work.(['DOXY_OFFSETCORR_' CORR.whichO2quantity]) = doxy_cte;                
%         % Transformation en ml/l pour le calcul de %SAT corrige
%         doxy_ml = DOXY_convert(doxy_cte,CORR.unit,'mL/L',argoWork.an_dens.data);
%         %   Saturation et saturation percentage
%         Work.(['PSAT_OFFSETCORR_' CORR.whichO2quantity]) = (doxy_ml./argoWork.sat.data)*100;
        % Transformation en mumol/l pour le calcul de %SAT corrige
        doxy = DOXY_convert(doxy_cte,CORR.unit,'mumol/L',argoWork.an_dens.data);
        %   Saturation et saturation percentage
        Work.(['PSAT_OFFSETCORR_' CORR.whichO2quantity]) = O2ctoO2s(doxy,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);
       
        
    case 'PSAT'
        if ~isfield(CORR,'no_log')
            fprintf('\t offset of %2.3f \n', dec);
        end
        
        Psat_cte = CORR.psat + dec;
        
%         %Transformation en ml/l pour le calcul de %SAT corrige
%         doxy_ml = (Psat_cte.* argoWork.sat.data)/100;
%         if strcmp(CORR.unit,'mL/L')
%             doxy_cte = doxy_ml;
%         else
%             doxy_cte = DOXY_convert(doxy_ml,'mL/L',CORR.unit,argoWork.an_dens.data);
%         end
%         Work.(['DOXY_OFFSETCORR_' CORR.whichO2quantity]) = doxy_cte;
%         Work.(['PSAT_OFFSETCORR_' CORR.whichO2quantity]) = Psat_cte;               
        doxy_cte = O2stoO2c(Psat_cte,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);
        doxy_cte = DOXY_convert(doxy_cte,'mumol/L',CORR.unit,argoWork.an_dens.data);
        Work.(['OFFSET_' CORR.whichO2quantity]) = dec;
        Work.(['DOXY_OFFSETCORR_' CORR.whichO2quantity]) = doxy_cte;
        Work.(['PSAT_OFFSETCORR_' CORR.whichO2quantity]) = Psat_cte;               
end
