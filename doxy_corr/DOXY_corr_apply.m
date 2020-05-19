% DOXY_corr_apply applies correction to DOXY data.
%
% SYNTAX
% [CORR] = DOXY_corr_apply(CORR,argoWork,m,b)
%
% DESCRIPTION
% DOXY_corr_apply applies the DOXY (or PSAT) correction and
% propagates it to the multiple profiles from a single cycle.
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
%   argoWork (struct)     Intermediate working argo data structure
%
%   m, b (double)         Coefficient of linear correction (m*x + b),
%                       computed in DOXY_corr_compute.m
% OUTPUT
%   CORR (struct)        Correction structure, with main information and
%                        data from the main profile (choosen vertical
%                        schemes). Completed with the corrected data (DOXY
%                        and PSAT).
%
% CALL : DOXY_convert, O2ctoO2s, O2stoO2c, O2ptoO2c
%
% SEE ALSO
%   DOXY_corr_main, DOXY_corr_compute, DOXY_corr_apply_main

% HISTORY
%   $created: 12/01/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v1.1 25/05/2016   Emilie Brion, ALTRAN OUEST
%                         debugging and improving
%       v2   26/01/2017   Emilie Brion, Altran Ouest
%                         argo - ameliorations 2017 phase 1
%
function [CORR] = DOXY_corr_apply(CORR,argoWork,m,b)


if CORR.presEff == 1, P = argoWork.pres_adjusted.data; else, P = 0; end

switch CORR.whichO2quantity
    case 'DOXY'
        CORR.doxy_corr = m * CORR.(CORR.strField) + b;
        % O2 converted in mL/L for the PSAT computation
        if strcmp(CORR.unit,'mumol/L')
            doxy = CORR.doxy_corr;
        else
            doxy = DOXY_convert(CORR.doxy_corr,CORR.unit,'mumol/L',argoWork.an_dens.data);
        end
        CORR.psat_corr = O2ctoO2s(doxy,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);
        
    case 'PSAT'
        CORR.psat_corr = m * CORR.(CORR.strField) + b;
        CORR.doxy_corr = O2stoO2c(CORR.psat_corr,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);
        CORR.doxy_corr = DOXY_convert(CORR.doxy_corr,'mumol/L',CORR.unit,argoWork.an_dens.data);
        
    case 'PPOX'
        CORR.ppox_corr = m * CORR.(CORR.strField) + b;
        CORR.doxy_corr = O2ptoO2c(CORR.ppox_corr,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);
        CORR.psat_corr = O2ctoO2s(CORR.doxy_corr,argoWork.temp_adjusted.data,argoWork.psal_adjusted.data,P);
        CORR.doxy_corr = DOXY_convert(CORR.doxy_corr,'mumol/L',CORR.unit,argoWork.an_dens.data);
end
