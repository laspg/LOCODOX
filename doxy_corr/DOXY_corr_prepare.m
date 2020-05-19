% DOXY_corr_prepare prepare data for correction applying.
%
% SYNTAX
% [CORR2, Work2, argo2Work] = DOXY_corr_prepare(CORR, Work, argoStruct, argoTrajWork)
%
% DESCRIPTION
% DOXY_corr_prepare prepare data for correction to be applied.The
% correction has be computed only with ascending profile, but has to be
% applied both over ascending and descending profiles. The programm gets
% these data, prepares it and compute drift from the drift computed over
% ascending profiles only.
%
% INPUT
%   CORR (struct)        Correction structure, with main information and
%                        data from the main profile (ascending only,
%                        vertical schemes used for the correction).
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
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA),
%                        completed with the annual drift.
%                        Example:
%                        Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                                     DRIFT: 4.8400
%
%   argoStruct (struct)  float structure from the main profile (choosen
%                        vertical schemes for the correction)
%                        Example :
%                         argo1Struct =
%
%                                 argo: [1x1 struct]
%                             argoWork: [1x1 struct]
%                                 Work: [1x1 struct]
%                         with:
%                         argo1Struct.argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                                         ...
%
%                         argo1Struct.argoWork =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%                                       sat: [1x1 struct]
%                                      psat: [1x1 struct]
%
%                          argo1Struct.Work =
%                                    readme: [1x1 struct]
%                                      unit: 'mumol/kg'
%                                       wmo: 1901205
%                                 whichCorr: 'WOA'
%                                  DOXY_RAW: [165x117 single]
%                                     timar: [165x20 char]
%                                     datat: [165x1 double]
%                                   DENSITY: [165x117 single]
%                                     DEPTH: [165x117 double]
%                                   CAPTEUR: 'Optode'
%                                   DOXY_QC: [165x117 char]
%                                      DENS: [165x117 single]
%                                           ...
%  argoTrajWork (structure)     float data structure (argoTrajWork) choose 
%                               for the Doxy correction.
%                           Example: argoTrajWork
%                          float_serial_no: [1x1 struct]
%                            wmo_inst_type: [1x1 struct]
%                                     juld: [1x1 struct]
%                                 latitude: [1x1 struct]
%                                longitude: [1x1 struct]
%
% OUTPUT
%   CORR2 (struct)       Correction structure, with main information and
%                        data from the secondary profile.
%
% CALL :
%O2ctoO2p, DOXY_convert, O2ctoO2s, O2ptoO2c
% SEE ALSO
%   DOXY_compute_corr, DOXY_apply_corr, DOXY_corr_main

% HISTORY
%   $created: 14/01/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%              v2 26/01/2017   Emilie Brion, Altran Ouest
%                              argo - ameliorations 2017 phase 1
%              v2.1 09/11/2017 Emilie Brion, Altran Ouest
%                              drift correction is now re-computed using
%                              polyval or feval, instead of a default case
%                              ycorrected = raw - a*daydiff => becomming
%                              ycorrected = raw - fitRegression and
%                              fitRegression = polyval(driftCoef,daydiff)
%                              or feval(function_fit,daydiff).
%                 26/10/2018   Marine GALLIAN, Altran Ouest :
%                              Add the possibility to apply a 
%                              constant drift from a certain day
%                 25/03/2019   Marine GALLIAN, Altran Ouest 
%                              Fix unit probleme for applying drift : 
%                              Convert PPOX data when drift is computed on
%                              WOA 
%                05/02/2020    Add conversion from mumol/L to mumol/kg when 
%                              drift is computed on PPOX. V. Thierry

function [CORR2] = DOXY_corr_prepare(CORR, Work, argoStruct, argoTrajWork)

% =========================================================================
%% Initialisation
% =========================================================================
% arguments
minArgs = 3;
maxArgs = 4;
narginchk(minArgs,maxArgs);

% Initialize
Work2 = argoStruct.Work;
argo2Work = argoStruct.argoWork;
argo2 = argoStruct.argo;
Work2.savePlot = Work.savePlot;

CORR2.presOrDens = CORR.presOrDens;
CORR2.cmpl = CORR.cmpl;
CORR2.whichCorr = CORR.whichCorr;
CORR2.whichO2quantity = CORR.whichO2quantity;
if isfield(Work,'refCyc') || strcmp(Work.whichCorr,'REF')
    CORR2.refCyc = CORR.refCyc;
    CORR2.nbprof = CORR.nbprof;
else
    okdim = strcmp(argo2Work.pres_adjusted.dim,'N_PROF');
    dim = size(argo2Work.pres_adjusted.data);
    nbprof = dim(okdim);
    CORR2.refCyc = 1:nbprof;    
    CORR2.nbprof = nbprof; %--MG
end

CORR2.unit = CORR.unit;
CORR2.presEff = CORR.presEff;
CORR2.strField = CORR.strField;
field = CORR2.strField;
if CORR.presEff == 1, P = argo2Work.pres_adjusted.data; else, P = 0; end
argo2Work.raw.doxy.data = argo2.doxy.data;
argo2Work.raw.doxy.units = argo2.doxy.units;
tmpDoxy = DOXY_convert(argo2Work.raw.doxy.data,argo2Work.raw.doxy.units,'mumol/L',argo2Work.an_dens.data);
argo2Work.raw.ppox.data = O2ctoO2p(tmpDoxy,argo2Work.temp_adjusted.data,argo2Work.psal_adjusted.data,P);
argo2Work.raw.psat.data=O2ctoO2s(tmpDoxy,argo2Work.temp_adjusted.data,argo2Work.psal_adjusted.data,P); %marine 12/06/19

% pressure effect taking into account or not
if Work.presEff == 1, P = argo2Work.pres_adjusted.data; else, P = 0; end

% Nb of day after deployment
dayjul = argo2.juld.data - argo2.juld.data(1);

% =========================================================================
%% Compute drift if needed
% The value comes from the Work struct. For Doxy correction, the drift has
% been computed from ascending profiles only.
% =========================================================================
if Work.DODRIFT==1
    
    % --------------------------------------------------------------------- 
    % Prepare
    % ---------------------------------------------------------------------  
    if isfield(Work,'isargo1Struct')
        warning('The drift computed on the ascending profiles of cycle N is applied on both descending and ascending profiles of cycle N. The right thing to do is to applied the computed drift to ascending profile of cycle N and to descending profile of cycle N+1');
    end
    % initialize
    okdim = strcmp(argo2Work.pres_adjusted.dim,'N_PROF');
    dim = size(argo2Work.pres_adjusted.data);
    CORR2.(field) = NaN(size(argo2Work.pres_adjusted.data));
    
    % the drift regression coefficient have been computed using DOXY.
    % For correction with PSAT, convert DOXY corrected to PSAT
    if ~strcmp(Work.whichDrift,'NCEP')
        O2DriftCoef = 'DOXY';
    else
        O2DriftCoef = 'PPOX';        
    end
    if Work.drift_spec == 0
        fitRegression = polyval(Work.([O2DriftCoef '_DRIFTCORR_COEF']),double(dayjul));
        if iscolumn(fitRegression)
            fitRegression = fitRegression';
        end
    else
        listCoef = [];
        for i=1:length(Work.([O2DriftCoef '_DRIFTCORR_COEF']))
            listCoef = [listCoef,sprintf('Work.([O2DriftCoef ''_DRIFTCORR_COEF''])(%d),',i)];
        end
        listCoef(end) = '';
        eval(['c = cfit(Work.drift_fittype,' listCoef ');']);
        fitRegression = feval(c,double(dayjul));
    end
    
    %Compute constant drift
    if isfield(Work,'ind_drift_stop')
            for i=1:length(dayjul) %find indice to stop drift 
                if dayjul(i)>Work.ind_drift_stop(2) && Work.ind_drift_stop(2)>-1
                    Work.ind_drift_stop(1)=i;
                    break
                end
            end
        fitRegression(Work.ind_drift_stop(1):end)=fitRegression(Work.ind_drift_stop(1));
    end
    CORR2.drift=fitRegression;
    
    % --------------------------------------------------------------------- 
    % Apply the drift 
    % ---------------------------------------------------------------------     
    if strcmp(Work.whichDrift,'WOA') %drift in DOXY
        for i = 1:dim(okdim)
            CORR2.doxy(i,:) = argo2Work.raw.doxy.data(i,:) - fitRegression(i);
        end
        tmpDoxy = DOXY_convert(CORR2.doxy,CORR2.unit,'mumol/L',argo2Work.an_dens.data);
        CORR2.psat = O2ctoO2s(tmpDoxy,argo2Work.temp_adjusted.data,argo2Work.psal_adjusted.data,P);        
        CORR2.ppox = O2ctoO2p(tmpDoxy,argo2Work.temp_adjusted.data,argo2Work.psal_adjusted.data,P);  
        
    elseif strcmp(Work.whichDrift,'NCEP') %drift in PPOX
        for i = 1:dim(okdim)
            CORR2.ppox(i,:) = argo2Work.raw.ppox.data(i,:) - fitRegression(i);
        end
        %VT 05022020
        % correction d'un probleme de conversion d'unite en mumol/L et
        % mumol/kg
       tmpDoxy = O2ptoO2c(CORR2.ppox,argo2Work.temp_adjusted.data,argo2Work.psal_adjusted.data,P);
       CORR2.doxy = DOXY_convert(tmpDoxy,'mumol/L',CORR2.unit,argo2Work.an_dens.data);
       CORR2.psat = O2ctoO2s(tmpDoxy,argo2Work.temp_adjusted.data,argo2Work.psal_adjusted.data,P);
       % ancien code
%        CORR2.doxy = O2ptoO2c(CORR2.ppox,argo2Work.temp_adjusted.data,argo2Work.psal_adjusted.data,P);
%        tmpDoxy = DOXY_convert(CORR2.doxy,CORR2.unit,'mumol/L',argo2Work.an_dens.data);
%        CORR2.psat = O2ctoO2s(tmpDoxy,argo2Work.temp_adjusted.data,argo2Work.psal_adjusted.data,P);
        %VT
    end
    
else
    CORR2.(field) = argo2Work.raw.(field).data;
    CORR2.drift=[];
end
