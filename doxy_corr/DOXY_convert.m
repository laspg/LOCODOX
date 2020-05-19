%  DOXY_convert converts the oxygene data from a unit to another.
%
% SYNTAX
% [valout]=DOXY_convert(valin,uin,uout,pdensAno)
%
% DESCRIPTION
% DOXY_convert converts the oxygene data from a unit to another
%
% INPUT
%   valin (double)     oxygen data (in uin unit) to be converted
%
%   uin (string)       orginal unit of valin
%
%   uout (string)      new unit to convert valin
%
%   pdensAno (double)  (OPTIONNAL) potential density reference anomaly
%                      anomalie de densité potentielle référence à 0 de
%                      l'eau considérée, si aucun valeur n'est entrée, on
%                      considère qu'on a de l'eau douce et l'anomalie de
%                      densité est nulle
% OUTPUT
%   valout (double)    oxygen data (in uin unit) to be converted
%
% REMARK
%   Valid units :
%         mL/L: milliliter per liter
%         mmol/m3: millimole per m3
%         mumol/L: micromole per liter
%         mg/L : milligram per liter
%         mumol/kg : micromole per kilo
%         mL/L * 44.6596 = mmol/m3 = mumol/L
%         mL/L*1.42903 = mg/L
%         mg/L * 44.66/1.42903 = mumol/L = mmol/m3
%         mumol/L = (rho/1000)* mumol/kg
%
% OUTPUT
%   valout (?)          oxygen data converted in uout unit.
%
% CALL 
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: //2009 $author: ?, LPO, CNRS
%   $Revision: version $Date: $author:
%       v2 21/11/2015   Emilie Brion, ALTRAN OUEST
%                       adapted to be shared to the O2 community

function [valout] = DOXY_convert(valin,uin,uout,pdensAno)

% =========================================================================
%% Initialisation
% =========================================================================

if nargin == 3
    rho=1000;
elseif nargin == 4
    rho=pdensAno+1000;
end

coef1=44.6596;
coef2=1.42903;
coef3=1000./rho;

valout=[];

% Checks unit in and out in possible name, and define the different cases
% -------------------------------------------------------------------------
% Possible name for doxy units
doxyUnits.molPerVol = {'mumol/L'  'mumol/m3' 'mmol/L'  'mmol/m3' 'micromole per liter' 'micromole per m3' 'micromole/L'  'micromole/m3'};
doxyUnits.molPerWeight = {'mumol/kg' 'mmol/kg' 'micromole per kilogram' 'micromole/kg'};
doxyUnits.volPerVol = {'mL/L' 'milliliter per liter' 'milliliter/L'};
doxyUnits.weightPerVol = {'mg/L' 'milligram per liter' 'milligram/L'};
doxyUnitFields = fieldnames(doxyUnits);

% Find uin and uout in doxyUnits
for f = 1:length(doxyUnitFields)
    if any(strcmpi(uin,doxyUnits.(doxyUnitFields{f})))
        ucaseIn = doxyUnitFields{f};
    end
    if any(strcmpi(uout,doxyUnits.(doxyUnitFields{f})))
        ucaseOut = doxyUnitFields{f};
    end
end

% =========================================================================
%% Unit conversion
% =========================================================================

switch ucaseIn
    case 'molPerVol'
        switch ucaseOut
            case 'molPerVol'
                valout = valin;                
            case 'molPerWeight'
                valout=valin.*coef3;
            case 'volPerVol'
                valout=valin/coef1;
            case 'weightPerVol'
                valout=valin*coef2/coef1;
        end
   case 'molPerWeight'
        switch ucaseOut
            case 'molPerVol'
                valout=valin./coef3;
            case 'molPerWeight'
                valout = valin;                
            case 'volPerVol'
                valout=valin./(coef1.*coef3);
            case 'weightPerVol'
                valout=valin*coef2./(coef3*coef1);
        end
   case 'volPerVol'
        switch ucaseOut
            case 'molPerVol'
                valout=valin*coef1;
            case 'molPerWeight'
                valout=valin*coef1.*coef3;
            case 'volPerVol'
                valout = valin;                
            case 'weightPerVol'
                valout=valin*coef2;
        end
   case 'weightPerVol'
        switch ucaseOut
            case 'molPerVol'
                valout=valin*coef1/coef2;
            case 'molPerWeight'
                valout=valin.*coef3*coef1/coef2;
            case 'volPerVol'
                valout=valin/coef2;
            case 'weightPerVol'
                valout = valin;                
        end
end

