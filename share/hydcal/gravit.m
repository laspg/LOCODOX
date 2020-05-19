%-------------------------------------------------------------------------------
%
% gravit        - Calcul de l'acceleration de gravite en fonction
%                 de la latitude en degres decimaux. (formule GRS-80)
%
%-------------------------------------------------------------------------------
% Version:
% -------
%  1.01 Création (d'après gravi8, chaine hydro)          14/06/94 F. Gaillard
% 
%-------------------------------------------------------------------------------


function [g] = gravit(xlat);
degrad = pi/180.0;
alat   = xlat*degrad;
slat   = sin(alat);
slat2  = sin(alat*2.00);
ge     = 9.780318;
b2     = 0.530244d-02;
b4     = -0.585d-05;
g = ge*(1.00 +b2*slat*slat+b4*slat2*slat2);
