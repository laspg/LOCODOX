%-------------------------------------------------------------------------------
%
% zenprs        - Conversion profondeur en pression 
%                 appel: function [p] = zenprs(z, sigma, xlat)
%
%-------------------------------------------------------------------------------
%
%     description :
%     -----------
%                                                  
%     entree : 
%     ------   
%        z     : tableau de profondeurs (en metres origine en surface
%                axe positif vers le haut)
%        sigma : tableau des anomalies de densite in situ (kg/m**3)
%        n     : nombre de points
%        xlat  : latitude du point (degres decimaux)
%  
%     sortie : 
%     ------   
%        p     : tableau des pressions en db
%
%     sous-programmes appeles : 
%     -----------------------  
%     gravi8
%
% Version:
% -------
%  1.01 Création (d'après zenprs, chaine hydro)          14/03/95 F. Gaillard
% 
%-------------------------------------------------------------------------------

function [p] = zenprs(z, sigma, xlat)


g0    = gravit(xlat);
gamma = -0.2226e-05;

n = length(z);
p = zeros(n,1);

fi = (1000.+sigma).*( g0+gamma*z);
fi = (fi(1:n-1) + fi(2:n))*0.5e-04;
dz = z(1:n-1) - z(2:n);
dzfi = dz.*fi;

%      
p(1) = -z(1)*(1000.+sigma(1))*(g0+0.5*gamma*z(1))*1.0e-04;

for i = 2:n,
   p(i) = p(i-1) + dzfi(i-1);
end;

clear dz fi dzfi;
