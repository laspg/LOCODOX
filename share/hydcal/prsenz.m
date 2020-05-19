%-------------------------------------------------------------------------------
%
% prsenz        - Conversion pression en profondeur
%
%-------------------------------------------------------------------------------
% Version:
% -------
%  1.01 Création (d'après prsenz, chaine hydro)       	14/06/94 F. Gaillard
%  1.02 test orientation des vecteurs          		24/08/00 F. Gaillard
% 
%-------------------------------------------------------------------------------
%
%
%     description :
%     -----------
%                                                  
%     entree : 
%     ------   
%        p     : tableau des pressions en db
%        sigma : tableau des anomalies de densite in situ (kg/m**3)
%        n     : nombre de points
%        xlat  : latitude du point (degres decimaux)
%  
%     sortie : 
%     ------   
%        z     : tableau de profondeurs (en metres origine en surface
%                axe positif vers le haut)
%
%     sous-programmes appeles : 
%     -----------------------  
%     gravi8
%-------------------------------------------------------------------------------

function [z] = prsenz(p, sigma, xlat)
%

g0    = gravit(xlat);
gamma = -0.2226e-05;
fac   = gamma*0.5/g0;

[np,mp] = size(p);
if mp>np
   p = p';
   sigma = sigma';
end

[n,m] = size(p);


z    = zeros(n,1);
fi   = ones(n,1)./(ones(n,1)+sigma/1000.);
fim  = 0.5*(fi(1:n-1) + fi(2:n));
dp   = p(1:n-1) - p(2:n);

z(1) = -10*p(1)*fi(1);
for i = 2:n
   z(i) = z(i-1) + 10*dp(i-1)*fim(i-1);
end;
z = z./(g0*ones(n,1) + fac*z);
   
