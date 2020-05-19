%-------------------------------------------------------------------------------
%
% fbruva        - Calcul de la frequence de Brunt-Vaisala
%
%-------------------------------------------------------------------------------
%
%  appel:   function [Nsq, N, Nsf] = fbruva(Sal, Temp, Pres,  xlat, idec)
%  -----  
%
% Remarque: le tableau doit etre echantillonne regulierement en Pression
% --------  fonctionne pour des vecteurs, pas des matrices
% 	
%  1.01 Création (d'après fbruva.f, chaine hydro) 	14/06/94 F. Gaillard
%  1.02 ajout du decalage  idec    			15/01/99 F. Gaillard
%       - si idec est inferieur a 1, il est mis automatiquement a 1
%       - pour les calculs - idec optionnel
%  1.03 integration dans la chaine           		20/01/99 Cathy
%  1.04 suppression de la boucle pour appel tetai       13/09/99 (C. Lagadec)  
%-------------------------------------------------------------------------------

function [Nsq, fbv, Nsf] = fbruva(Sal, Temp, Pres,  xlat, idec)

% idec est optionnel => si le nombre d'arguments est < 5,
%                       on suppose que idec est manquant,
%                       on le force a 1
% idec doit etre au moins egal a 1

if  nargin < 5 
     idec = 1;
else
 if  idec < 1 
     idec = 1;
 end;
end;

%  Calcule les constantes geophysiques locales
%  -------------------------------------------
g0 = gravit(xlat);
gdeux = g0*g0*1.0d-4;
fcor  = pi*sin(pi*xlat/180.0)/(6*3600);

nz = length(Pres);
i0 = 1 + idec;
in = nz - idec;

%  Calcule la densite in situ:
%  --------------------------
[alpha, Sig] = swstat90 (Sal, Temp, Pres);


Psup = NaN*ones(nz,1);
Pinf = NaN*ones(nz,1);
Tsup = NaN*ones(nz,1);
Tinf = NaN*ones(nz,1);
Sigsup = NaN*ones(nz,1);
Siginf = NaN*ones(nz,1);


%  Deplace la particule de idec niveau vers le haut et vers le bas:
%  ---------------------------------------------------------------

Psup(i0:in) = Pres(i0 - idec:in - idec);
Pinf(i0:in) = Pres(i0 + idec:in + idec);

% -----------------------------------------------------------------------
% modif 13/09/99 (C. Lagadec) : suppression de la boucle pour appel tetai
% -----------------------------------------------------------------------

%for iz = i0:in
%   Tsup(iz) =  tetai(Pres(iz), Temp(iz), Sal(iz), Psup(iz));
%   Tinf(iz) =  tetai(Pres(iz), Temp(iz), Sal(iz), Pinf(iz));
%end

 Tsup =  tetai(Pres, Temp, Sal, Psup);
 Tinf =  tetai(Pres, Temp, Sal, Pinf);


[alpha, Sigsup(i0:in)] = swstat90 (Sal(i0:in), Tsup(i0:in), Psup(i0:in));
[alpha, Siginf(i0:in)] = swstat90 (Sal(i0:in), Tinf(i0:in), Pinf(i0:in));

%  Calcule N**2:
%  ------------
%  Nsq = -(g/rho)* drho/dz = 10**-4 * g*g * dsig/dp

Nsq = NaN*ones(nz,1);
fbv = NaN*ones(nz,1);
%dsup = (Sigsup(i0+1:in-1) - Sig(i0:in-2))./(Psup(i0+1:in-1) - Pres(i0+1:in-1));
%dinf = (Siginf(i0+1:in-1) - Sig(i0+2:in))./(Pinf(i0+1:in-1) - Pres(i0+1:in-1));
%Nsq(i0+1:in-1) = gdeux*(dsup+dinf)/2;

dsup = (Sigsup(i0:in) - Sig(i0-idec:in-idec))./(Pres(i0:in) - Psup(i0:in));
dinf = (Siginf(i0:in) - Sig(i0+idec:in+idec))./(Pres(i0:in) - Pinf(i0:in));
Nsq(i0:in) = gdeux*(dsup+dinf)/2;


%  Calcule Fbv et N/f:
%  ------------------

isnotok = find(Nsq<0);
bid = Nsq;
bid(isnotok) = -bid(isnotok);
fbv = sqrt(bid);
fbv(isnotok) = -fbv(isnotok);
Nsf = fbv/fcor;

clear bid Psup Pinf Tsup Tinf Sigsup Siginf alpha;

