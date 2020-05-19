%------------------------------------------------------------------------------
%
% dynam0        - calculates vertical modes of an arbitrary N**2 profile
%                 with free surface 
%
%                !!! Attention: Utiliser des N**2 non nuls pour le moment !!!
%                !!! Attention: verifier condition surface libre !!!
%
%  call: function [rossby, fmod, gmod] = dynam(bv2, zz, ylat, nmod, iop)
%
%------------------------------------------------------------------------------
%
% Version:
% -------
%  1.01 Création (d'après dynam, chaine hydro)          21/08/96  F. Gaillard
% 
%
%  Method :
%  ------
%  Applying the change of variable: g = fz/bv2
%  converts the problem to a matrix eigenvalue problem:
%                                  gzz/bv2    = ev*g
%  subject to the boundary conditions:  g(-d) = 0
%                          a)     -gz(0)/grav = ev*g(0)
%                          b)            g(0) = 0
%   where   p(x,y,z,t) = a(x,y,t)*f(z) et f(z) = gz(z)/ev
%    f   is the stream function vertical mode
%    g   is the mode for vertical displacement
%    ev  is the eigenvalue
%    bv2 is the buoyancy frequency squared
%   grav is the gravitational accel  
%
%
%** approximates second partial with respect to z by drawing a
%   parabola through three points and taking its curvature as the
%   second partial at the center point.
%
%** approximates the first partial with respect to z (which appears
%   in the boundary condition) as a simple finite difference:
%   gz = (g(2)-g(1)) / (z(2)-z(1))
%   
%** Dans les tranches ou N**2 est nulle, le mode sera lineaire de i-1 a j
%       g(k) = a(k)*g(j) + b(k)*g(i-1)
%       where       a(k) = (z(k)-z(i-1))/(z(j)-z(i-1))
%                   b(k) = (z(j)-z(k))/(z(j)-z(i-1))
%   on peut eliminer g(k), i<= k <= j-1
%   donc
%           -l(i-1)*g(i-2) + d(i-1)*g(i-1) - u(i)*g(i) = ev*g(i-1)
%   la derniere equation avant le premier niveau n = 0 devient:
%     -l(i-1)*g(i-2)+[d(i-1)u(i)*b(i)]*g(i-1)-u(i)*a(i)*g(i)=ev*g(i-1)
%   et de meme pour la j-eme equation
%
%  donnees en entree:
%  -----------------
%     bv2  : Brunt-Vaisala frequency squared
%     zz   : vertical profile, first point is assumed to be z = 0 (surface)
%    ylat  : latitude en degres decimaux
%    nmod  : nombre de modes a calculer
%     iop  : si = 1 impose la condition de surface libre
%            sinon g(0) = 0
%
%  sorties:
%  -------
%  rossby   : rayon de deformation de rossby des modes ,
%             calcule a partir de la valeur propre ev:
%             rossby(i) = 1/(f0*sqrt(ev(i)))
%    fmod   : modes pour la fonction courant
%    gmod   : modes normalises pour le deplacement vertical
%
%------------------------------------------------------------------------------
                                              
function [rossby, fmod, gmod] = dynam(bv2, zz, ylat, nmod, iop)

%
%  Constantes:
bvmin   = 1.0e-12;

grav    = gravit(ylat);
fcor    = pi*sin(pi*ylat/180.0)/(6*3600);

%  S'assure que les tableaux sont verticaux:
[n,m] = size(bv2);
if m>n,
   bv2 = bv2';
   zz  = zz';
   nz = m;
else
   nz = n;
end;

%  repere les niveaux ou N**2 est nulle:
%  ------------------------------------
bvnul = find(bv2<bvmin);
nzeros = length(bvnul);
if nzeros>0
   fprint (1, '  Valeurs nulles de N**2 rencontrees attention \n');
   itrou = 1;
   i     = 1;
   while i<nzeros+1
      ideb(itrou) = bvnul(i) - 1;
      if i == nzeros
         ifin(itrou) = bvnul(i) + 1;
      else
         if bvnul(i+1)==bvnul(i) + 1
            i = i + 1;
         else
            ifin(itrou) = bvnul(i) + 1;
            itrou = itrou + 1;
         end
      end
   end
else
   ntrou = 0;
end

%  Cas general flag mis a 1:
%  ------------------------
iflag        = ones(nz,1);
if nzeros>0,
   iflag(bvnul) = zeros(size(bvnul));
end
iflag(1)     = 0;
iflag(nz)    = 0;
isok         = find(iflag);

%  ==========================================
%     Rempli les matrices tridiagonales:
%  ==========================================
%    tridiag(i,1)   = diagonale inferieure
%    tridiag(i,2)   = diagonale principale
%    tridiag(i,3)   = diagonale superieure

tridiag = nan*ones(nz,3);
dz1     = zz(2:nz) - zz(1:nz-1);
dz2     = [0 ; zz(3:nz) - zz(1:nz-2)];
bvinv   = -ones(size(isok))./bv2(isok);

%  calcule les diagonales quand N**2 est non nulle (et hors des bornes)
%  -------------------------------------------------------------------
%   approximates second partial with respect to z by drawing a
%   parabola through three points and taking its curvature as the
%   second partial at the center point.
tab1 =  dz1(isok  ).*dz2(isok);
tab2 = -dz1(isok-1).*dz1(isok);
tab3 =  dz1(isok-1).*dz2(isok);
tridiag(isok,1) = 2*bvinv./tab1;
tridiag(isok,2) = 2*bvinv./tab2;
tridiag(isok,3) = 2*bvinv./tab3;


%  !!!!!!!!!!!!!!!!!!!!   Partie a verifier - debut !!!!!!!!!!!!!!!!!!!!!!!!!!!

%  calcule les diagonales aux limites de N*2 nulle   
%  -----------------------------------------------
%  transforme l'equation ideb:
if nzeros>0,
   deltaz = zz(ifin) - zz(ideb);
   tridiag(ideb,1) = tridiag(ideb,1) ...
                  - tridiag(ideb,3).*(zz(ifin) - zz(ideb+1))./deltaz;
   tridiag(ideb,3) = tridiag(ideb,3).*(zz(ideb+1) - zz(ideb))./deltaz;

%  transforme l'equation ifin:
   tridiag(ifin,1) = tridiag(ifin,1) ...
                  - tridiag(ifin,3).*(zz(ifin-1) - zz(ideb-1))./deltaz;
   tridiag(ideb,3) = tridiag(ifin,3).*(zz(ifin) - zz(ifin-1))./deltaz;
end

%  !!!!!!!!!!!!!!!!!!!!   Partie a verifier - fin !!!!!!!!!!!!!!!!!!!!!!!!!!!!

%
%  =======================================================
%      Calcul des valeurs propres et vecteurs propres
%  =======================================================

%  Methode bestiale: Fournit la matrice complete:

nzok  = length(isok)+2;
A_mat = zeros(nzok,nzok);

%  Niveaux internes:
%  ----------------
for i = 2:nzok-1
   A_mat(i,i-1) = tridiag(isok(i-1),1);
   A_mat(i,i)   = tridiag(isok(i-1),2);
   A_mat(i,i+1) = tridiag(isok(i-1),3);
end

%  Condition limite au fond:
%  ------------------------
%   gn = 0 ==> on suprime dernier ligne et derniere colonne:
nzok  = nzok - 1;
A_mat = A_mat(1:nzok,1:nzok);


%  Condition limite en surface:
%  ---------------------------
if iop == 0
   A_mat = A_mat(2:nzok,2:nzok);
   nzok  = nzok - 1;
else
%  !!!!!!!!!!!!!!!!!!!!   Partie a verifier - debut !!!!!!!!!!!!!!!!!!!!!!!!!!!
   A_mat(1,1) = 1.0/(grav*dz1(1));
   A_mat(1,2) = -A_mat(1,1);
%  !!!!!!!!!!!!!!!!!!!!   Partie a verifier - fin   !!!!!!!!!!!!!!!!!!!!!!!!!!!
end

[RR, RI, VR, VI, niter, ifail] = f02agf(A_mat);
[EV,II] = sort(RR);
Isel = II(1:nmod);
Esel = EV(1:nmod);
clear A_mat


%  Extrait les modes de deplacement vertical:
%  -----------------------------------------
if iop == 0
   gmod = [zeros(1,nmod); VR(:,Isel); zeros(1,nmod)] ;
   z0   = [zz(1); zz(isok); zz(nz)];
else
   gmod = [VR(:,Isel); zeros(1,nmod)];
   z0   = [zz(isok); zz(nz)];
end

%  Repasse sur la grille complete:
%  ------------------------------
if nzeros>0,
   bid = interp1(z0,gmod,zz);
   gmod = bid;
end
clear bid VR


%  Calcule les modes fonction courant: f(z) = -gz(z)/ev
%  ----------------------------------------------------
fmod = zeros(nz,nmod);
der = gmod(2:nz,:) - gmod(1:nz-1,:);
der = der./(dz1*ones(1,nmod));
der = (der(1:nz-2,:) + der(2:nz-1,:))*0.5;
fmod(2:nz-1,:) = -der./(ones(nz-2,1)*Esel');
fmod(1,:)      = -((gmod(2,:)  - gmod(1,:))/dz1(1))./Esel';
fmod(nz,:)     = -((gmod(nz,:) - gmod(nz-1,:))/dz1(nz-1))./Esel';

%  Calcule les rayons de Rossby:
%  ----------------------------
rossby = (1/fcor)*ones(nmod,1)./sqrt(Esel);



