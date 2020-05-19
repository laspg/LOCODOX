%
%
%------------------------------------------------------------------
 
%     subroutine vorpot(bv2,n,xlat,vrp)
 
%------------------------------------------------------------------
%
%     creation : Fevrier 90                      auteur : D.Jacolot
%     --------                                   ------
%
%     objet : Calcul de la vorticite potentielle GE  
%     -----  
% 
%     description :
%     -----------
%                                                  
%     entree : 
%     ------   
%        bv2  : tableau des frequences de Brunt-Vaisala au 
%               carre rad**2/s**2
%        n    : nombre de points 
%        xlat : latitude du point ou est calculee la vorticite
%               potentielle  (degres decimaux)
%  
%     sortie : 
%     ------   
%        vrp  : tableau des vorticites potentielles (10J/kg)
%
%     sous-programmes appeles : 
%     -----------------------  
%     gravi8, f08
%
%-------------------------------------------------------------------
%
%  Modif C.Lagadec  (mars 2000) 
%  test si la valeur de bv2 est à 1.e+36 => pas de calcul de vrp
%
%-------------------------------------------------------------------

   function vrp = vorpot(bv2,xlat)

   valdef = 1.e+36;

   coriolis = f08(xlat);
   gsf = coriolis/ gravit(xlat);

% modif bizarre mais le == ne marche pas !!!!
   inan = find(bv2>1.e+35);
%  inan = find(bv2==valdef);

   vrp = gsf * bv2;
   vrp(inan) = valdef;

   clear coriolis gsf inan valdef;

