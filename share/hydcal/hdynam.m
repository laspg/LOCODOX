      function [hdy, messhdy] =  hdynam(p,sig,si35)

% -----------------------------------------------------------------
%     Creation : Fevrier 90                      auteur : D.Jacolot
%     --------                                   ------
%
%     Modification :                             par :
%     ------------                               ---
%     - corrections                        09/10/90    D.Jacolot
%
%     - correction :                       15/09/92    C.Lagadec
%       Ajout du signe negatif
%       pour le calcul de fideb et fi
%
%     Remarque : 
%     ----------
%      Pb rencontre avec donnees de CONVIV :
%      forte salinite => sigma > sigma 35 
%                     => signe de hdy inverse
%
%     Objet : Calcul des hauteurs dynamiques  
%     -----  
% 
%     Description :
%     -----------
%                                                  
%     entree : 
%     ------   
%        p    : tableau des pressions (db)
%        sig  : tableau des anomalies de densite in situ (kg/m**3)
%        si35 :  tableau des anomalies de densite in situ 
%                pour t = 0 et s = 35 (kg/m*3)
%  
%     sortie : 
%     ------   
%        hdy  : tableau des hauteurs dynamiques en mdyn (10J/kg)
%
% -----------------------------------------------------------------
%    Version MATLAB :
%    ----------------
%
%  Calcul de la hauteur dynamique
%  fonction Matlab créée à partir du sous-programme fortran HDYNAM 
%  de la chaine Hydro
% 
%  C.Lagadec - janv.99
% 
%   messhdy est initialise quand probleme de calcul
%
%------------------------------------------------------------------
 
% pression de reference 
 global xpr;

       messhdy  = [];
       hdy      = [];


       n    = length(p);
       if  p(n) < xpr
             mess= sprintf('%s%d%s',' La pression maximale du fichier (',p(n),'), inférieure a la pression de référence, sera la pression en compte pour les calculs. ');
             h = warndlg(mess,'Attention !');
             waitfor(h);
       end;
       xpr  = min(p(n),xpr);
       iref = find(p==xpr);

       if   (isempty(iref))
           messhdy = sprintf('%d%s',xpr,' : Pression de référence inconnue dans le tableau des pressions');
           hdy = [];

               else

           hdy(iref) = 0.;
           fideb = - (si35(iref)-sig(iref))/(1000.+sig(iref)+si35(iref)+1.e-03*sig(iref)*si35(iref));
           fi = fideb;

% ------   calcul de la partie supérieure ----

           for i = iref+1:n
              fim1 = fi;
              fi   = - (si35(i)-sig(i))/(1000. + sig(i) + si35(i) + 1.e-03*sig(i)*si35(i));

              hdy(i) = hdy(i-1) + (fi+fim1)*(p(i)-p(i-1))*0.5;
           end

% ------   calcul de la partie inférieure ----
          
           fi = fideb;
  
           for  i = iref-1:-1:1
              fip1 = fi;
              fi   = - (si35(i)-sig(i))/(1000. + sig(i) + si35(i) + 1.e-03*sig(i)*si35(i));
              hdy(i) = hdy(i+1) - (fi+fip1)*(p(i+1)-p(i))*0.5;
           end
   end

          
