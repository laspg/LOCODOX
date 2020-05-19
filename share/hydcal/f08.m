      function coriolis = f08(xlat);

% ----------------------------------------------------------------
%
%    Fonction f08 : C.Lagadec - Janv.99
%
%    objet : Calcul du parametre de Coriolis.
%             reecrit a partir du SP f08 de la chaine Hydro(armen)
%     xlat : Latitude du point ou l'acceleration de gravite
%                sera calculee (degres decimaux)

      degrad = pi/180.0;
      alat = xlat*degrad;
      slat = sin(alat);
      coriolis = slat*degrad/120.0;

