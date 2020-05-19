%
% fonction de conversion
% 
% C. Lagadec - mars 2002
%
% fonction écrite à partir du sous-programme Fortran
% de la chaîne Hydro
%
%
 function [ax,ay] = conversion(xlono,xlato,xlon1,xlat1);

 pi=4.*atan(1.);

 ty=(xlato+xlat1)/2.;

 ax=(xlon1-xlono)*cos(pi*ty/180.)*60.*1.852;

 ay=(xlat1-xlato)*60.*1.852;






