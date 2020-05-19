function [y]=lanczos(x,fn,nj)
 
% Filtre de lanczos :
 
% y  : output filtered data
% x  : input data
% fn : the cut off frequency (any k/n between 1/n , 0.5)
% nj : (number of points in the filter-1)/2 ; transition width of the 
%      transfer function is 1/nj.
 
% The filter uses lanczos smoothing of a square box transfert function.
 
% Auteur : C. Hemon a partir source fortran Herle Mercier.


% -----------------------------------------------------------
% version modifiee pour application Bathysonde par C.Lagadec
% en oct. 98 (non inversion de tableaux)
% -----------------------------------------------------------
 

 
% Test si bon appel a la fonction.
 
	if nargin ~= 3
		error('Il faut 3 parametres en entree de la fonction lanczos')
	end

	nk = nj+1;

	n=length(x);
 
% Pre-allocation memoire.
 
	y=zeros(1,n);
	a=zeros(1,nk);
	b=zeros(1,nk);
	ex=zeros(1,nk);
	ey=zeros(1,nk);
	pds=zeros(1,nk);
 
% Inverse vecteur en entree (on travaille sur vecteur colonne).

% modif C.Lagadec (oct 98)
% pas d'inversion car vecteur x deja vecteur colonne
%	x=x' ;
 
% Calcul des poids du filtre.

	pds(1)= 2*fn;
	i=2:nk;
	a(i)=pi*(i-1);
	b(i)=2*a(i)*fn;
	ex(i)=(sin(b(i))./a(i));
	ey(i)=a(i)/nj;
	pds(i)=ex(i).*sin(ey(i))./ey(i);
	
 
% Somme des poids.
 
	ss=pds(1)+sum(2*pds(2:nk));

%
% Filtrage des nj premiers points.
%
	for i=1:nj
		den=0;
		k2=2*i-1;
		j=1:k2;
		val=abs(i*ones(size(j))-j)+1;
		y(i)=sum(pds(val).*x(j));
		den=sum(pds(val));
		y(i)=y(i)./den;
	end

%
% Filtrage du centre de la serie.
%
	for i=nk:n-nj
		k1=i-nj;
		k2=i+nj;
		j=k1:k2;
		val=abs(i*ones(size(j))-j)+1;
		y(i)=sum(pds(val).*x(j));
		y(i)=y(i)/ss;
	end

%
% Filtrage des nj derniers points de la serie.
%
	for i=n-nj+1:n
		k1=2*i-n;
		j=k1:n;
		den=0;
		val=abs(i*ones(size(j))-j)+1;
		y(i)=sum(pds(val).*x(j));
		den=sum(pds(val));
		y(i)=y(i)./den;
	end

%
% Inverse vectuer en sortie (pour obtenir vecteur ligne).
% Modif C.Lagadec (oct. 98) :
% Pas d'inversion de matrice
%	y=y';	
	
