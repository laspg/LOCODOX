function [refID,lon,lat,sta_ctd,juld,presi,tempi,theta0i,psali,sig0i,doxyi,o2ok] = rd_O2data_hydro_LOPS(ficcamp,datatype )
%
% Fonction qui lit les donnees oxygene calibrees (campagne et bouteille)
% pour calibrer les donnees flotteurs
% Input:
%        ficcamp : nom du fichier HYDRO (avec chemin complet)
%        datatype = type de donnees lues
%                 = 'ctd' ou 'btl'
%

presi = 0:10:4000;
        
        if strcmp(datatype,'ctd') == 1
            % Lecture des donnees campagne CTD
            o2ok = 1;

            [refID, lat,lon,juld,platform,temp,theta0,psal,doxy,pres,sig0,sta_ctd] = ...
                read_nc_base_hydro(ficcamp,'TEMP','TPOT','PSAL','OXYK','PRES','SIG0','STATION_NUMBER');
            temp(temp == -9999) = NaN;
            theta0(theta0 == -9999) = NaN;
            psal(psal == -9999) = NaN;
            doxy(doxy == -9999) = NaN;
            pres(pres == -9999) = NaN;
            sig0(pres == -9999) = NaN;
            sta_ctd = sta_ctd';
            if ~isempty(lon)
               for i = 1:length(lon)
                 iok = find(isfinite(pres(i,:)));
                 tempi(i,:) = interp1(pres(i,iok),temp(i,iok),presi);
                 theta0i(i,:) = interp1(pres(i,iok),theta0(i,iok),presi);
                 psali(i,:) = interp1(pres(i,iok),psal(i,iok),presi);
                 sig0i(i,:) = interp1(pres(i,iok),sig0(i,iok),presi);
                 doxyi(i,:) = interp1(pres(i,iok),doxy(i,iok),presi);
               end
            end
           
        elseif strcmp(datatype,'btl') == 1
            % Lecture des donnees campagne bouteilles
            o2ok = 0;
            pres = NaN;
            temp = NaN;
            theta0 = NaN;
            psal = NaN;
            sig0 = NaN;
            oxyl = NaN;
            doxy = NaN;
        end
        
    


