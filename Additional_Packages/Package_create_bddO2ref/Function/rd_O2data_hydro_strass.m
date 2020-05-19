function [refID,lon,lat,sta_ctd,juld,pres,temp,theta0,psal,sig0,doxy,o2ok] = rd_O2data_ref_strass(campyear,datatype)
%
% Fonction qui lit les donnees oxyg??ne de la campagne strasse
% name='strasse'
%datatype='btl'
%
        o2ok = 1;
        % chargement des donnees realisees au deploiement
        strasse = load('/Users/vthierry/matlab/ARGO_DEEP/DEEP_ARVOR/CTDeduitup.017');
        strasse = flipud(strasse);
        pres = strasse(:,2)';
        temp = strasse(:,3)';
        psal = strasse(:,4)';
        doxy = strasse(:,5)';
        %theta0 = sw_ptmp(psal,temp,pres,0);
        %sig0 = sw_pden(psal,temp,pres,0)-1000;
        lon = -35.700;
        lat = 25.98;
        sta_ctd = 1;
        fprintf('pas de jour julien pour les donnees strasse\n');
        juld = datenum(2012,06,15)-datenum(1950,1,1,0,0,0);
        refID = sprintf('%s_1',campyear);
        
        fprintf('correction en salinte (-0.004) et pression (-1.6) des donnees strasse\n')
        presstrassei = pres -1.6;
        tempstrassei = temp;
        psalstrassei = psal - 0.004;
        doxystrassei = doxy;
        
        % Reduction des donnees ?? 5db
        pas = 5;
        pres = [5:pas:4000];
        temp = NaN(1,length(pres));
        psal = NaN(1,length(pres));
        doxy = NaN(1,length(pres));
        for ipres = 1:length(pres)
            iok = (presstrassei > pres(1,ipres)-0.5*pas & presstrassei <=  pres(1,ipres)+0.5*pas);
            %refpres(ipres,1) = refpres(ipres);%mean(presstrassei(iok));
            temp(1,ipres) = mean(tempstrassei(iok));
            psal(1,ipres) = mean(psalstrassei(iok));
            doxy(1,ipres) = mean(doxystrassei(iok));
        end
        
        theta0 = sw_ptmp(psal,temp,pres,0);
        sig0 = sw_pden(psal,temp,pres,0)-1000;
        
   
        