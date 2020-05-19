function [refID,lon,lat,sta_ctd,juld,pres,temp,theta0,psal,sig0,doxy,o2ok] = rd_O2data_ref_ovid11( campyear,datatype )
%
% Fonction qui lit les donnees oxyg??ne calibrees (campagne et bouteille)
% pour calibrer les donnees flotteurs
% Input:
%        campyear = annee des donnees de reference
%                 = ov10th/ov11me/ov11di/gh10md/ov12sg/geov14
%        datatype = type de donnees lues
%                 = 'ctd' ou 'btl'
%
switch campyear
    
case 'ov11me'
        if strcmp(datatype,'ctd') == 1
            % Lecture des donnees campagne CTD
            o2ok = 1;
            
            tabnum = [1901209;
                1901210;
                1901211;
                1901212;
                1901213;
                1901214;
                1901215;
                1901217;
                1901218];
            
            
            tabo2file = ['m852_003.ctd';
                'm852_027.ctd';
                'm851_014.dat';
                'm851_018.dat';
                'm851_026.dat';
                'm851_079.dat';
                'm851_082.dat';
                'm851_095.dat';
                'm851_106.dat'];
            
            nfile = size(tabo2file,1);
            inpathlpo = '/home4/homedir4/perso/vthierry/';
            inpathmac = '/Users/vthierry/';
            if exist(inpathlpo,'dir') == 7
                inpath = fullfile(inpathlpo,'PROJETS/OXYGEN/ACHATS/2010/OXYBTL4CALIB/');
            elseif exist(inpathmac, 'dir') == 7
                inpath = '/Users/vthierry/PROJETS/OXYGEN/ACHATS/2010/OXYBTL4CALIB/';
            else
                disp('Pas de donnees de reference');
                stop
            end
            lat = NaN(nfile,1);
            lon = NaN(nfile,1);
            juld = NaN(nfile,1);
            pres = NaN(nfile,2100);
            temp = NaN(nfile,2100);
            psal = NaN(nfile,2100);
            doxy = NaN(nfile,2100);
            sig0 = NaN(nfile,2100);
            sta_ctd = (1:nfile)';
            refID = cell(nfile,1);
            
            for ifile = 1:nfile
                fname = tabo2file(ifile,:);
                filename = fullfile(inpath,num2str(tabnum(ifile)),fname);
                switch fname(1:4)
                    case 'm852'
                        [latf,lonf,juldf,presf,dephf,tempf,psalf,doxyovf,sig0f] = rd_m852(filename);
                        refID{ifile} = sprintf('ov11me_%d',ifile);
                        lat(ifile) = latf;
                        lon(ifile) = lonf;
                        juld(ifile) = juldf;
                        pres(ifile,1:min(length(presf),2100)) = presf(1:min(length(presf),2100));
                        %deph(ifile,1:min(length(presf),2100)) = dephf(1:min(length(presf),2100));
                        temp(ifile,1:min(length(presf),2100)) = tempf(1:min(length(presf),2100));
                        psal(ifile,1:min(length(presf),2100)) = psalf(1:min(length(presf),2100));
                        doxy(ifile,1:min(length(presf),2100)) = doxyovf(1:min(length(presf),2100));
                        sig0(ifile,1:min(length(presf),2100)) = sig0f(1:min(length(presf),2100));
%                         theta0 = pottemp(psal,temp,pres,0);
                        theta0 = sw_ptmp(psal,temp,pres,0);
                        
                        
                    case 'm851'
                        [latf,lonf,juldf,presf,tempf,theta0f,psalf,doxyovf,sig0f] = rd_m851(filename);
                        refID{ifile} = sprintf('ov11me_%d',ifile);
                        lat(ifile) = latf;
                        lon(ifile) = lonf;
                        juld(ifile) = juldf;
                        pres(ifile,1:min(length(presf),2100)) = presf(1:min(length(presf),2100));
                        temp(ifile,1:min(length(presf),2100)) = tempf(1:min(length(presf),2100));
                        psal(ifile,1:min(length(presf),2100)) = psalf(1:min(length(presf),2100));
                        doxy(ifile,1:min(length(presf),2100)) = doxyovf(1:min(length(presf),2100));
                        sig0(ifile,1:min(length(presf),2100)) = sig0f(1:min(length(presf),2100));
%                         theta0 = pottemp(psal,temp,pres,0);
                        theta0 = sw_ptmp(psal,temp,pres,0);
                        
                    case 'disc'
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
        
    case 'ov11di'
        if strcmp(datatype,'ctd') == 1
            o2ok = 1;
            % Lecture des donnees campagne CTD
            inpath = '/home4/homedir4/perso/vthierry/PROJETS/OXYGEN/ACHATS/2010/DISCOVERY/DATA/';
            pres = NaN(29,6000);
            presqc = NaN(29,6000);
            temp = NaN(29,6000);
            tempqc = NaN(29,6000);
            psal = NaN(29,6000);
            psalqc = NaN(29,6000);
            doxy = NaN(29,6000);
            doxyqc = NaN(29,6000);
            sta_ctd = NaN(29,1);
            lat = NaN(29,1);
            lon = NaN(29,1);
            juld = NaN(29,1);
            refID = cell(29,1);
            
            for iprf = 1:29
                fname = sprintf('a16n2011_000%2.2d_00001_ct1.csv',iprf);
                fid = fopen(fullfile(inpath,fname),'r');
                sta_ctd(iprf) = iprf;
                refID{iprf} = sprintf('ov11di_%d',iprf);
                nl = 1;
                
                % Get metadata
                headers = cell(60,1);
                headers{nl} = fgetl(fid);
                while ischar(headers{nl})
                    nl = nl+1;
                    tmp = fgetl(fid);
                    if tmp(1) == ' '
                        break
                    else
                        headers{nl} = tmp;
                    end
                end
                headers = headers(~cellfun('isempty',headers));
                
                % get useful variables
                myvar = {'DATE','TIME', 'LATITUDE','LONGITUDE'};
                for v= 1:length(myvar)
                    tmp = headers{~cellfun('isempty',strfind(headers,myvar{v}))};
                    H.(myvar{v}) = char(regexp(tmp,'(?:-?\d+.\d*)','match'));
                end
                juld(iprf) = datenum(str2double(H.DATE(1:4)),str2double(H.DATE(5:6)),str2double(H.DATE(7:8)));
                if juld(iprf) >= 700000
                    juld(iprf) = juld(iprf) - datenum([1950 01 01 00 00 00]);
                end
                lat(iprf)  = str2double(H.('LATITUDE'));
                lon(iprf) = str2double(H.('LONGITUDE'));
                
                % get data
                frewind(fid)
                data = textscan(fid,'%f %d %f %d %f %d %f %d %f','headerlines',15, 'delimiter',',');
                pres(iprf,1:length(data{1})) = data{1};
                presqc(iprf,1:length(data{1})) = data{2};
                temp(iprf,1:length(data{1})) = data{3};
                tempqc(iprf,1:length(data{1})) = data{4};
                psal(iprf,1:length(data{1})) = data{5};
                psalqc(iprf,1:length(data{1})) = data{6};
                doxy(iprf,1:length(data{1})) = data{7};
                doxyqc(iprf,1:length(data{1})) = data{8};                
            end
            pres(pres == -999.0 | presqc == 9) = NaN;
            temp(temp == -999.0 | presqc == 9) = NaN;
            psal(psal == -999.0 | presqc == 9) = NaN;
            doxy(doxy == -999.0 | presqc == 9) = NaN;
            theta0 = sw_ptmp(psal,temp,pres,0);
            sig0 = sw_pden(psal,temp,pres,0);

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

   
end


