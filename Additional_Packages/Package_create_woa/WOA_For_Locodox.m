clear all;
close all;

%%% création de "WOA_monthly_5500_1deg.nc" pour LOCODOX
%%% depuis 2018 WOA fournit ses données d'oxygène en µmol/kg, la conversion
%%% n'est donc plus nécessaire. 
%%% pour matcher avec LOCODOX mais pas remplies
%%% MENSUEL : depth 1:1500 puis ANNUEL : 1:5500
%%% Virginie Racape:    13/01/2020
%%% Thierry Reynaud     22/03/2020

% Define directories and file name structures
DIR=define_dir;

% Reading Annual files ==> 0 - 5500m
CLIM_A=read_annual_clim(DIR);

% Reading Monthly files ==> 0 - 1500m
CLIM_M=read_monthly_clim(DIR);

% Build the monthly climatology 0 - 5500m
CLIM_FM=build_full_monthly(CLIM_A,CLIM_M);
clear CLIM_A CLIM_M

% Compute the monthly Pressure + potential density Climatologies 0 - 5500m 
CLIM_FM=cal_pressure(CLIM_FM);

% Create netcdf outputfile
create_woa_netcdf(DIR,CLIM_FM)

%<===============================================>
function create_woa_netcdf(DIR,CLIM_FM)

% creation of the WOA 2018 netcdf file for Locodox

% Remove existing file
dest=DIR.dest;
if exist(dest)
  delete(dest);
end

% variables de sortie

lonwoa=CLIM_FM.lonwoa;
ibeg=find(lonwoa>=0);
iend=find(lonwoa<0);

nLon=CLIM_FM.nLon;
nLat=CLIM_FM.nLat;
nZ=CLIM_FM.nZ;
nT=CLIM_FM.nT;

% Writting latwoa
latwoa=CLIM_FM.latwoa;
nccreate(dest,'latitude','Format','netcdf4_classic','Dimensions',{'lat' nLat});
ncwrite(dest,'latitude', latwoa);
ncwriteatt(dest,'latitude','long_name','latitude');
ncwriteatt(dest,'latitude','units','^(o) N');
clear latwoa

% Writting lonwoa2
lonwoa2=[lonwoa(ibeg);360+lonwoa(iend)];
nccreate(dest,'longitude','Format','netcdf4_classic','Dimensions',{'lon' nLon});
ncwrite(dest,'longitude', lonwoa2);
ncwriteatt(dest,'longitude','long_name','longitude');
ncwriteatt(dest,'longitude','units','^(o) E');
clear lonwoa lonwoa2

Depth_a=CLIM_FM.Depth;
% Writting Depth_a : 102 vertical levels
nccreate(dest,'Depth','Format','netcdf4_classic','Dimensions',{'Depth' nZ});
ncwrite(dest,'Depth', Depth_a);
ncwriteatt(dest,'Depth','long_name','Depth');
ncwriteatt(dest,'Depth','units','m');
clear Depth_a


% Writting timewoa2 : day numbers of the current year
timewoa2=datenum(2019,1:12,15)-datenum(2019,1,1)+1;
timewoa2(3:end)=timewoa2(3:end)+1;
nccreate(dest,'time','Format','netcdf4_classic','Dimensions',{'time' nT});
ncwrite(dest,'time', timewoa2);
ncwriteatt(dest,'time','long_name','time');
ncwriteatt(dest,'time','units','day');
clear timewoa2

% Writting O2_woa_full : Dissolved Oxygen in mumol/Kg
O2_woa_full=CLIM_FM.O2_woa_full;
tmp=NaN*ones(nLon,nLat,nZ,12);
tmp(1:size(ibeg,1),:,:,:)=O2_woa_full(ibeg,:,:,:);
tmp(size(ibeg,1)+1:end,:,:,:)=O2_woa_full(iend,:,:,:);
nccreate(dest,'doxywoa','Format','netcdf4_classic','Dimensions',{'lon' nLon 'lat' nLat 'Depth' nZ 'time' nT});
ncwrite(dest,'doxywoa', tmp);
ncwriteatt(dest,'doxywoa','long_name','DISSOLVED OXYGEN (mumol/kg) from WOA 2018');
ncwriteatt(dest,'doxywoa','units','mumol/kg');
clear tmp O2_woa_full

% nccreate(dest,'doxywoa_sd','Format','netcdf4_classic','Dimensions',{'lon' size(m_an,1) 'lat' size(m_an,2) 'Depth' size(m_an,3) 'time' size(m_an,4)});
% ncwrite(dest,'doxywoa_sd', m_sd);
% ncwriteatt(dest,'doxywoa_sd','long_name','DISSOLVED OXYGEN (umol/kg) from WOA 2018');
% ncwriteatt(dest,'doxywoa_sd','units','umol/kg');

% Writting Psat_woa_full : Oxygen Saturation in percent
Psat_woa_full=CLIM_FM.Psat_woa_full;
tmp=NaN*ones(nLon,nLat,nZ,12);
tmp(1:size(ibeg,1),:,:,:)=Psat_woa_full(ibeg,:,:,:);
tmp(size(ibeg,1)+1:end,:,:,:)=Psat_woa_full(iend,:,:,:);
nccreate(dest,'Psatwoa','Format','netcdf4_classic','Dimensions',{'lon' nLon 'lat' nLat 'Depth' nZ 'time' nT});
ncwrite(dest,'Psatwoa', tmp);
ncwriteatt(dest,'Psatwoa','long_name','%saturation from WOA 2018');
ncwriteatt(dest,'Psatwoa','units','%');
clear tmp Psat_woa_full

% nccreate(dest,'Psatwoa_sd','Format','netcdf4_classic','Dimensions',{'lon' size(m_an,1) 'lat' size(m_an,2) 'Depth' size(m_an,3) 'time' size(m_an,4)});
% ncwrite(dest,'Psatwoa_sd', Psatwoa_sd);
% ncwriteatt(dest,'Psatwoa_sd','long_name','%saturation from WOA 2018');
% ncwriteatt(dest,'Psatwoa_sd','units','%');

% Writting Pdens_woa_full : Surface Referenced Potential Density
Pdens_woa_full=CLIM_FM.Pdens_woa_full;
tmp=NaN*ones(nLon,nLat,nZ,12);
tmp(1:size(ibeg,1),:,:,:)=Pdens_woa_full(ibeg,:,:,:);
tmp(size(ibeg,1)+1:end,:,:,:)=Pdens_woa_full(iend,:,:,:);
nccreate(dest,'density','Format','netcdf4_classic','Dimensions',{'lon' nLon 'lat' nLat 'Depth' nZ 'time' nT});
ncwrite(dest,'density', tmp);
ncwriteatt(dest,'density','long_name','density from WOA 2018');
ncwriteatt(dest,'density','units','kg/m3');
clear tmp Pdens_woa_full

% Writting P_woa_full : Pressure
P_woa_full=CLIM_FM.P_woa_full;
tmp=NaN*ones(nLon,nLat,nZ,12);
tmp(1:size(ibeg,1),:,:,:)=P_woa_full(ibeg,:,:,:);
tmp(size(ibeg,1)+1:end,:,:,:)=P_woa_full(iend,:,:,:);
nccreate(dest,'preswoa','Format','netcdf4_classic','Dimensions',{'lon' nLon 'lat' nLat 'Depth' nZ 'time' nT});
ncwrite(dest,'preswoa',tmp);
ncwriteatt(dest,'preswoa','long_name','pressure from WOA 2018 ');
ncwriteatt(dest,'preswoa','units','dbar');
clear tmp P_woa_full

%Writting S_woa_full : Salinity
S_woa_full=CLIM_FM.S_woa_full;
tmp=NaN*ones(nLon,nLat,nZ,12);
tmp(1:size(ibeg,1),:,:,:)=S_woa_full(ibeg,:,:,:);
tmp(size(ibeg,1)+1:end,:,:,:)=S_woa_full(iend,:,:,:);
nccreate(dest,'PSAL_WOA','Format','netcdf4_classic','Dimensions',{'lon' nLon 'lat' nLat 'Depth' nZ 'time' nT});
ncwrite(dest,'PSAL_WOA',tmp);
ncwriteatt(dest,'PSAL_WOA','long_name','PRATICAL SALINITY from WOA 2018 ');
ncwriteatt(dest,'PSAL_WOA','units','psu');
clear tmp S_woa_full

% Writting T_woa_full : Temperature
T_woa_full=CLIM_FM.T_woa_full;
tmp=NaN*ones(nLon,nLat,nZ,12);
tmp(1:size(ibeg,1),:,:,:)=T_woa_full(ibeg,:,:,:);
tmp(size(ibeg,1)+1:end,:,:,:)=T_woa_full(iend,:,:,:);
nccreate(dest,'TEMP_WOA','Format','netcdf4_classic','Dimensions',{'lon' nLon 'lat' nLat 'Depth' nZ 'time' nT});
ncwrite(dest,'TEMP_WOA',tmp);
ncwriteatt(dest,'TEMP_WOA','long_name','SEA TEMPERATURE from WOA 2018 ');
ncwriteatt(dest,'TEMP_WOA','units','degree_Celsius');
clear tmp T_woa_full
clear nLat nLon nT nZ

return
end

function CLIM_FM=cal_pressure(CLIM_FM)
% Calculating the Pressure profile for the whole water column
% from individual S and T profiles (0 to 5500m).
% The potential density fields are here calculated

nLon=CLIM_FM.nLon;
nLat=CLIM_FM.nLat;
nZ=CLIM_FM.nZ;
nT=CLIM_FM.nT;

P_woa_full=NaN*ones(nLon,nLat,nZ,nT);%Pressure
Pdens_woa_full=NaN*ones(nLon,nLat,nZ,nT);%Pressure

S_woa_full=CLIM_FM.S_woa_full;
T_woa_full=CLIM_FM.T_woa_full;
Depth=CLIM_FM.Depth;
latwoa=CLIM_FM.latwoa;

tic;
for i=1:nLon
    for j=1:nLat
        for k=1:nT
            if i==50 && j==50
                display(['Computing pressure field for month: ',num2str(k)]);
            end
            iter=0;
            crit=100;
            s=squeeze(S_woa_full(i,j,:,k));
            t=squeeze(T_woa_full(i,j,:,k));
            idx=~isnan(s);
            z=Depth(idx);s=s(idx);t=t(idx);
            if sum(idx) >0
                while crit>0.0001
                    iter=iter+1;
                    %display(['Iter= ',num2str(iter)]);
                    if (iter==1)
                        p=abs(z);
                    end
                    sigma=sw_dens(s,t,p)-1000;
                    p_old=p;
                    p=zenprs(-abs(z),sigma,latwoa(j));
                    crit=max(abs(p-p_old));
                end
                P_woa_full(i,j,1:length(p),k)=p;
                rp=sw_pden(s,t,p,0);
                Pdens_woa_full(i,j,1:length(rp),k)=rp;
            end
        end
        
    end
end
toc;
CLIM_FM.P_woa_full=P_woa_full;
CLIM_FM.Pdens_woa_full=Pdens_woa_full;

return
end

function CLIM_FM=build_full_monthly(CLIM_A,CLIM_M)

% Building the whole water column monthly climatologies
% 0 to 5500m from the monthly (0 - 1500m) and annual (0 to 5500m)
% climatologies.

nLon=CLIM_A.nLon;
nLat=CLIM_A.nLat;
nZ=CLIM_A.nZ;
nZ_top=CLIM_M.nZ_top;
nT=CLIM_M.nT;

CLIM_FM.latwoa=CLIM_A.latwoa;
CLIM_FM.lonwoa=CLIM_A.lonwoa;
CLIM_FM.Depth=CLIM_A.Depth_a;
CLIM_FM.nLon=nLon;
CLIM_FM.nLat=nLat;
CLIM_FM.nT=nT;
CLIM_FM.nZ=nZ;

CLIM_FM.T_woa_full=NaN*ones(nLon,nLat,nZ,nT);%Temperature Insitu
CLIM_FM.S_woa_full=NaN*ones(nLon,nLat,nZ,nT);%Salinity
CLIM_FM.Pdens_woa_full=NaN*ones(nLon,nLat,nZ,nT);%Insitu Density

CLIM_FM.O2_woa_full=NaN*ones(nLon,nLat,nZ,nT);%O2
CLIM_FM.O2_sd_woa_full=NaN*ones(nLon,nLat,nZ,nT);%O2 standard deviation

CLIM_FM.Psat_woa_full=NaN*ones(nLon,nLat,nZ,nT);%O2 sat
CLIM_FM.Psat_sd_woa_full=NaN*ones(nLon,nLat,nZ,nT);%O2 sat standard deviation


for i=1:12
    CLIM_FM.T_woa_full(:,:,1:nZ_top,i)=CLIM_M.T_woa(:,:,1:nZ_top,i);
    CLIM_FM.T_woa_full(:,:,nZ_top+1:nZ,i)=CLIM_A.T_woa_a(:,:,nZ_top+1:nZ);
    
    CLIM_FM.S_woa_full(:,:,1:nZ_top,i)=CLIM_M.S_woa(:,:,1:nZ_top,i);
    CLIM_FM.S_woa_full(:,:,nZ_top+1:nZ,i)=CLIM_A.S_woa_a(:,:,nZ_top+1:nZ);
    
    CLIM_FM.O2_woa_full(:,:,1:nZ_top,i)=CLIM_M.oxywoa(:,:,1:nZ_top,i);
    CLIM_FM.O2_woa_full(:,:,nZ_top+1:nZ,i)=CLIM_A.oxywoa_a(:,:,nZ_top+1:nZ);
    
    CLIM_FM.O2_sd_woa_full(:,:,1:nZ_top,i)=CLIM_M.oxywoa_sd(:,:,1:nZ_top,i);
    CLIM_FM.O2_sd_woa_full(:,:,nZ_top+1:nZ,i)=CLIM_A.oxywoa_sd_a(:,:,nZ_top+1:nZ);
    
    CLIM_FM.Psat_woa_full(:,:,1:nZ_top,i)=CLIM_M.psatwoa(:,:,1:nZ_top,i);
    CLIM_FM.Psat_woa_full(:,:,nZ_top+1:nZ,i)=CLIM_A.psatwoa_a(:,:,nZ_top+1:nZ);
    
    CLIM_FM.Psat_sd_woa_full(:,:,1:nZ_top,i)=CLIM_M.psatwoa_sd(:,:,1:nZ_top,i);
    CLIM_FM.Psat_sd_woa_full(:,:,nZ_top+1:nZ,i)=CLIM_A.psatwoa_sd_a(:,:,nZ_top+1:nZ);
end


return
end

function CLIM_M=read_monthly_clim(DIR)

% Read the WOA 2018 Monthly climatologies (0 to 1500m)

oxyfile=DIR.oxyfile;
tfile=DIR.tfile;
sfile=DIR.sfile;
satfile=DIR.satfile;

% % Variables atlas WOA18 mensuel    360   180  57 12
depth_m = double(ncread([oxyfile,'01_01.nc'],'depth'));

for i = 1:12
    mois=sprintf('%2.2i',i);
    
    tmp=double(ncread([tfile,mois,'_01.nc'],'time'));
    scale_units = ncreadatt([tfile,mois,'_01.nc'],'time','units');
    if strfind(scale_units,'months')
        ref=sscanf(scale_units, 'months since %f-%f-%f %f:%f:%f');
        ref=datenum(ref');
        time_t=ref+365.25*tmp/12;
        display(['Extracting month number: ',mois ]);
        display(['Temperature file time:        ',datestr(time_t,'dd/mm/yyyy HH:MM:SS')]);
    else
        display(['Error decoding Temp time:     ',scale_units,'--> ','months']);
    end
    
    tmp=double(ncread([sfile,mois,'_01.nc'],'time'));
    scale_units = ncreadatt([tfile,mois,'_01.nc'],'time','units');
    if strfind(scale_units,'months')
        ref=sscanf(scale_units, 'months since %f-%f-%f %f:%f:%f');
        ref=datenum(ref');
        time_s=ref+365.25*tmp/12;
        display(['Salinity file time   :        ',datestr(time_s,'dd/mm/yyyy HH:MM:SS')]);
    else
        display(['Error decoding Salinity time: ',scale_units,'--> ','months']);
    end
    
    tmp=double(ncread([oxyfile,mois,'_01.nc'],'time'));
    scale_units = ncreadatt([oxyfile,mois,'_01.nc'],'time','units');
    if strfind(scale_units,'months')
        ref=sscanf(scale_units, 'months since %f-%f-%f %f:%f:%f');
        ref=datenum(ref');
        time_o2=ref+365.25*tmp/12;
        display(['Oxygen file time     :        ',datestr(time_o2,'dd/mm/yyyy HH:MM:SS')]);
    else
        display(['Error decoding O2 time:       ',scale_units,'--> ','months']);
    end
    
    tmp=double(ncread([satfile,mois,'_01.nc'],'time'));
    scale_units = ncreadatt([satfile,mois,'_01.nc'],'time','units');
    if strfind(scale_units,'months')
        ref=sscanf(scale_units, 'months since %f-%f-%f %f:%f:%f');
        ref=datenum(ref');
        time_psat=ref+365.25*tmp/12;
        display(['Oxygen Saturation file time:  ',datestr(time_psat,'dd/mm/yyyy HH:MM:SS')]);
    else
        display(['Error decoding PSAT time:     ',scale_units,'--> ','months']);
    end  
    
    timewoa(i,1) = time_t;
    
    oxywoa(:,:,:,i) = ncread([oxyfile,mois,'_01.nc'],'o_an');
    oxywoa_sd(:,:,:,i) = ncread([oxyfile,mois,'_01.nc'],'o_sd');  
    psatwoa(:,:,:,i) = ncread([satfile,mois,'_01.nc'],'O_an');
    psatwoa_sd(:,:,:,i) = ncread([satfile,mois,'_01.nc'],'O_sd');  
    T_woa(:,:,:,i) = ncread([tfile,mois,'_01.nc'],'t_an'); 
    S_woa(:,:,:,i) = ncread([sfile,mois,'_01.nc'],'s_an');
    display(['<=====================>']);
end
[nLon,nLat,nZ_top,nT] = size(oxywoa);

% Variables communes
latwoa=double(ncread([oxyfile,'01_01.nc'],'lat'));
lonwoa=double(ncread([oxyfile,'01_01.nc'],'lon'));
Depth_m=ncread([sfile,'01_01.nc'],'depth');


CLIM_M.latwoa=latwoa;
CLIM_M.lonwoa=lonwoa;
CLIM_M.Depth_m=Depth_m;

CLIM_M.T_woa=T_woa;
CLIM_M.S_woa=S_woa;
CLIM_M.oxywoa=oxywoa;
CLIM_M.oxywoa_sd=oxywoa_sd;
CLIM_M.psatwoa=psatwoa;
CLIM_M.psatwoa_sd=psatwoa_sd;

CLIM_M.time_s=time_s;
CLIM_M.time_t=time_t;
CLIM_M.time_o2=time_o2;
CLIM_M.time_psat=time_psat;

CLIM_M.nLon=nLon;
CLIM_M.nLat=nLat;
CLIM_M.nZ_top=nZ_top;
CLIM_M.nT=nT;

end

function CLIM_A=read_annual_clim(DIR)

% Read the WOA 2018 Annual climatologies (0 to 5500m)

oxyfile=DIR.oxyfile;
tfile=DIR.tfile;
sfile=DIR.sfile;
satfile=DIR.satfile;

% Variables atlas WOA 2018 annuel(00)    360   180   102
mois='00';

filename=[oxyfile,mois,'_01.nc'];
display(['Reading Annual mean field: ',filename]);
Depth_a=double(ncread(filename,'depth'));

filename=[tfile,mois,'_01.nc'];
display(['Reading Annual mean field: ',filename]);
T_woa_a=ncread(filename,'t_an');

tmp=double(ncread(filename,'time'));
scale_units = ncreadatt(filename,'time','units');
if strfind(scale_units,'months')
    ref=sscanf(scale_units, 'months since %f-%f-%f %f:%f:%f');
    ref=datenum(ref');
    time_t=ref+365.25*tmp/12;
    display(['Temperature file time   : ',datestr(time_t,'dd/mm/yyyy HH:MM:SS')]);
else
    display(['Error decoding Temperature time: ',scale_units,'--> ','months']);
end
    
filename=[sfile,mois,'_01.nc'];
display(['Reading Annual mean field: ',filename]);
S_woa_a=ncread(filename,'s_an');

tmp=double(ncread(filename,'time'));
scale_units = ncreadatt(filename,'time','units');
if strfind(scale_units,'months')
    ref=sscanf(scale_units, 'months since %f-%f-%f %f:%f:%f');
    ref=datenum(ref');
    time_s=ref+365.25*tmp/12;
    display(['Salinity file time   : ',datestr(time_s,'dd/mm/yyyy HH:MM:SS')]);
else
    display(['Error decoding Salinity time: ',scale_units,'--> ','months']);
end

filename=[oxyfile,mois,'_01.nc'];
display(['Reading Annual mean field: ',filename]);
oxywoa_a=ncread(filename,'o_an');%mumole/kg
oxywoa_sd_a=ncread(filename,'o_sd');%mumole/kg
tmp=double(ncread(filename,'time'));
scale_units = ncreadatt(filename,'time','units');
if strfind(scale_units,'months')
    ref=sscanf(scale_units, 'months since %f-%f-%f %f:%f:%f');
    ref=datenum(ref');
    time_o2=ref+365.25*tmp/12;
    display(['Oxygen file time   : ',datestr(time_o2,'dd/mm/yyyy HH:MM:SS')]);
else
    display(['Error decoding Oxygen time: ',scale_units,'--> ','months']);
end

filename=[satfile,mois,'_01.nc'];
display(['Reading Annual mean field: ',filename]);
psatwoa_a=ncread(filename,'O_an');% percent
psatwoa_sd_a=ncread(filename,'O_sd');% percent
tmp=double(ncread(filename,'time'));
scale_units = ncreadatt(filename,'time','units');
if strfind(scale_units,'months')
    ref=sscanf(scale_units, 'months since %f-%f-%f %f:%f:%f');
    ref=datenum(ref');
    time_psat=ref+365.25*tmp/12;
    display(['Oxygen file time   : ',datestr(time_psat,'dd/mm/yyyy HH:MM:SS')]);
else
    display(['Error decoding Psat time: ',scale_units,'--> ','months']);
end
display(['<=====================>']);
[nLon,nLat,nZ] = size(oxywoa_a);

% Variables communes
latwoa=double(ncread([oxyfile,'00_01.nc'],'lat'));
lonwoa=double(ncread([oxyfile,'00_01.nc'],'lon'));

CLIM_A.latwoa=latwoa;
CLIM_A.lonwoa=lonwoa;
CLIM_A.Depth_a=Depth_a;

CLIM_A.T_woa_a=T_woa_a;
CLIM_A.S_woa_a=S_woa_a;
CLIM_A.oxywoa_a=oxywoa_a;
CLIM_A.oxywoa_sd_a=oxywoa_sd_a;
CLIM_A.psatwoa_a=psatwoa_a;
CLIM_A.psatwoa_sd_a=psatwoa_sd_a;

CLIM_A.time_s=time_s;
CLIM_A.time_t=time_t;
CLIM_A.time_o2=time_o2;
CLIM_A.time_psat=time_psat;

CLIM_A.nLon=nLon;
CLIM_A.nLat=nLat;
CLIM_A.nZ=nZ;

return
end

function DIR=define_dir

% Path, directories and file names are here define.

% Sea Water Library
SeaWaterLibrary='/Users/thierry_reynaud/IFREMER/MATLAB/LOCODOX/LOCODOX3.4/share/seawater/seawater_330_its90_lpo';
addpath(SeaWaterLibrary);

dirhome = '/Users/thierry_reynaud/IFREMER/MATLAB/LOCODOX/WOA/';
addpath(dirhome)
DIR.dirhome=dirhome;
% sauvegarde
dest =[dirhome,'WOA2018_LOCO/WOA2018_monthly_5500_1deg.nc'];
DIR.dest=dest;

% localisation fichiers woa mensuels et annuels
%dirwoa = [dirhome,'WOA2018_ORI/A5B7/'];%2005-2017
%pre='A5B7';%2005-2017

dirwoa = [dirhome,'WOA2018_ORI/DECAV/'];%1981-2010
pre='decav';%1981-2010

oxyfile = [dirwoa,'oxygen/woa18_all_o'];
satfile = [dirwoa,'o2sat/woa18_all_O'];
tfile = [dirwoa,'temperature/woa18_',pre,'_t'];
sfile = [dirwoa,'salinity/woa18_',pre,'_s'];

DIR.dirwoa=dirwoa;
DIR.pre=pre;
DIR.oxyfile=oxyfile;
DIR.satfile=satfile;
DIR.tfile=tfile;
DIR.sfile=sfile;
return
end
