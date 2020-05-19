% =========================================================================
% Create reference database (structure REF) for O2 correction with LOCODOX
% Create the structure .mat used as input in the software  
%
% adaptation by M.GALLIAN (ALTRAN) in february 2019 of the function gen_bdd_O2ref_all.m 
% (Virginie Thierry and Emilie Brion (ALTRAN)) 
% 
% Explaination :
% The user should change the "path_file" and set the corresponding function to decode the netcdf file ("name_func"). 
% Run this code as many netcdf files you have to add to the database, and select "extend" after the first run.
% Be carefull to not write 2 times the same file in the reference matrix, it will induce an error in LOCODOX

% =========================================================================

 clear all;
 close all;

% =========================================================================
%% PARAMETERS
% =========================================================================
% data type : 
% - Hydro LOPS : name_func=rd_O2data_hydro_LOPS
% - Strass : name_func=rd_O2data_hydro_strass
% - Ovide2011me : name_func= rd_O2data_hydro_ovid11
% - Ovide2011di : name_func= rd_O2data_hydro_ovid11

%Define the function to use to decode the netcdf file and its directory
%The function has to return : refId,lon,lat,juld,pres,temp,psal,sig0,doxy
%Path of the netcdf file to decode or name (=input 1 of the function)

path_func='Function/';

for i=1:9

    if i==1; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/ovid10_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
    if i==2; path_file='//Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/cat12_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
    if i==3; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/geov_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
    if i==4; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/bo16_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
    if i==5; path_file='/Volumes/qlpo5/HYDROCEAN/OVIDE18_TEMP/MLT/ov18_d_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
    if i==7; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/RREX/RREX15/rr15_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
    if i==8; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/RREX/RREX17/rr17_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
    if i==9; path_file='strass';name_func='rd_O2data_hydro_strass';end
    if i==10; path_file='ov11me';name_func= 'rd_O2data_hydro_ovid11';end
    if i==11; path_file='ov11di';name_func= 'rd_O2data_hydro_ovid11';end


%Path of the .mat structure to update, or create
BDD_O2_REF= 'bddo2ref_avecov18_temp.mat';


% =========================================================================
%% Initialisation
% =========================================================================

%Add path of the function
addpath(path_func)

%Initialise the REF structure
REF.id   =   [];
presi = 0:10:4000;
REF.lon  =   [];
REF.lat  =   [];
REF.juld =   [];
REF.temp =   [];
REF.psal =   [];
REF.sig0 =   [];
REF.doxy =   [];
REF.pres =   [];
REF.psat =   [];

% =========================================================================
%% Load reference matrix or overwrite
% =========================================================================

if exist(BDD_O2_REF,'file')    
    answer = questdlg('Reference database exists ! Do you want to overwrite or extend it ?','Refrence database','Overwrite','Extend','Extend');
    if strcmp(answer,'Overwrite')
        eval(['!\rm ' BDD_O2_REF]);
        REF.id   =   [];
        presi = 0:10:4000;
        REF.lon  =   [];
        REF.lat  =   [];
        REF.juld =   [];
        REF.temp =   [];
        REF.psal =   [];
        REF.sig0 =   [];
        REF.doxy =   [];
        REF.pres =   [];
        REF.psat =   [];
    else
        load(BDD_O2_REF) ;
    end
end


% =========================================================================
%% Create reference database
% =========================================================================

%Extract data
[refId,lon,lat,~,juld,pres,temp,~,psal,sig0,doxy,~] = feval(name_func,path_file,'ctd');

%Interpolate data on pressure
if ~isempty(lon) && ~strcmp(name_func,'rd_O2data_hydro_LOPS')
    for i = 1:length(lon)
        iok = find(isfinite(pres(i,:)));
        tempi(i,:) = interp1(pres(i,iok),temp(i,iok),presi);
        psali(i,:) = interp1(pres(i,iok),psal(i,iok),presi);
        sig0i(i,:) = interp1(pres(i,iok),sig0(i,iok),presi);
        doxyi(i,:) = interp1(pres(i,iok),doxy(i,iok),presi);
    end
elseif ~isempty(lon) && strcmp(name_func,'rd_O2data_hydro_LOPS')
    tempi=temp;
    psali=psal;
    sig0i=sig0;
    doxyi=doxy;
end

%Add reference data to the REF structure
if ~isempty(lon)
    REF.lon = [REF.lon;lon];
    REF.lat = [REF.lat;lat];
    REF.juld = [REF.juld;juld];
    REF.temp = [REF.temp;tempi];
    REF.psal = [REF.psal;psali];
    REF.sig0 = [REF.sig0;sig0i];
    REF.doxy = [REF.doxy;doxyi];
    REF.id = [REF.id;refId];
    REF.pres = ones(length(REF.lon),1)*presi;
    tmpsat = sw_satO2(REF.psal,REF.temp);
    sat = DOXY_convert(tmpsat,'mL/L','mumol/kg',REF.sig0);
    REF.psat = 100*REF.doxy./sat;
    
    % Save the data
    save(BDD_O2_REF,'REF');
end


    end
