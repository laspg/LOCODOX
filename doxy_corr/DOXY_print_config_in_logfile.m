% DOXY_print_config_in_logfile
%
% SYNTAX
% [] = DOXY_print_config_in_logfile(log_summ,CONFIG,wmo)
%
% DESCRIPTION
% DOXY_print_config_in_logfile print the parameters defined in the config file in the log file
%
% INPUT
%       log_summ (double)     File identifier 
%
%       CONFIG (struct)       Configuration structure with data path,
%                             operator choices, ...
%                             CONFIG =
%                           DataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/data_coriolis/INAIR/'
%                       NCEPDataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/Data_NCEP'
%                       saveDataDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/data/'
%                       savePlotDir: '/home6/pharos/ebrion/LOCODOX/SOFT_2017_09_04/plots/'
%                                ...
%                                ...
%                       inAirFormat: 0
%                           inAirMC: 1090
%                                ...
%
%       wmo (double)           WMO of the float
%
% OUTPUT
%
% CALL :
%
% SEE ALSO
% DOXY_corr_main
%
% HISTORY
%   $created: 01/04/2019 $author:  Marine GALLIAN, Altran Ouest 
%   $Revision: version $Date: $author:


function [] = DOXY_print_config_in_logfile(log_summ,CONFIG,wmo)

%% Heading
fprintf(log_summ,'================================================================================\n');
fprintf(log_summ,'\n');
fprintf(log_summ,'               SUMMARY OF THE OPTIONS APPLIED FOR THE CORRECTION\n');
fprintf(log_summ,'\n');
fprintf(log_summ,sprintf('                             Float %d                 \n',wmo));
fprintf(log_summ,'================================================================================\n');
fprintf(log_summ,'\n');
fprintf(log_summ,'\n');
fprintf(log_summ,'\n');

%% Configuration file
fprintf(log_summ,'================================================================================\n');
fprintf(log_summ,'---------------------    PART 1 : CONFIGURATION FILE     -----------------------\n');
fprintf(log_summ,'================================================================================\n');
fprintf(log_summ,'\n');

%% Directories and files
fprintf(log_summ,'---------------->    DIRECTORIES AND FILES       \n');
fprintf(log_summ,'\n');

%Directories
fprintf(log_summ,'Directories : \n');
fprintf(log_summ,sprintf('- Main directory of LOCODOX : %s\n',CONFIG.LocodoxMainDir));
fprintf(log_summ,sprintf('- Directory of the Argo NetCDF data : %s\n',CONFIG.DataDir));
fprintf(log_summ,sprintf('- Directory of the NCEP data : %s\n',CONFIG.NCEPDataDir));
fprintf(log_summ,sprintf('- Directory for saving : %s\n',CONFIG.resultsDir));
fprintf(log_summ,'\n');

%Files
fprintf(log_summ,'Files : \n');
fprintf(log_summ,sprintf('- Climatology file : %s\n',CONFIG.WOAfile));
fprintf(log_summ,sprintf('- Reference In-Situ Database : %s\n',CONFIG.bddFile));
fprintf(log_summ,sprintf('- RefArgoFile, gathering the wmo of the floats and the associated in-situ reference profiles : %s\n',CONFIG.RefArgoFile));
fprintf(log_summ,sprintf('- Mask file : %s\n',CONFIG.maskFile));
fprintf(log_summ,'\n');

%NCEP informations
fprintf(log_summ,'NCEP files and FTP address : \n');
if CONFIG.ncepDoUpdate==1
    fprintf(log_summ,sprintf('- NCEP data has been updated \n'));
else
    fprintf(log_summ,sprintf('- NCEP data hasn''t been updated \n'));
end
fprintf(log_summ,sprintf('- Ftp website of NCEP : %s\n',CONFIG.ncepFtp));
fprintf(log_summ,sprintf('- Path where to find the NCEP data in the ftp website : %s\n',CONFIG.ncepFtpDir));
fprintf(log_summ,sprintf('- Sub directory where to find the NCEP data in the ftp website : %s\n',strjoin(CONFIG.ncepFtpSubDir)));
fprintf(log_summ,sprintf('- NCEP files to be read : %s\n',strjoin(CONFIG.ncepFiles)));
fprintf(log_summ,sprintf('- NCEP data is read for years : %s\n',num2str(CONFIG.ncepYears)));
fprintf(log_summ,'\n');
fprintf(log_summ,'\n');


%%  Parameter for correction
fprintf(log_summ,'---------------->    PARAMETERS FOR CORRECTION  \n');
fprintf(log_summ,'\n');

fprintf(log_summ,sprintf('- Reference data unit : %s\n',CONFIG.refUnit));
if CONFIG.presEff==0
    fprintf(log_summ,sprintf('- Pressure effect correction for conversion DOXY/PSAT/PPOX : unactivated (presEff = 0)\n'));
else
    fprintf(log_summ,sprintf('- Pressure effect correction for conversion DOXY/PSAT/PPOX : activated (prefEff = 1)\n'));
end
if CONFIG.isokC==0
    fprintf(log_summ,sprintf('- Carry over parameter : In-Air oxygen measurements are supposed not to be biaised by splash of water (isokC = 0)\n'));
else
    fprintf(log_summ,sprintf('- Carry over parameter : In-Air oxygen measurements are supposed to be biaised by splash of water (isokC = 1)\n'));
end
fprintf(log_summ,'\n');

%Data mode selection 
fprintf(log_summ,sprintf('Data mode selection :\n'));
if CONFIG.DM_pres==1
    fprintf(log_summ,sprintf('- "Delayed Mode" field is privileged for Pressure (--> PRES_ADJUSTED))\n'));
else
    fprintf(log_summ,sprintf('- "Real Time" field is used for Pressure (--> PRES))\n'));
end
if CONFIG.DM_temp==1
    fprintf(log_summ,sprintf('- "Delayed Mode" field is privileged for Temperature (--> TEMP_ADJUSTED))\n'));
else
    fprintf(log_summ,sprintf('- "Real Time" field is used for Temperature (--> TEMP))\n'));
end
if CONFIG.DM_psal==1
    fprintf(log_summ,sprintf('- "Delayed Mode" field is privileged for Salinity (--> PSAL_ADJUSTED))\n'));
else
    fprintf(log_summ,sprintf('- "Real Time" field is used for Salinity (--> PSAL))\n'));
end
fprintf(log_summ,'\n');

%Qc selection
fprintf(log_summ,sprintf('QC selection :\n'));
fprintf(log_summ,sprintf('- Only the DOXY data whose QC equal to [%s] are retained for the calculation of the correction.\n',num2str(CONFIG.QC_O)));
fprintf(log_summ,sprintf('- Only the PRES data whose QC equal to [%s] are retained for the calculation of the correction.\n',num2str(CONFIG.QC_P)));
fprintf(log_summ,sprintf('- Only the TEMP data whose QC equal to [%s] are retained for the calculation of the correction.\n',num2str(CONFIG.QC_T)));
fprintf(log_summ,sprintf('- Only the PSAL data whose QC equal to [%s] are retained for the calculation of the correction.\n',num2str(CONFIG.QC_S)));
fprintf(log_summ,'\n');

%Drift parameter
fprintf(log_summ,'Drift parameter : \n');
if CONFIG.drift_spec==0
    fprintf(log_summ,sprintf('- The data drift is computed using the polynomial fitting of order : %s \n',num2str(CONFIG.drift_fitPolynomialDegree)));
else
    fprintf(log_summ,sprintf('- The data drift is computed using the following equation : %s \n',formula(CONFIG.drift_fittype)));
end
fprintf(log_summ,sprintf('- Minimum depth necessary for calculating the drift : %s m \n', num2str(CONFIG.min_drift_depth)));
fprintf(log_summ,'\n');

%WOA or REF
fprintf(log_summ,sprintf('For WOA of REF correction :\n'));
fprintf(log_summ,sprintf('- The linear correction R2 coefficient (threshold under which LOCODOX suggests to apply constant correction) is equal to : %s \n',num2str(CONFIG.R2threshold)));
fprintf(log_summ,'\n');

%IN AIR 
fprintf(log_summ,sprintf('For IN AIR correction :\n'));
fprintf(log_summ,sprintf('- Measurement code for In-Air samples : %s \n',num2str(CONFIG.inAirMC)));
fprintf(log_summ,sprintf('- Measurement code for In-Water samples : %s \n',num2str(CONFIG.inWaterMC)));
fprintf(log_summ,sprintf('- The maximum depth for TEMP and PSAL selection as "surface temperature" and "surface salinity" is : %s m\n',num2str(CONFIG.inAirMaxPresForTS)));
fprintf(log_summ,'\n');

%For netcdf and history 
fprintf(log_summ,'For netcdf writing and for history : \n');
fprintf(log_summ,sprintf('- Software name : %s \n',CONFIG.history_software)); 
fprintf(log_summ,sprintf('- Software version used to correct the data : %s\n', CONFIG.history_software_release));
fprintf(log_summ,sprintf('- Reference data set : %s \n', CONFIG.history_reference));


%blank
fprintf(log_summ,'\n');
fprintf(log_summ,'\n');
fprintf(log_summ,'\n');
























end
