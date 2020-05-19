function get_profiler(wmo)
% This function requieres 
% 1- Pulse Secure to be switched on or to be connected through the Ifremer internet network (IP with 134.246...)
% 2- An Ifremer account onto a unix work station
% 3- the installation of SSH Keys between the laptop and the Ifremer account onto a unix work stations
% Furthermore:
% Data are copied from: /home/oo26/coriolis/co05/co0508/dac/coriolis/
% Data are locally stored in: /Users/thierry_reynaud/IFREMER/MATLAB/LOCODOX/RAW_DATA_CORIOLIS/

%lineps='/Applications/Pulse\ Secure.app/Contents/MacOS/Pulse\ Secure';
%unix(lineps);


cwmo=num2str(wmo);
space=' ';
line1='rsync -avz -e ssh treynaud@pelican:/home/oo26/coriolis/co05/co0508/dac/coriolis/';
line2='/Users/thierry_reynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_FLOAT_DATA/';
line=[line1,cwmo,space,line2,'.'];
display(line);
unix(line);
return