%  check_profile_qc updates the profile_param_qc respect to the
%  param_adjusted_qc.
%
% SYNTAX
% [Co,GlobQC,ValnumQC] = check_profile_qc(Co,param,fillval)
%
% DESCRIPTION
%  check_profile_qc updates the profile_param_qc respect to the
%  param_adjusted_qc.
%
% INPUT
%   Co (structure)      Structure filled with all variables in a NetCDF
%                       Argofile. Made with NCR_file.m.
%                       Example:
%                       Co = 
%                             dc_reference: [1x1 struct]
%                                     juld: [1x1 struct]
%                                  juld_qc: [1x1 struct]
%                            juld_location: [1x1 struct]
%                                 latitude: [1x1 struct]
%                                longitude: [1x1 struct]
%                          profile_doxy_qc: [1x1 struct]
%                                     doxy: [1x1 struct]
%                                  doxy_qc: [1x1 struct]
%                            doxy_adjusted: [1x1 struct]
%                         doxy_adjusted_qc: [1x1 struct]
%                      doxy_adjusted_error: [1x1 struct]
%                                         ...
%                                         ...
%
%   param (cell)       (OPTIONNAL) name of the variable(s) for whom the
%                       profile_param_qc has to be updated. 
%                       By Default, take all the variables that have an
%                       _adjusted_qc field.
%                       Example : param = {'PSAL','TEMP'};
%                                 the programm will update profile_psal_qc
%                                 and profile_temp_qc by taking into
%                                 account psal_adjusted_qc and
%                                 temp_adjusted_qc.
%
%   fillval (scalar)    (OPTIONNAL) the value of the fillvalue of
%                        param_adjusted_qc. If not given, the program reads
%                        Co.(param).FillValue_.
%
% OUTPUT
%   Co (structure)      The input structure with updated profile_param_qc
%                       field.
%
%   GlobQC (array)      Contains the profile_param_qc updated. String
%                       array.
%                       Example: 
%                             GlobQC =
%                             B
%                             F
%
%   ValnumQC (double)   The proportion of "Good QC" all over the profile of
%                       param. GlobQC is the rating corresponding to this
%                       value.
%                       Example:
%                             ValnumQC =
%                                 0.9725
%                                      0
%
% CALL :
%   check_FirstDimArray, construct_best_param, format_flags_char2num
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: //2009 $author: Cecile Cabanes, IFREMER
%   $Revision: version $Date: $author:
%       v2 22/04/2016   Emilie Brion, ALTRAN OUEST
%                       adapted to be shared to the O2 community

function [Co,GlobQC,ValnumQC] = check_profile_qc(Co,param,fillval)

ValnumQC = [];
GlobQC = [];

if nargin<=2
    fillval = 'FillValue_';
end
if nargin == 1
    fillval = 'FillValue_';
    thefields = fieldnames(Co);
    ischarcell = strfind(thefields,'_adjusted_qc');
    ischar =~ cellfun('isempty',ischarcell);
    param = strrep(thefields(ischar),'_adjusted_qc','');
else
    param = lower(param);
    if ~iscell(param)
        param = {param};
    end
end

INITFIRSTDIM = [];
if isfield(Co,'firstdimname')
    INITFIRSTDIM = Co.firstdimname;
else
    INITFIRSTDIM = 'N_HISTORY';
end
Co = check_FirstDimArray(Co,'N_PROF');

for k = 1:length(param)
    onechamp = param{k};
    Co_best = construct_best_param(Co,{onechamp},[],0);
    Co_best = format_flags_char2num(Co_best);
    
    if isfield (Co_best, [onechamp '_best_qc'])        
        qc = Co_best.([onechamp '_best_qc']).data;        
        if ~isempty(Co_best.([onechamp '_best_qc']).(fillval))   
            isfill = (qc == Co_best.([onechamp '_best_qc']).(fillval)|qc == 9);
        else
            isfill = (qc == 9);
        end
        
        qc(qc == 1|qc == 2|qc == 5|qc == 8) = 1;
        qc(qc == 3|qc == 4|qc == 0) = 2;
        A = sum(qc == 1,2);
        B = sum(~isfill,2);
        nf = B ~= 0;
        ValnumQC = 999*ones(size(qc,1),1);
        GlobQC = repmat(' ',[size(qc,1),1]);
        
        ValnumQC(nf) = A(nf)./B(nf);
        
        GlobQC(ValnumQC == 0) = 'F';
        GlobQC(ValnumQC > 0) = 'E';
        GlobQC(ValnumQC >= 0.25) = 'D';
        GlobQC(ValnumQC >= 0.5) = 'C';
        GlobQC(ValnumQC >= 0.75) = 'B';
        GlobQC(ValnumQC == 1) = 'A';
        GlobQC(ValnumQC == 999) = ' ';
        
        Co.(['profile_' onechamp '_qc']).data = GlobQC;        
    end
end

Co = check_FirstDimArray(Co,INITFIRSTDIM);

