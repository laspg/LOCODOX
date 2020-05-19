% fillValue_2_nan replace NaN value by fill
%
% SYNTAX
% [Co] = fillValue_2_nan(Co,fillvalName)
%
% DESCRIPTION
% fillValue_2_nan replace NaN value by fill if the variable is numeric
% (usefull before creating a file )
%
% INPUT
%     Co (structure)        float structure with the fillvalue field
%                           required for each variable.
%                           Example :
%                               Co.temp.FillValue_ = 99999
%
%     fillvalName (string)  (OPTIONNAL) The name of the fillvalue,
%                           corresponding to a filed of the structure Co
%
%     verbose (boolean)     (OPTIONNAL) ([0/1]) makes the program verbose
%
% OUTPUT
%     Co (structure)        float structure with the fillvalue field
%                           required for each variable.
%                           Example :
%                               Co.temp.FillValue_ = 99999
%
% CALL :
%
% SEE ALSO
%   create_netcdf_allthefile

% HISTORY
%   $created: //2009 $author: Cecile Cabanes, LPO, CNRS
%   $Revision: version $Date: $author:
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       adapted to be shared to the O2 community

function [Co] = fillValue_2_nan(Co,fillvalName,verbose)

% =========================================================================
%% Initialisation
% =========================================================================

if nargin==1
    fillvalName = 'FillValue_';
    verbose = 0;
elseif nargin == 2
    verbose = 0;
end

% =========================================================================
%% Check for the fillnan
% =========================================================================
myfields = fieldnames(Co);    %myfields={'psal','psalqc','psalad',....}
Nbfields = length(myfields);
gotoend = 0;
% if isfield(Co, 'fillisnan')
%     if Co.fillisnan==0
%         gotoend = 1;
%     end
% end

% Loop over all the fields (variables)
if gotoend==0
    for k=1:Nbfields
        oneChamp = myfields{k};
        if isfield(Co.(oneChamp),'data')
            if ~isempty(Co.(oneChamp).data)
                if isnumeric(Co.(oneChamp).data)
                    if isfield(Co.(oneChamp),fillvalName)
                    %if ~isempty(Co.(oneChamp).(fillvalName))         
                        %selec_fill = isnan(Co.(oneChamp).data);
                        selec_fill = Co.(oneChamp).data == Co.(oneChamp).(fillvalName);
                        if Co.(oneChamp).type == 4 || Co.(oneChamp).type == 2
                            % type = int32 or string : NaN could not exist : keep
                            % fillValue as is                            
                        else
                            Co.(oneChamp).data(selec_fill) = NaN;
                            % Co.(oneChamp).data(selec_fill) = Co.(oneChamp).(fillvalName);
                        end
                        Co.fillisnan = 1;
                    else
                        if verbose
                            warning([oneChamp ': Did not find a fillvalue attribute'])
                        end
                    end
                end
            end
        end
    end
end
