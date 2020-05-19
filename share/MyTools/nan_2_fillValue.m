% nan_2_fillValue replace fillvalue by NaN
%
% SYNTAX
% [Co] = nan_2_fillValue(Co,fillvalName)
%
% DESCRIPTION
% nan_2_fillValue replace fillValue by NaN if the variable is numeric
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
% OUTPUT
%     Co (structure)        float structure with the fillvalue field
%                           required for each variable.
%                           Example :
%                               Co.temp.FillValue_ = 99999
%
% CALL :
%
% SEE ALSO
%   NCR_file

% HISTORY
%   $created: //2009 $author: Cecile Cabanes, LPO, CNRS
%   $Revision: version $Date: $author:
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       adapted to be shared to the O2 community

function [Co] = nan_2_fillValue(Co,fillvalName)

% =========================================================================
%% Initialisation
% =========================================================================

if nargin==1
    fillvalName='FillValue_';
end

% =========================================================================
%% Check for the fillnan
% =========================================================================
myfields = fieldnames(Co);    %myfields={'psal','psalqc','psalad',....}
Nbfields = length(myfields);
gotoend = 0;

% Loop over all the fields (variables)
if gotoend==0
    for k=1:Nbfields
        oneChamp = myfields{k};
        if isfield(Co.(oneChamp),'data')
            if ~isempty(Co.(oneChamp).data)
                if isnumeric(Co.(oneChamp).data)
                    if isfield(Co.(oneChamp),fillvalName)
                        if ~isempty(Co.(oneChamp).(fillvalName))
                            selec_fill = isnan(Co.(oneChamp).data);
                            Co.(oneChamp).data(selec_fill) = Co.(oneChamp).(fillvalName);
                            Co.fillisnan = 0;
                        else
                            warning([oneChamp ': Did not find a fillvalue attribute'])
                        end
                    end
                end
            end
        end
    end
end
