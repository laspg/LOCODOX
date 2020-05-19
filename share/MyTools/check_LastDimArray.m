% check_LastDimArray checks that the last dimension of all the array in the
% structure Co are the same as DIMNAME
%
% SYNTAX
% [Co] = check_LastDimArray(Co,DIMNAME)
%
% DESCRIPTION
% check_LastDimArray check that the last dimension of all the array in the
% structure Co are the same as DIMNANE and change dimension order if
% necessary
%
% INPUT
%   Co (structure)      Float structure
%   DIMNAME (char)      Fimension name. Ex: 'N_PROF'
%
% OUTPUT
%   Co (structure)      Float structure with the DIMNAME dimension in the
%                       last place for all the arrays. Required field
%                       for each variable : lastdimname (double)
%                       Example : 
%                           Co.lastdimname = 'N_PROF'
%
% CALL :
%
% SEE ALSO
%   NCR_file
%

% HISTORY
%   $created: //2009 $author: Cecile Cabanes, LPO, CNRS
%   $Revision: version $Date: $author:
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       adapted to be shared to the O2 community

function [Co] = check_LastDimArray(Co,DIMNAME)


% =========================================================================
%% Initialisation
% =========================================================================

lowname=lower(DIMNAME);


gotoend=0;
if isfield(Co,'firstdimname')
    if strcmp(Co.firstdimname,DIMNAME)
        gotoend = 1;
    end
end

% =========================================================================
%% Check for the last dimension for all fields.
% =========================================================================
myfields = fieldnames(Co);    %myfields={'psal','psalqc','psalad',....}
Nbfields = length(myfields);
testonechange = 0;
Co.(lowname) = 0;

% Loop over all the fields (variables)
for k=1:Nbfields
    oneChamp = myfields{k};
    if isfield(Co.(oneChamp),'data')
        if isfield(Co.(oneChamp),'dim')
            isthedim=strcmp(Co.(oneChamp).dim,DIMNAME);
%             if sum(isthedim)==1
            if any(isthedim)
                testonechange = 1;
                if gotoend==0
                    vecdim = 1:length(isthedim); 
                    vecdim_sauv = vecdim;
                    vecdim = circshift(vecdim,[1,vecdim(isthedim)-1]);
                    Co.(oneChamp).dim = circshift(Co.(oneChamp).dim,[1,vecdim_sauv(isthedim)-1]);
                    
                    if length(vecdim)>1
                        Co.(oneChamp).data = permute(Co.(oneChamp).data,vecdim);
                    else
                        if size(Co.(oneChamp).data,2)>1
                            Co.(oneChamp).data = Co.(oneChamp).data';
                        end
                    end
                    
                    % Save the last dimension
                    if ~isempty(Co.(oneChamp).data)
                        Co.(lowname) = size(Co.(oneChamp).data,1);
                    end
                else
                    lowname = lower(DIMNAME);
                    if ~isempty(Co.(oneChamp).data)
                        Co.(lowname) = size(Co.(oneChamp).data,1);
                    end
                end
            end
        end
    end
end

% =========================================================================
%% Change the dimension order if DIMNAME is not the last dimension
% =========================================================================
if ~testonechange
    disp(['Does not find this dimension: ',DIMNAME])    
else
    Co.lastdimname = DIMNAME;
    if isfield(Co,'firstdimname')
       Co = rmfield(Co,'firstdimname');
    end
end