% extract_profile_dim extract one or more profiles or one or more levels
%
% SYNTAX
% [Coi,Dimi] = extract_profile_dim(Co,Dim,DIMNAME,iprofiles)
%
% DESCRIPTION
% extract_profile_dim extract one or more profiles or one or more levels
%
% INPUT
%     Co (structure)        float structures where to extract profiles.
%                           Structures format: see NCR_file.m
%
%     Dim (structure)       float Dimension structures. Structures
%                           format: see NCR_file.m
%
%     DIMNAME (string)      Dimension name (ex: 'N_PROF' or 'N_LEVELS')
%
%     iprofiles (double)    profiles  to extract ex:[1:3,5]
%
% OUTPUT
%     Co (structure)        float structure data concatenated from Coi and
%                           Cor, along the dimension DIMNAME
%
%     Dim (structure)       float structure Dimension for the structure Co
%
% CALL : check_FirstDimArray
%
% SEE ALSO
%   NCR_file

% HISTORY
%   $created: //2009 $author: Cecile Cabanes, LPO, CNRS
%   $Revision: version $Date: $author:
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       adapted and corrected, to be shared to the O2 community

function [Coi,Dimi] = extract_profile_dim(Co,Dim,DIMNAME,iprofiles)



if isempty(strfind(Co.obj,'ObsInSitu'))
    error('extract_profile not define for this type of structure')
else
    INITFIRSTDIM=[];
    if isfield(Co,'firstdimname')
        INITFIRSTDIM=Co.firstdimname;
    else
        INITFIRSTDIM='N_HISTORY';
    end
    Co = check_FirstDimArray(Co,DIMNAME);
    Coi=Co;
    Dimi=Dim;
    
    champs = fieldnames(Co);    %champs={'psal','psalqc','psalad',....}
    Nbfields = length(champs);
    
    for k=1:Nbfields            % boucle sur toutes les variables
        oneChamp=champs{k};
        if isfield(Co.(oneChamp),'data')
            if isempty(Co.(oneChamp).data)==0
                isthedim=strcmp(Co.(oneChamp).dim,DIMNAME);
                if sum(isthedim)==1
                    nbdim = length(size(Coi.(oneChamp).data));
                    ap='';
                    if nbdim>1
                        ap=repmat(',:',[1,nbdim-1]);
                    end
                    expre=['Coi.(oneChamp).data=Co.(oneChamp).data(iprofiles' ap ');'];
                    eval(expre);
                    %Coi.(oneChamp).data=Co.(oneChamp).data(iprofiles,:);
                end
            end
        end
    end
    dimnamei=fieldnames(Dimi);
    for k=1:length(dimnamei)
        if strcmpi(dimnamei{k},DIMNAME) % dimension selon laquelle les champs sont extraits
            Dimi.(dimnamei{k}).dimlength = sum(iprofiles);
        end
    end
    Coi = check_FirstDimArray(Coi,DIMNAME);
    if isempty(INITFIRSTDIM)==0
        Coi = check_FirstDimArray(Coi,INITFIRSTDIM);
    end
end
