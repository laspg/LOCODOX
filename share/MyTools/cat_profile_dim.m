% cat_profile_dim cat two float structure along one dimension
%
% SYNTAX
% [Co,Dim] = cat_profile_dim(Coi,Cor,Dimi,Dimr,DIMNAME)
%
% DESCRIPTION
% cat_profile_dim cat two float structure along one dimension
%
% INPUT
%     Coi, Cor (structure)    float structures to be concatenated along the
%                             dimension DIMNAME. Structures format: see
%                             NCR_file.m
%
%     Dimi, Dimr (structure)  float Dimension structures. Structures
%                             format: see NCR_file.m
%
%     DIMNAME (string)        Dimension along which the structures should be
%                             concatenated
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

function [Co,Dim] = cat_profile_dim(Coi,Cor,Dimi,Dimr,DIMNAME)

% if isempty(strfind(Coi.obj,'ObsInSitu'))
%     error('replace_profile not define for this type of structure')
% else

INITFIRSTDIM = [];
if isfield(Coi,'firstdimname')
    INITFIRSTDIM = Coi.firstdimname;
end
%      keyboard
Cor = check_FirstDimArray(Cor,DIMNAME);
if ~isempty(Coi)
    Coi = check_FirstDimArray(Coi,DIMNAME);
end

% keyboard
if ~isempty(Coi)
    champsi = fieldnames(Coi);    %champs = {'psal','psalqc','psalad',....}
else
    champsi = fieldnames(Cor);
end
champsr = fieldnames(Cor);


% =========================================================================
%% Empty fields management in the structures to be concatenated:
% The empty fields are filled by FillValues
% =========================================================================

% -------------------------------------------------------------------------
% Empty fields in Coi  => fillvalues
% -------------------------------------------------------------------------
if ~isempty(Coi)
    Nbfields = length(champsi);
    for k = 1:Nbfields
        oneChamp = champsi{k};
        if isfield(Coi.(oneChamp),'data')
            if isempty(Coi.(oneChamp).data)
                clear siz
                for l = 1:length(Coi.(oneChamp).dim)
                    siz(l) = Dimi.(lower(Coi.(oneChamp).dim{l})).dimlength;
                    if siz(l) == 0
                        siz(l) = 1;
                    end
                end
                siz(length(Coi.(oneChamp).dim)+1) = 1;
                if isfield(Coi.(oneChamp),'FillValue_')
                    Coi.(oneChamp).data = repmat(Coi.(oneChamp).FillValue_,siz);
                else
                    disp(oneChamp);
                    error('Did not find FillValue for this field')
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% Empty fields in Cor  => fillvalues
% -------------------------------------------------------------------------
Nbfields = length(champsr);
for k = 1:Nbfields
    oneChamp = champsr{k};
    if isfield(Cor.(oneChamp),'data')
        if isempty(Cor.(oneChamp).data)
            clear siz
            for l = 1:length(Cor.(oneChamp).dim)
                siz(l) = Dimr.(lower(Cor.(oneChamp).dim{l})).dimlength;
                if siz(l) == 0
                    siz(l) = 1;
                end
            end
            siz(length(Cor.(oneChamp).dim)+1) = 1;
            if isfield(Cor.(oneChamp),'FillValue_')
                Cor.(oneChamp).data = repmat(Cor.(oneChamp).FillValue_,siz);
            else
                disp(oneChamp);
                error('Did not find FillValue for this field')
            end
        end
    end
end



% =========================================================================
%% Unik values management in the structures:
% For a field present in only one structure, the field is created in the
% other structure and filled with fillValues.
% =========================================================================

% -------------------------------------------------------------------------
% Fields only in Coi  => field created in Cor and filled with fillvalues
% -------------------------------------------------------------------------
champs_uni = setdiff(champsi,champsr);
Nbfields = length(champs_uni);
for k = 1:Nbfields
    oneChamp = champs_uni{k};
    if isfield(Coi.(oneChamp),'data')
        Cor.(oneChamp) = Coi.(oneChamp);
        % Find the dimension
        clear siz
        for l = 1:length(Cor.(oneChamp).dim)
            if isfield(Dimr,lower(Cor.(oneChamp).dim{l}))
                siz(l) = Dimr.(lower(Cor.(oneChamp).dim{l})).dimlength;
                if siz(l) == 0
                    siz(l) = 1;
                end
            else
                siz(l) = Dimi.(lower(Cor.(oneChamp).dim{l})).dimlength;
                %Dimr.(lower(Cor.(oneChamp))) = Dimi.(lower(Cor.(oneChamp)));% cree la dimension pour Cor si elle n'existe pas
                Dimr.(lower(Cor.(oneChamp).dim{l})) = Dimi.(lower(Cor.(oneChamp).dim{l})); %--MG correction de bug
                
                
                if siz(l) == 0
                    siz(l) = 1;
                end
            end
        end
        siz(length(Coi.(oneChamp).dim)+1) = 1;
        
        % Create new array in Cor
        if isfield(Cor.(oneChamp),'FillValue_')
            Cor.(oneChamp).data = repmat(Cor.(oneChamp).FillValue_,siz);
        else
            disp(oneChamp);
            error('Did not find FillValue for this field')
        end
    end
end

% -------------------------------------------------------------------------
% Fields only in Cor  => field created in Coi and filled with fillvalues
% -------------------------------------------------------------------------
champs_unr = setdiff(champsr,champsi);
Nbfields = length(champs_unr);
for k = 1:Nbfields
    oneChamp = champs_unr{k};
    if isfield(Cor.(oneChamp),'data')
        Coi.(oneChamp) = Cor.(oneChamp);
        % Find the dimension
        clear siz
        for l = 1:length(Coi.(oneChamp).dim)
            if isfield(Dimi,lower(Coi.(oneChamp).dim{l}))
                siz(l) = Dimi.(lower(Coi.(oneChamp).dim{l})).dimlength;
                %disp('toto')
                if siz(l) == 0
                    siz(l) = 1;
                end
            else
                siz(l) = Dimr.(lower(Coi.(oneChamp).dim{l})).dimlength;
                Dimi.(lower(Coi.(oneChamp).dim{l})) = Dimr.(lower(Coi.(oneChamp).dim{l})); % cree la dimension pour Cor si elle n'existe pas
                if siz(l) == 0
                    siz(l) = 1;
                end
            end
        end
        siz(length(Coi.(oneChamp).dim)+1) = 1;
        
        % Create new array in Coi
        if isfield(Coi.(oneChamp),'FillValue_')
            %siz
            %size(Coi.(oneChamp).FillValue_)
            Coi.(oneChamp).data = repmat(Coi.(oneChamp).FillValue_,siz);
            %size(Coi.(oneChamp).data)
        else
            disp(oneChamp);
            error('Did not find FillValue for this field')
        end
    end
end

% =========================================================================
%% Common fields management
% =========================================================================
if ~isempty(Coi)
    champsi = fieldnames(Coi);    %champs = {'psal','psalqc','psalad',....}
else
    champsi = fieldnames(Cor);
end
%     champsi = fieldnames(Coi);    %champs = {'psal','psalqc','psalad',....}
champsr = fieldnames(Cor);
champs = union(champsi,champsr);
Nbfields = length(champs);
if ~isempty(Coi)
    Co = Coi;
else
    Co = Cor;
end

% Loop over all the variables
for k = 1:Nbfields
    oneChamp = champs{k};
    if isfield(Coi,oneChamp) & isfield(Cor,oneChamp)
        if isfield(Coi.(oneChamp),'data') & isfield(Cor.(oneChamp),'data')
            if ~isempty(Coi.(oneChamp).data) & ~isempty(Cor.(oneChamp).data)
                
                if sum(strcmp(Coi.(oneChamp).dim ,Cor.(oneChamp).dim)) == length(Coi.(oneChamp).dim)  % fields with same dim names
                    isthedim_i = strcmp(Coi.(oneChamp).dim,DIMNAME);
                    isthedim_r = strcmp(Cor.(oneChamp).dim,DIMNAME);
                    
                    if sum(isthedim_i) == 1 & sum(isthedim_r) == 1 % if DIMNAME exists in both fields
                        otherdim_i = (isthedim_i == 0);          % look to the other dimension for the i field
                        otherdim_r = (isthedim_r == 0);          % look to the other dimension for the r field
                        fullsizei = size(Coi.(oneChamp).data);
                        if numel(fullsizei) == numel(isthedim_i) - 1
                            fullsizei(end+1) = 1;
                        end
                        sizei = fullsizei(otherdim_i);
                        if isempty(sizei);sizei = 1;end;
                        fullsizer = size(Cor.(oneChamp).data);
                        if numel(fullsizer) == numel(isthedim_r) - 1
                            fullsizer(end+1) = 1;
                        end
                        sizer = fullsizer(otherdim_r);
                        if isempty(sizer);sizer = 1;end;
                        maxsize = max(sizei,sizer);
                        nbdim_i = length(sizei);
                        nbdim_r = length(sizer);
                        
                        if isequal(sizei,sizer) == 1
                            ap = '';
                            for m = 1:nbdim_i
                                ap = [ap ',1:sizei(' num2str(m) ')'];
                            end
                            %themax = [fullsizei(1),sizei];
                            if isfield(Coi.(oneChamp),'FillValue_')
                                %Co.(oneChamp).data = repmat(Coi.(oneChamp).FillValue_,themax);
                                expre = ['Co.(oneChamp).data(:' ap ') = Coi.(oneChamp).data;'];
                                %expre
                                eval(expre)
                            end
                        end
                        if isequal(sizei,maxsize) == 0
                            ap = '';
                            for m = 1:nbdim_i
                                ap = [ap ',1:sizei(' num2str(m) ')'];
                            end
                            themax = [fullsizei(1),maxsize];
                            if isfield(Coi.(oneChamp),'FillValue_')
                                Co.(oneChamp).data = repmat(Coi.(oneChamp).FillValue_,themax);
                                expre = ['Co.(oneChamp).data(:' ap ') = Coi.(oneChamp).data;'];
                                %expre
                                eval(expre)
                            end
                        end
                        
                        if isequal(sizer,maxsize) == 0
                            ap = '';
                            for m = 1:nbdim_r
                                ap = [ap, ',1:sizer(' num2str(m) ')'];
                            end
                            themax = [fullsizer(1),maxsize];
                            tempo = Cor.(oneChamp).data;
                            if isfield(Cor.(oneChamp),'FillValue_')
                                Cor.(oneChamp).data = repmat(Cor.(oneChamp).FillValue_,themax);
                                expre = ['Cor.(oneChamp).data(:' ap ') = tempo;'];
                                %expre
                                eval(expre)
                                clear tempo
                            end
                        end
                        
                        nbdim = length(size(Co.(oneChamp).data));
                        iprofiles = size(Co.(oneChamp).data,1)+1;
                        rprofiles = size(Co.(oneChamp).data,1)+size(Cor.(oneChamp).data,1);
                        api = '';
                        apr = '';
                        if nbdim>1
                            api = repmat(',:',[1,nbdim-1]);
                            apr = repmat(',:',[1,nbdim-1]);
                        end
                        expre = ['Co.(oneChamp).data(' num2str(iprofiles) ':' num2str(rprofiles)  api ') = Cor.(oneChamp).data(:' apr ');'];
                        %expre
                        eval(expre)
                        
                    end
                end
            end
        end
    end
end
if ~isempty(Dimi)
    dimnamei = fieldnames(Dimi);
else
    dimnamei = [];
end
dimnamer = fieldnames(Dimr);
if isempty(Dimi)
    Dimi = Dimr;
end
Dim = Dimi;
for k = 1:length(dimnamer)
    
    if sum(strcmp(dimnamei,dimnamer{k})) == 1
        if strcmp(dimnamer{k},lower(DIMNAME))
            Dim.(dimnamer{k}).dimlength = Dimi.(dimnamer{k}).dimlength+Dimr.(dimnamer{k}).dimlength;
        else
            Dim.(dimnamer{k}).dimlength = max(Dimi.(dimnamer{k}).dimlength,Dimr.(dimnamer{k}).dimlength);
            
        end
    else
        Dim.(dimnamer{k}) = Dimr.(dimnamer{k});
    end
end
Co = check_FirstDimArray(Co,DIMNAME);
if isempty(INITFIRSTDIM) == 0
    Co = check_FirstDimArray(Co,INITFIRSTDIM);
end