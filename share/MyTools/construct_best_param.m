% construct_best_param replaces param by adjusted parameter if they exist
% for allparam
%
% SYNTAX
% [Ci]  =  construct_best_param(Co,allparam,Cp,verbose)
%
% DESCRIPTION
% construct_best_param replaces param by adjusted parameter if they exist
% for allparam. Keep param if not.
%
% INPUT
%     Co (structure)        float structures where to read param.
%                           Structures format: see NCR_file.m
%
%     allparam (cellstring) parameters used. Ex: {'temp','psal','pres'}
%
%     Cp (structure)        (optionnal) Output from a previous call of the
%                           function construct_best_param
%
%     verbose (boolean)     (optionnal, 0/1) Default 0. Make the program
%                           verbose.
%
% OUTPUT
%     Ci (structure)        float structure with the _adjusted parameter
%                           kept if exists.
%
% CALL : check_FirstDimArray
%
% SEE ALSO
%   NCR_file

% HISTORY
%   $created: //2012 $author: Cecile Cabanes, LPO, CNRS
%   $Revision: version $Date: $author:
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       adapted and corrected, to be shared to the O2 community


function [Ci]  =  construct_best_param(Co,allparam,Cp,verbose)

% =========================================================================
%% Initialisation
% =========================================================================
if nargin < 3
    Ci  =  [];
    verbose  =  1;
else
    Ci  =  Cp;
    if nargin < 4
        verbose  =  1;
    end
end

INITFIRSTDIM  =  [];
if isfield(Co,'firstdimname')
    INITFIRSTDIM  =  Co.firstdimname;
else
    INITFIRSTDIM  =  'N_HISTORY';
end
Co  =  check_FirstDimArray(Co,'N_PROF');
fillval  =  'FillValue_';

isADJprofile = true(Co.n_prof,1);

% =========================================================================
% Test if the profiles have adjusted values for all the parameter in
% allparam.
% Condition: 
% teste si les profils ont des valeurs ajustees pour tous les parametres allparam
% condition: il ne faut pas que les param_adjusted_qc soient  a FillValue
% si ce n'est pas la cas pour au moins un parametre alors isADJprofile = 0 pour le profil
% =========================================================================
for k = 1:length(allparam)
    param =  allparam{k};  
    if isfield(Co,param)        
        paramad = [param '_adjusted'];
        paramadqc = [paramad '_qc'];
        paramaderror = [paramad '_error'];
        
        if ~isfield(Co,paramad)
            isADJprofile = false(Co.n_prof,1); % test failed si il manque le parametre ajuste d'un des parametres allparam
            disp(['No adjusted parameter  for ' param])
        end
        
        [Test] = check_isfillval_prof(Co,paramadqc);

        if ~isempty(Test.(paramadqc))
            isADJprofile  =  isADJprofile & ~Test.(paramadqc);
        else
            disp(['No adjusted parameter QC for ' param])
            isADJprofile = false(Co.n_prof,1); % test fail si il manque le parametre ajusted_qc d'un des parametres allparam
        end
    else
        isADJprofile = false(Co.n_prof,1); % test fail si un parametre n'existe pas
    end
    
end

% verifie si c'est consistant avec le data mode

if isfield(Co,'data_mode')
    isADJfromDmode  =  Co.data_mode.data == 'D'|Co.data_mode.data == 'A';
end

if isequal(isADJprofile,isADJfromDmode) == 0
    if verbose == 1
        disp('******** Warning')
        disp(['Number of profiles with adjusted param well filled: ', num2str(sum(isADJprofile))])
        disp(['Number of profiles with data mode  =  A or D: ', num2str(sum(isADJfromDmode))])
        disp('Adjusted parameters may not be filled consistently with the data mode')
        disp('We require that both adjusted param are well filled and data mode = A or D')
        disp('********')
    end
end
disp(' ')

isADJprofile  =  isADJprofile&isADJfromDmode;

if sum(isADJprofile) ~= 0
    texte = [ 'In best: ' num2str(sum(isADJprofile)) ,' profiles are replaced by adjusted profiles for: '];
    if verbose == 1
        disp((texte))
        disp(allparam')
    end
    for k = 1:length(allparam)
        param =  allparam{k};
        if isfield(Co,param)
            paramqc = [param '_qc'];
            paramad = [param '_adjusted'];
            paramadqc = [paramad '_qc'];
            paramaderror = [paramad '_error'];
            
            paramnew = [param '_best'];
            paramnewqc = [param '_best_qc'];
            paramnewerror = [param '_best_error'];

            Ci.(paramnew) = Co.(paramad);
            Ci.(paramnewqc) = Co.(paramadqc);
            Ci.(paramnewerror) = Co.(paramaderror);
            Ci.(paramnew).name = upper(paramnew);
            Ci.(paramnewqc).name = upper(paramnewqc);
            Ci.(paramnewerror).name = upper(paramnewerror);
            % on prend les parametres TR si pas d'ajustement A ou D
            Ci.(paramnew).data(isADJprofile == 0,:)  =  Co.(param).data(isADJprofile == 0,:);
            Ci.(paramnewqc).data(isADJprofile == 0,:)  =  Co.(paramqc).data(isADJprofile == 0,:);
        end
    end
else
    for k = 1:length(allparam)
        param =  allparam{k};
        if isfield(Co,param)
            if verbose == 1
                disp(' In best: all profiles are RT profiles for :')
                disp(param)
            end
            paramqc = [param '_qc'];
            paramad = [param '_adjusted'];
            paramadqc = [paramad '_qc'];
            paramaderror = [paramad '_error'];
            
            paramnew = [param '_best'];
            paramnewqc = [param '_best_qc'];
            paramnewerror = [param '_best_error'];

            Ci.(paramnew) = Co.(param);
            Ci.(paramnewqc) = Co.(paramqc);
            Ci.(paramnew).name = upper(paramnew);
            Ci.(paramnewqc).name = upper(paramnewqc);            
        end        
    end
end

Ci.isADJprofile = isADJprofile;

%Co  =  check_FirstDimArray_is(Co,INITFIRSTDIM);
