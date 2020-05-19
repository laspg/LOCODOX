% format_flags_char2num changes flag char strings to numerical vectors
%
% SYNTAX
% Ci  =  format_flags_char2num(Co)
%
% DESCRIPTION
% format_flags_char2num changes flag char strings to numerical vectors,
% in "_QC" fields.
% Example:'11111441111 ' -> [1 1 1 1 1 4 4 1 1 1 1 999]
%
% INPUT
%     Co (structure)        float structure for which to replace char
%                           strings to numerical vectors.
%                           Structures format: see NCR_file.m
%
% OUTPUT
%     Co (structure)        float structure with "_QC" fields for which the
%                           char string have been replaced by numerical
%                           vectors.
%
% CALL
%
% SEE ALSO
%   NCR_file

% HISTORY
%   $created: //2009 $author: Cecile Cabanes, LPO, CNRS
%   $Revision: version $Date: $author:
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       adapted and corrected, to be shared to the O2 community

function Ci  =  format_flags_char2num(Co)


fillval = 'FillValue_';
champs  =  fieldnames(Co);    %champs = {'psal','psalqc','psalad',....}
Nbfields  =  length(champs);
Ci = Co;

for k = 1:Nbfields            % boucle sur toutes les variables
    oneChamp = champs{k};
    if isempty(findstr(oneChamp,'_qc')) == 0            % test si il y a 'qc' dans le nom du champ
        if isempty(Ci.(oneChamp).data) == 0            % test si le champ est rempli
            Ci.(oneChamp).ischar2num = 0;  % variable logique qui identifie
            % si le tableau de flag a ete transforme char-> num ( = 1) ou non ( = 0)
            if ischar(Ci.(oneChamp).type) == 1
                if strcmp(Ci.(oneChamp).type, 'char') == 1
                    thisiscar = 1;
                else
                    thisiscar = 0;
                end
            else
                if Ci.(oneChamp).type == 2
                    thisiscar = 1;
                else
                    thisiscar = 0;
                end
            end
            
            if thisiscar == 1      % test si c'est une chaine de caractère
                if length(size(Ci.(oneChamp).data))>2       % test si c'est un tableau de 2 dim
                   % warning('Does not accept array of dim >2')
                else
                    % test si le tableau de flag est alphanumerique charactere ie '0' '1' '2' ...'9' seulement
                    % (les fillvalues ne sont pas prisent en compte, elles peuvent ne pas etre alphanumeriques ex ' ')
                    poubflag  =  Ci.(oneChamp).data;  % tableau temporaire
                    poubflag(poubflag == Ci.(oneChamp).(fillval)) = '0'; % On remplace les fillvalue par un caratere alphanumerique ('0') pour le test seulement.
                    poubflag(isstrprop(poubflag,'cntrl')) = '0';
                    isalphanum = numel(isstrprop(poubflag,'digit')) == numel(poubflag); % = 1 si tableau alphanumerique
                    clear poubflag

                    if isalphanum == 1  % test si tableau alphanumerique
                        % on transforme le tableau de flag en tableau numerique
                        Ci.(oneChamp) = rmfield(Ci.(oneChamp),'data');
                        Ci.(oneChamp).data = single(999*ones(size(Co.(oneChamp).data))); % fillval numeriques
                        Ci.(oneChamp).type = 5; % single precision
                        Ci.(oneChamp).(fillval) = single(999);
                        size1 = size(Ci.(oneChamp).data,1);
                        %keyboard
                        %size1
                        %oneChamp
                        for i = 1:size1                              % boucle sur chaque profil                            
                            flagstr = Co.(oneChamp).data(i,:)' ;
                            isnofill = (flagstr ~= Co.(oneChamp).(fillval)&~isstrprop(flagstr,'cntrl'))';
                            flagnum = str2num(flagstr(isnofill))';
                            if sum(isnofill) ~= 0
                                Ci.(oneChamp).data(i,isnofill) = single(flagnum);
                            end
                        end
                        Ci.(oneChamp).ischar2num = 1;
                        Co.(oneChamp).ischar2num = 1;
                        Ci.(oneChamp)  =  orderfields(Ci.(oneChamp),Co.(oneChamp)); % ordonne les champs comme dans la structure initiale
                    end
                end
            end
        end
    end
end
