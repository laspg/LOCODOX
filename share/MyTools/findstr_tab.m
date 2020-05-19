% findstr_tab finds a string in a string array or a cell array
%
% SYNTAX
% [isfound] = findstr_tab(tab,pattern)
%
% DESCRIPTION
% findstr_tab finds a string in a string array or a cell array
%
% INPUT
%     tab (string array)  String array or cell array of string
%         (cell array)  
%
%     pattern (string)        String to look for in tab
%
% OUTPUT
%     isfound (double)    indicates the starting index of each occurence of
%                         the pattern in tab.
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

function [isfound] = findstr_tab(tab,pattern)


if ~iscell(tab)
    % make a cell string from tab
    tabcell = cellstr(tab);
else
    tabcell = tab;
end

% Find the pattern in tab
ischarcell = strfind(tabcell,pattern );

% Look for ischarcell not empty
isfound =~ cellfun('isempty',ischarcell);   
