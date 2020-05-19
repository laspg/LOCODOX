% str2numQC convert string to numbers, with NaN instead of blancks
%
% SYNTAX
% [QC_num] = str2numQC(QC_str)
%
% DESCRIPTION
% str2numQC convert string to numbers, with NaN instead of blancks
% (usefull for Argo QC)
%
% INPUT
%     QC_str (string)       dimension (char): nc x nz 
%
% OUTPUT
%     QC_num (double)       same fields of QC_str but in double format
%                           dimension (char): nc x nz
%
% SEE ALSO
%   str2num
%
% HISTORY
%   $created: 
%   $Revision: version $Date: $author:
%

function[QC_num] = str2numQC(QC_str)

    QC_num = NaN(size(QC_str,1),size(QC_str,2));
    for ip = 1:size(QC_num,1)
        for im = 1:size(QC_num,2)
            if ~strcmp(QC_str(ip,im),' ') 
                QC_num(ip,im) = str2num(QC_str(ip,im));
            end
        end
    end
    
end