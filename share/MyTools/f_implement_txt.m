function file2=f_implement_txt(file1,svg2)

%%% =====================================================================
%%% file2=f_implement_txt(file1,svg2) >> function to update txt file
%%% INPUT: 
%%%     file1 = INPUT file to update
%%%     svg2 =  string to update in file1
%%% OUTPUT:
%%%     file2 = file1 + svg2
%%% vr, August 2018
%%% ======================================================================

% open file 1
if ~(exist(file1,'file'))
    fid = fopen(file1,'w');
    fclose(fid)
end
fid = fopen(file1,'r');

% get lines from file 1 
ligne = fgetl(fid); % L1
count = 1; 
while ischar(ligne)
    Tligne{count} = ligne;
    ligne = fgetl(fid); % following line
    count = count + 1;
end

% close file1
fclose(fid);

% create new file2 with name of file1
file2 = file1;
fid2 = fopen(file2,'w');
% save all lines from file 1
if exist ('Tligne')
    for i = 1:length(Tligne)
        line2 = Tligne{i};
        fprintf(fid2,'%s\n',line2);
    end
end

% save new line
fprintf(fid2,'%s\n',svg2);

%close file2
fclose(fid2);

