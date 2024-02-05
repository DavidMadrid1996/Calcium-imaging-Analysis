function [ValidPathdataCalcium, ValidPathdatawaves,Animal_List_unique] = File_Location(Animal_List,format)

% For excel files --> '\**\*.xlsx'
% For matlab files --> '\**\*.mat'
clear ValidPathdataCalcium ValidPathdatawaves Animal_List_unique

mainPath = "E:\EXPERIMENTAL\CALCIUM_IMAGING\DAVID\WHOLE_BRAIN_PROJECT_SWIMMING\VALID_RECONSTRUCTION_Final_Script_to_score_cells";
files = dir(strcat(mainPath,format));
NumberFiles = size(char(files.name),1);
NameFiles = char(files.name);
Path2Files = char(files.folder);
Animal_List_unique = Animal_List;

for i = 1:size(Animal_List ,1)
    temp = strcat(Animal_List(:,1),Animal_List(:,2));

    for ii = 1:size(Path2Files,1)
        clear out out1 out2
        out=regexp(Path2Files(ii,:),'\','split');  % Separating in folders and subfolders
        %         out1 = regexprep(out(end), '\W',''); % Removing empty spaces
        out1 =regexprep(out(end), '[^A-Za-z0-9_.]', ''); % Name of the trial and animal
        out2 = regexprep(out(end-1), '[^A-Za-z0-9_.]', '');  % Folder dataCalcium or data wave


        if isequal(out1,temp(i)) && isequal(out2,{'dataCALCIUM'})
            ValidPathdataCalcium{i,:} = strrep(Path2Files(ii,:), ' ','');
        end
        if isequal(out1,temp(i)) && isequal(out2,{'datawaves'})
            ValidPathdatawaves{i,:} = strrep(Path2Files(ii,:), ' ','');
        end
    end
    Animal_List_unique{i,4} = ValidPathdataCalcium{i,:};
    Animal_List_unique{i,5} = ValidPathdatawaves{i,:};
end


end






