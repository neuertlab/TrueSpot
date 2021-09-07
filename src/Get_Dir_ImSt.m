function [images1,prefixes2] = Get_Dir_ImSt(Search_dir,Ywin)

%% The below code uses dir to list the contents of a directory and then populates 
% the cell arrays folders1 with the contents of the directory, then images1
% with the contents of the directories as long as they start with the
% directory name
% This is how stack images are usually stored
% Search_dir = 'C:\Users\keslerbk\Desktop\Microscopy\20180619';   %Directory to be searched
listing = dir(Search_dir);           %Finds contents of directory
total_folders=length(listing);      %Total files in directory (plus 2)
clear folders1 images1 prefixes1 prefixes2
folders1 = {};      %contents of Search_dir will be listed here (image folders)
images1 = {};        %contents of directories in Search_dir will be listed here (usually images)
prefixes1 = {};        %prefixes for each subdirectory
prefixes2 = {};        %prefixes for each image
counter = 1;        %counter for images1
counter1 = 1;       %counter for folders1 (not used at the moment) 
for i=3:total_folders                               %loop through contents of search directory (first two entries are not actually content)
    if Ywin
        folders1{i-2}=strcat(Search_dir,'\',listing(i).name);    %gives one content of the search directory to folders
    else
        folders1{i-2}=strcat(Search_dir,'/',listing(i).name);    %gives one content of the search directory to folders
    end
    prefixes1{i-2} = listing(i).name;
    d1=dir(folders1{i-2});                           %Finds contents of the content or subdirectory found in previous line
    total_files2=length(d1);                         %length of contents in subdirectory
    for j=3:length(d1)   %loop through contents of subdirectory (first two entries are not actually content)                   
        if size(listing(i).name,2) < size(d1(j).name,2)     %Check if the name of the subdirectory is shorter than its contents (true for stack images)
            if strcmp(d1(j).name(1:size(listing(i).name,2)-1),listing(i).name(1:size(listing(i).name,2)-1))    % Only add files in the folder that start with the same name as the folder (done with stack images)
                if Ywin
                    images1{counter}=strcat(folders1{i-2},'\',d1(j).name);   %Add the entire path for each image
                else
                    images1{counter}=strcat(folders1{i-2},'/',d1(j).name);   %Add the entire path for each image
                end
                prefixes2{counter} = prefixes1{i-2};
                counter = counter+1;
            end
        end
    end
end