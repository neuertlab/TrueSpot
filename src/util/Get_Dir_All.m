function [images1,folders1] = Get_Dir_All(Search_dir,Ywin)

%% The below code uses dir to list the contents of a directory and then populates 
% the cell arrays folders1 with the contents of the directory, then images1
% with the contents of the directories as long as they start with the
% directory name
% This is how stack images are usually stored
% Search_dir = 'C:\Users\keslerbk\Desktop\Microscopy\20180619';   %Directory to be searched
temp_fold = {Search_dir};
clear folders1 images1
folders1 = {Search_dir};      %folders will be listed here. Different rows are different layers inside, and columns are the numbers within a specific layer
images1 = {};        %contents of directories in Search_dir will be listed here (usually images)
counter1 = 1;       %counter for number of folders needed to get through
counter_dir = 1; %counter for folders
counter_file = 1; %counter for files
folder_layer = 1;   %counter for the folder layer 
while counter1
    counter_dir_temp = 1;   
    clear temp_fold1
    temp_fold1 = {};         %Folders to loop through for this specific layer of folders
    for j = 1:size(temp_fold,2)
        temp_fold{j};
        listing = dir(temp_fold{j});           %Finds contents of directory
        total_content=length(listing);      %Total files in directory (plus 2)
        for i=3:total_content                               %loop through contents of search directory (first two entries are not actually content)
            listing(i).name;
            if listing(i).isdir
                counter1 = counter1+1;          %Keep track that another folder was found
                if Ywin
                    folders1{folder_layer+1,counter_dir_temp}=strcat(temp_fold{j},'\',listing(i).name);    %Store entire path of directory
                    temp_fold1{counter_dir_temp} = strcat(temp_fold{j},'\',listing(i).name);              %Store entire path of directory
                else
                    folders1{folder_layer,counter_dir_temp}=strcat(temp_fold{j},'/',listing(i).name);    %Store entire path of directory
                    temp_fold1{counter_dir_temp} = strcat(temp_fold{j},'/',listing(i).name);    %Store entire path of directory
                end
                counter_dir_temp = counter_dir_temp+1;
                counter_dir = counter_dir+1;
            else
                if Ywin
                    images1{counter_file}=strcat(temp_fold{j},'\',listing(i).name);    %gives one content of the search directory to folders
                else
                    images1{counter_file}=strcat(temp_fold{j},'/',listing(i).name);    %gives one content of the search directory to folders
                end
                counter_file = counter_file+1;
            end
        end
        counter1 = counter1-1;
    end
    temp_fold = temp_fold1;
    folder_layer = folder_layer+1;
end
%[indx,tf] = listdlg('ListString',images1,'ListSize',[1000 500],'PromptString','Select which images to process') % creates a modal dialog box that allows the user to select one or more items from the specified list.