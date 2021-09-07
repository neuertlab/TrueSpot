function [seg_files] = Get_Dir(Seg_dir,Pref_check)

%% The below code looks for files in the directory (Seg_dir) that start with the prefix Pref_check, and then returns a cell array (seg_files) with the complete path of all the files

% Seg_dir = 'C:\Users\keslerbk\Dropbox (VU Basic Sciences)\Vanderbilt Computer\2018-07-17 RNA Quant (20170607 Timecourse)';   %Directory to be searched
% Pref_check = 'mRNA';                     %What prefix to look for when storing file names
listing = dir(Seg_dir);           %Finds contents of directory
total_folders=length(listing);      %Total files in directory (plus 2)
clear seg_files
seg_files = {};                 %this will have the names of the seg files
counter = 1;
for i=3:total_folders  
    if size(listing(i).name,2) >= size(Pref_check,2)                %Check to see if content name is at least as long as the prefix (otherwise would error)
        if strcmp(listing(i).name(1:size(Pref_check,2)),Pref_check) %Check to see if it starts with the prefix
            seg_files{counter}  = listing(i).name;                  %Add to list of files
            counter = counter+1;
        end
    end
end

    
