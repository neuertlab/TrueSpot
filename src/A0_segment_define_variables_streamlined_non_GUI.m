%% This code is to define the variables needed to run segmentation or RNA threshold determination.
function A0_segment_define_variables_streamlined_non_GUI(node_num) %A0_segment_define_variables_streamlined_Jurkat

Diffstack = 0;
Seg = 1;                                                                    %Set to 1 if you want to threshold the images (must be done before RNA)
Yth = 0;                                                                    %set to 1 if just want thresholding and not segmentation
RNA_thres = 0;                                                              %Set to 1 if you want to determine RNA spots (must be done after segmentation)
max_int_thres = 1;                                                          %Set to 1 if you want to use the maximum intensity projection for RNA thresholding visualization
max_int_spot_th = 1;                                                        % set to 1 if you want the spot detection to be carried out on the maximum intensity projection during thresholding 
Ych = [0 1 0 1 0 1 0 0 1];                                                  % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
yeast_seg = 0;                                                                             % This is needed to compute the number of images                                                                                                                                                 
Ywin = 0;                                                                  %specifies if using windows in order to use "\" or "/". Also runs in parallel if == 0.
Yim = 0;                                                                   % Show images of nucleus and trans segmentation
S = 0;                                                                     % S = 1 save files; S = 0 do not save files; 
file_type = 'multiple' ;                                                    % 'single' or ' multiple'   % Defines if images have been saved as single files or if all files have been saved in a single file
%Only need to change if different organism. Below is for MEFs
%% Jurkat 100x
min_nucleus_size = 7000;                                                   %MEF 10000
max_nucleus_size = 35000;                                                  %MEF 200000
min_cell_size = 15000;                                                       % 250; % minimum cell size
max_cell_size = 60000;                                                      % 1300;% maximum cell size

% %% Jurkat 20x
% min_nucleus_size = 300;                                                  
% max_nucleus_size = 4000; 
%                                         
% min_cell_size = 400;
% max_cell_size = 5000;   

dxy1 = round(sqrt(max_nucleus_size/pi));                                       %radius of square for nucleus determination
img_stacks = [45,75];                                                        %specify the image stacks used for segmentation. These are number of stacks above the plane with least contrast (in focus)

%% Directory for FISH images. Choose between local directory [ A ] or a directory on your cluster account [ B ]
if Ywin == 0;;
   fish_dir = strcat('/gpfs23/scratch/neuertg/BenK/2019_Diff_Timecourses/'); %[ B ]
%     fish_dir = strcat('/Users/gregorneuert/Dropbox (VU Basic Sciences)/Neuert lab/Benjamin Kesler/Segmentation GUI Generalized/Sample Images for GUI/Jurkat/100x/DAPI and TRANS in one stack (per image)');
    % Jurcat 100x
else
  %fish_dir = strcat('C:\Users\keslerbk\Dropbox (VU Basic Sciences)\Segmentation Paper\Segmentation GUI Generalized\Sample Images for GUI\Jurkat\100x\DAPI and TRANS in one stack (per image)'); %[A]
  %  fish_dir = strcat('D:\Dropbox (VU Basic Sciences)\Segmentation Paper\Segmentation GUI Generalized\Sample Images for GUI\Jurkat\100x\DAPI and TRANS in one stack (per image)'); %[A]
  fish_dir = 'D:\Dropbox (VU Basic Sciences)\Segmentation Paper\Folder for Github Upload\Sample Images for GUI\Jurkat\100x\DAPI and TRANS in one stack (per image)'     
%   fish_dir = strcat('F:\Dropbox (VU Basic Sciences)\Neuert lab\Segmentation Paper\Segmentation GUI Generalized\Sample Images for GUI\Jurkat\100x');
end

%% Directory to store segment information. choose between local directory [ A ] or a directory on your cluster account [ B ] for output of seg files 
% or where they are currently stored
if Ywin == 0;
    
%    outfile_prefix_seg=strcat('/gpfs23/scratch/neuertg/BenK/20180517/'); % [ B ]
    outfile_prefix_seg = strcat('/gpfs23/scratch/neuertg/BenK/2019_Diff_Timecourses/');
else
   outfile_prefix_seg=strcat('D:\Dropbox (VU Basic Sciences)\Ben Kesler\Segmentation GUI Generalized\Segmentation Output Files\'); %[ A ]
%    outfile_prefix_seg = strcat('F:\Dropbox (VU Basic Sciences)\Neuert lab\Segmentation Paper\Segmentation GUI Generalized\Individual Codes\Segmented Cells_fixed threshold\Jurkat_100x');
end
 
%% RNA threshold specific variables %%%%%%%%%%%%%%%%%%%%%%
% choose between local directory [ A ] or a directory on your cluster account [ B ] for output of RNA files
% or where they are currently stored
if Ywin == 0;;
    outfile_prefix_RNA=strcat('/gpfs23/scratch/neuertg/BenK/'); % [ B ]
else
   outfile_prefix_RNA=strcat('C:\Users\keslerbk\Documents\2017-10-19 RNA Thresholding (20171010)\'); %[ A ]
end

thA = [1,1,200]; % Thresholds for RNA spot intensities
ths  = [50 50 21 50 101 50 50 50 50]; % Thresholds for RNA spot intensities % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS; defines which fluorecnet dyes have been imaged; 

%% Segmentation specific variables
%%
Seg_thres = 0;  % 0                                                            %Set to 1 if want to segment using DAPI threshold instead of the adding method. Needed for manual thresholding
ManTh = 0;      %0                                                            %set to 1 if you want to manually choose a threshold (only works if Seg_thres == 1)                                                                 
     
% Chooses type of segmentation process
segment_mode = 'last5';                                                     %either "midplane" or "midplane2" or "max_cells" or "first5" or "last5"

[images1,folders1] = Get_Dir_All(fish_dir,Ywin)                                           %Finds the file names within the fish directory
if Ywin == 0
    images1 = {images1{node_num+1}}
    folders1 = {folders1{node_num+1}}
end
clear images2 im_prefixes  
 images2 = {};          %This will have the file names (not the entire directory and not the type of file)
 for i = 1:size(images1,2)
     if Ywin                     %The type of slash depends on whether it is windows
     beginning1 = strfind(images1{i},'\');      
     else
         beginning1 = strfind(images1{i},'/');
     end
     beginning1 = beginning1(end)+1;    %Find the part that is not directory
     if strfind(images1{i},'_MMStack')
         end1 = strfind(images1{i},'_MMStack'); %look for MMStack since it is often at the end of the file names
     else
         end1 = strfind(images1{i},'.');        %Look for a period otherwise as the end of the file name
     end
     end1 = end1(1)-1;
     images2{i} = images1{i}(beginning1:end1);  %Populate the cell array with the file names
 end 
 images2';  %For visualizing all files extracted (remove semicolon)
  handles.images1 = images1;
  handles.images2 = images2;
if Diffstack
    clear handles.images_DAPI handles.images_TRANS
    handles.images_DAPI = {};
    handles.images_DAPI = {};
  counter_DAPI = 1;             %counter for marker images found +1 (used to populate both DAPI and TRANS cell arrays)
  counter_TRANS = 0;            %counter for TRANS images found (used to determine if TRANS prefix was correct)
  handles.DAPI_prefix
  for i = 1:size(images2,2) 
      if strfind(images2{i},handles.DAPI_prefix) == 1
          handles.images_DAPI{counter_DAPI} = images1{i};
          im_prefixes{counter_DAPI} =  images2{i}(size(handles.DAPI_prefix,2)+1:end);   %populate imprefixes (used for later segmentation file names)
          for j = 1:size(images2,2)
              if strfind(images2{j},handles.TRANS_prefix) == 1 &...         %Look for TRANS prefix at beginning
                      strcmp(images2{j}(size(handles.TRANS_prefix,2)+1:end),... %Look for the rest of image name to match between TRANS
                      images2{i}(size(handles.DAPI_prefix,2)+1:end))        %and DAPI
                      handles.images_TRANS{counter_DAPI} = images1{j}; 
              end
          end
          counter_TRANS = counter_TRANS+1;                                  %Add to counter (even if match not found)
          counter_DAPI = counter_DAPI+1;                                    %Add to counter (even if match not found)
       if strfind(images2{i},handles.TRANS_prefix) %searching for TRANS prefix
%            handles.images_TRANS{counter_TRANS} = images1{i};
           counter_TRANS = counter_TRANS+1;
       end
      end
  end
  if counter_TRANS == 0
      msgbox('There were no images found with the boundary prefix')
  end
  if isempty(handles.images_DAPI)
      msgbox('There were no images found with the marker prefix')     
  elseif isempty(handles.images_TRANS)
      msgbox('There were no images paired with the prefixes specified')
  end
else
%[images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                %Finds the file names within the fish directory (old way)
im_prefixes = images2;
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
end
%% Old way of finding file names
% 
% %% Experiment date and note (need to be adjusted)
 % exp_date = {'20180619','20180620','20180621','20180622','20180623','20180626','20180628','20180629'};                                                        % Date of the experiment and folder name on the hard disk
% strain = 'F1-2-1';
%task_ids = [1:5];                                                          %Set to how many timepoints there are.
% osmo = '0.2MStep';
% genepair = {'Xist-CY5-Maoa-TMR','Xist-CY5-Pdk3-TMR','Xist-CY5-RepA-TMR','Xist-CY5-Tsix-TMR'};
% file_end = '_MMStack.ome.tif';
% 
% % Time points
% exp_names = {'0d' '6hr' '12hr','1d','2d','5d','9d','9d'};% '6min' '8min' '10min' '12min' '14min',...
%    % '15min' '18min' '20min' '22min' '25min' '40min' '50min' '120min' '122min' '129min'};
% 
% positions_touse = { ...
%  [1:11];... %  0min 4                                                        %Change each number to the range of image numbers each timepoint has
%  [1:11];... %  1min 4
%  [1:11];... %  2min 4
%   [1:11];... %  4min 4 
%  % [1:11];... %  6min 6
%  % [1:6];... %  8min 4
%  % [1:9];... %  10min
%  % [1:6];... %  12min
%  % [1:4];... %  14min
%  % [1:5];... %  16min
%  % [1:7];... %  18min
%  % [1:6];... %  20min 4
% %  [1:8];... %  22min
% %  [1:7];... %  25min 4
% %  %[1:4];... %  30min 4
% %  %[1:5];... %  35min 4
% %  [1:8];... %  40min 4
% %  %[1:4];... %  45min 4
% %  [1:6];... %  50min 4
% %  [1:8];... %  120min 4
% %  [1:5];... %  122min 4
% %  [1:4];... %  129min
%  };

% %% run this if you want to know how many total images there are. Set ntasks
% %to this number and array to 0-(1-number of images) in the cluster file
% if Ywin == 0 
%     Yim = 0;                                                                   % Cannot show images on cluster
%     img_pos = zeros(2,1);                                                       %An array for specifying which node will process which image.
%                                                                             %the first row indicates the timepoint and the second row is the image number
%                                                                             %The total number of columns should equal the number of images processed
%     %if size(positions_touse,1) == size(exp_names,2)
%         counter = 1;
%         for i = 1:size(exp_date,2)
%             for k = 1: size(genepair,2)
%             for j = 1:size(positions_touse{k},2)
%                 img_pos(1,counter) = i;
%                 img_pos(2,counter) = j;
%                 img_pos(3,counter) = k;
%                 counter = counter+1;
%             end
%             end
%         end
% %     else
% %         'Number of timepoints (exp_names) specified does not equal size of positions_touse'
% %         a = notavariable
% %     end
%     img_pos
% 
%     task_ids = img_pos(1,node_num+1);
%     img_num = img_pos(2,node_num+1);
%     pairnums = img_pos(3,node_num+1);
%     total_images = size(img_pos,2)
% else
%     img_num = 1;                                                            %if in windows, placeholder for image number until it is redefined in later functions
%     node_num = 1;
%     pairnums = 1:size(genepair,2);
% end




%% Returns first file name of microscopy image for comparison
%  exp_name = [exp_date '_' osmo '_' strain '_' exp_names{1}  '_' genepair '_img'];        
%         if Ywin == 0;;
%             filename = [fish_dir '/' exp_name '_1' '/' exp_name '_1' file_end];
%         else
%             filename = [fish_dir '\' exp_name '_1' '\' exp_name '_1' file_end];
%         end;

%% Run the segmentation or RNA threshold determination
 if Seg == 1
     A1_segment_predefined_variables_streamlined_nonGUI(Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_stacks,yeast_seg,images1,im_prefixes,dxy1)
 end
 if RNA_thres == 1
     Run_RNA_detection_predefined_variables(thA,ths,outfile_prefix_RNA,Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,img_num,node_num,max_int_thres,max_int_spot_th)
 end   
