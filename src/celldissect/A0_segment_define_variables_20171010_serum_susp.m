%% This code is to define the variables needed to run segmentation or RNA threshold determination.
function A0_segment_define_variables_no_probes(node_num)

Seg = 1;                                                                    %Set to 1 if you want to threshold the images (must be done before RNA)
Yth = 0;                                                                     %set to 1 if just want thresholding and not segmentation
RNA_thres = 0;                                                              %Set to 1 if you want to determine RNA spots (must be done after segmentation)
max_int_thres = 1;                                                          %Set to 1 if you want to use the maximum intensity projection for RNA thresholding visualization
max_int_spot_th = 1;                                                        % set to 1 if you want the spot detection to be carried out on the maximum intensity projection during thresholding 
Ych = [0 1 1 1 1 0 0 0 1];                                                  % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
                                                                             % This is needed to compute the number of images                                                                                                                                                 
Ywin = 1;                                                                  %specifies if using windows in order to use "\" or "/". Also runs in parallel if == 0.
Yim = 0;                                                                   % Show images of nucleus and trans segmentation
task_ids = [1:2];                                                          %Set to how many timepoints there are.

%Only need to change if different organism. Below is for MEFs
min_nucleus_size = 4000;                                                   %MEF 10000
max_nucleus_size = 20000;                                                  %MEF 200000

min_cell_size = 8000;                                                       % 250; % minimum cell size
max_cell_size = 50000;                                                      % 1300;% maximum cell size
img_stacks = [45,60];                                                        %specify the image stacks used for segmentation. These are number of stacks above the plane with least contrast (in focus)


%% Experiment date and note (need to be adjusted)
file_type = 'multiple' ;                                                    % 'single' or ' multiple'   % Defines if images have been saved as single files or if all files have been saved in a single file
exp_date = '20171010';                                                        % Date of the experiment and folder name on the hard disk
strain = 'mESC_2i-serum-susp';
osmo = '0.2MStep';
genepair = 'Xist-CY5-Tsix-AF594-Jpx-TMR';
file_end = '_MMStack.ome.tif';

% Time points
exp_names = {'4d' '4dRA' '2dRA'};% '6min' '8min' '10min' '12min' '14min',...
   % '15min' '18min' '20min' '22min' '25min' '40min' '50min' '120min' '122min' '129min'};

positions_touse = { ...
 [1:6];... %  0min 4                                                        %Change each number to the range of image numbers each timepoint has
 [1:6];... %  1min 4
 [1:10];... %  2min 4
%  [1:4];... %  4min 4 
%  [1:4];... %  6min 6
%  [1:6];... %  8min 4
%  [1:9];... %  10min
%  [1:6];... %  12min
%  [1:4];... %  14min
%  [1:5];... %  16min
%  [1:7];... %  18min
%  [1:6];... %  20min 4
%  [1:8];... %  22min
%  [1:7];... %  25min 4
%  %[1:4];... %  30min 4
%  %[1:5];... %  35min 4
%  [1:8];... %  40min 4
%  %[1:4];... %  45min 4
%  [1:6];... %  50min 4
%  [1:8];... %  120min 4
%  [1:5];... %  122min 4
%  [1:4];... %  129min
 };

%run this if you want to know how many total images there are. Set ntasks
%to this number and array to 0-(1-number of images) in the cluster file
if Ywin == 0 
    Yim = 0;                                                                   % Cannot show images on cluster
    img_pos = zeros(2,1);                                                       %An array for specifying which node will process which image.
                                                                            %the first row indicates the timepoint and the second row is the image number
                                                                            %The total number of columns should equal the number of images processed
    if size(positions_touse,1) == size(exp_names,2)
        counter = 1;
        for i = 1:size(positions_touse,1)
            for j = 1:size(positions_touse{i},2)
                img_pos(1,counter) = i;
                img_pos(2,counter) = j;
                counter = counter+1;
            end
        end
    else
        'Number of timepoints (exp_names) specified does not equal size of positions_touse'
        a = notavariable
    end
    img_pos

    task_ids = img_pos(1,node_num+1);
    img_num = img_pos(2,node_num+1);
    total_images = size(img_pos,2)
else
    img_num = 1;                                                            %if in windows, placeholder for image number until it is redefined in later functions
    node_num = 1;
end


%% Directory for FISH images. Choose between local directory [ A ] or a directory on your cluster account [ B ]
if Ywin == 0;;
    fish_dir = strcat('/gpfs23/scratch/neuertg/BenK/',exp_date); %[ B ]
else
    fish_dir = strcat('C:\Users\keslerbk\Desktop\Microscopy\',exp_date); %[A]
end

%% Directory to store segment information. choose between local directory [ A ] or a directory on your cluster account [ B ] for output of seg files 
% or where they are currently stored
if Ywin == 0;;
    outfile_prefix_seg=strcat('/gpfs23/scratch/neuertg/BenK/',exp_date, '/'); % [ B ]
else
   outfile_prefix_seg=strcat('C:\Users\keslerbk\Documents\2017-10-04 mESC Segmentation (20171003 20171010)\'); %[ A ]
end

%% Returns first file name of microscopy image for comparison
%  exp_name = [exp_date '_' osmo '_' strain '_' exp_names{1}  '_' genepair '_img'];        
%         if Ywin == 0;;
%             filename = [fish_dir '/' exp_name '_1' '/' exp_name '_1' file_end];
%         else
%             filename = [fish_dir '\' exp_name '_1' '\' exp_name '_1' file_end];
%         end;

%% Segmentation specific variables
%%
Seg_thres = 0;                                                              %Set to 1 if want to segment using DAPI threshold instead of the adding method. Needed for manual thresholding
ManTh = 0;                                                                  %set to 1 if you want to manually choose a threshold (only works if Seg_thres == 1)                                                                 
     
% Chooses type of segmentation process
segment_mode = 'last5';                                                     %either "midplane" or "midplane2" or "max_cells" or "first5" or "last5"

%% RNA threshold specific variables %%%%%%%%%%%%%%%%%%%%%%
%% choose between local directory [ A ] or a directory on your cluster account [ B ] for output of RNA files
% or where they are currently stored
if Ywin == 0;;
    outfile_prefix_RNA=strcat('/gpfs23/scratch/neuertg/BenK/',exp_date, '/'); % [ B ]
else
   outfile_prefix_RNA=strcat('C:\Users\keslerbk\Documents\2017-10-19 RNA Thresholding (20171010)\'); %[ A ]
end

thA = [1,1,200]; % Thresholds for RNA spot intensities
ths  = [50 50 21 50 101 50 50 50 50]; % Thresholds for RNA spot intensities % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS; defines which fluorecnet dyes have been imaged; 


%% Run the segmentation or RNA threshold determination
 if Seg == 1
     A1_segment_predefined_variables(Yth,Ych,ManTh,Ywin,Yim,task_ids,file_type,exp_date,strain,osmo,genepair,file_end,exp_names,positions_touse,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_num,img_stacks)
 end
 if RNA_thres == 1
     Run_RNA_detection_predefined_variables(thA,ths,outfile_prefix_RNA,Yth,Ych,ManTh,Ywin,Yim,task_ids,file_type,exp_date,strain,osmo,genepair,file_end,exp_names,positions_touse,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,img_num,node_nummax_int_thres,max_int_spot_th)
 end   
