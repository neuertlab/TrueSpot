%% This code segment yeast cells using a stack of trans and DAPI images. 
%% The code can automatically segment the cell boundary in 2D.
%% The cell nucleus is segmented in 3D using three different threshold for the DAPI signas in each cell. 
%% This allows to estimate the change of the threshold by +/-10% and its effect on the nuclear size. 
%% We use the difference in the segmented nucleus to estimate the error in the quantification on the nuclear RNA spots.

function A1_segment_predefined_variables(Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_stacks,yeast_seg,images1,im_prefixes,dxy1,Diffstack,DAPI_images,TRANS_images)
% matlabpool
%change line 179 if you want to segment a different way (threshold vs adding
%together thresholds)
images1 
file_type
counter2 = 1                                                                    %counter for where to put in max nuclei                                                               % Show images of nucleus and trans segmentation
  max_nuclei = zeros(4,66);
fail_count = 0;         %Will count how many images failed
diff_channel_count = 0; %will count how many times the number of channels was incorrect  
if ManTh
    Seg_thres= 1;
end
for im_num = 1:size(images1,2)                                                            % Specify the job (experiment and time point) that needs to be processes. comment out if run in a loop or on the cluster
task_id = 1;
% try

%% Channels
%%%Old way that was not generalizable
%  chi = find(Ych == 1);                                                       % Determines how many channels have been imaged
% ch = size(chi,2);                                                           % Determines how many channels have been imaged
%%%

%%%generalizable way
%in Ych, first is the total number of channels, second is the marker slice, and third is the boundary slice
ch = Ych(1);                                                           % Determines how many channels have been imaged
%%%

%% When switching exp date in loop
% if Ywin == 0;;
%     fish_dir = strcat('/gpfs23/scratch/neuertg/BenK/',exp_date{task_id}); %[ B ]
% else
%     fish_dir = strcat('C:\Users\keslerbk\Desktop\Microscopy\',exp_date{task_id}); %[A]
% end
%% Chooses one genepair after another

% for pairnum = pairnums
% counter2

%% Chooses one image after another;
% counter = 1;
% if Ywin
%     img_num = positions_touse{task_id}                                      %gives loop all images if doing on windows (not in parallel)
% end
% for counter=img_num

    %% Defines the sample names (may be different at every experiment)
%     exp_name = [exp_date{task_id} '_' exp_names{task_id} '_' strain  '_' genepair{pairnum} '_img'];
% 
%     counter
%     tic
%     h =sprintf('%03d',counter);                                             % determine the image index
%     
%% Load micromanger images
if strcmp(file_type,'Multiple')                                        % load images that are saved as image stacks
    %filename = [fish_dir '/' exp_name '_' num2str(counter) '/' exp_name '_MMImages.ome.tif']
    %         if Ywin == 0;
    %             filename = [fish_dir '/' exp_name '_' num2str(counter) '/' exp_name '_' num2str(counter) file_end];
    %         else
    %             filename = [fish_dir '\' exp_name '_' num2str(counter) '\' exp_name '_' num2str(counter) file_end];
    %         end;
    %f = msgbox('Loading Image');
    if Diffstack
        DAPI_images{im_num}
        [stack, img_read] = tiffread2(DAPI_images{im_num});
        ImSt = img_read;
        DAPIim = [1:img_read];
        im_size = size(stack(1,DAPIim(1)).data);                        %BK 5/15/15
        DAPI_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/15/15 resets the DAPI image so previous ones do not get incorporated
        for i = 1:ImSt;
            DAPI_ims(:,:,i) = stack(1,DAPIim(i)).data;
        end
        TRANS_images{im_num}
        [stack, img_read] = tiffread2(TRANS_images{im_num});
        ImSt = img_read;
        TRANSim = [1:img_read];
        im_size = size(stack(1,TRANSim(1)).data);                        %BK 5/15/15
        TRANS_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 5/15/15 resets the TRANS image so previous ones do not get incorporated
        for i = 1:ImSt;
            TRANS_ims(:,:,i) = stack(1,TRANSim(i)).data;
        end;
        clear stack  
    else
        [stack, img_read] = tiffread2(images1{im_num});
        ImSt = img_read/ch;
        
        %         if Ych(2) == 1;                                                      % load DAPI    CHANGED TO 7 TEMPORARILY BECAUSE OF SOMETHING STRANGE
        %             ij = find(chi==2);                                                  %CHANGED TO 7 TEMPORARILY BECAUSE OF SOMETHING STRANGE
        if Ych(2);
            ij = Ych(2);
            DAPIim = [ij:ch:img_read];
            im_size = size(stack(1,DAPIim(1)).data);                        %BK 5/15/15
            DAPI_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/15/15 resets the DAPI image so previous ones do not get incorporated
            for i = 1:ImSt;
                DAPI_ims(:,:,i) = stack(1,DAPIim(i)).data;
            end;
        else
        end;
        
        %         if Ych(9) == 1;                                                     % load TRANS
        %             ij = find(chi==9);
        if Ych(3)                                                     % load TRANS
            ij = Ych(3);
            TRANSim = [ij:ch:img_read];
            im_size = size(stack(1,TRANSim(1)).data);                        %BK 5/15/15
            TRANS_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 5/15/15 resets the TRANS image so previous ones do not get incorporated
            for i = 1:ImSt;
                TRANS_ims(:,:,i) = stack(1,TRANSim(i)).data;
            end;
        else
        end;
    end
    
    
elseif strcmp(file_type,'single');                                      % load images that are saved as single image
    %  counter
    
    if Ych(8) == 1;                                                    % load DAPI
        filebase = 'img_000000000_1-DAPI_';                             % compare to file name;
        i = 1;
        for i = 1:50;
            h =sprintf('%03d',i-1);
            %fileN = [fish_dir '/' exp_name num2str(counter) '/' filebase h '.tif'];
            if Ywin == 0;;
                fileN = [fish_dir '/' exp_name num2str(counter) '/' filebase h '.tif'];
            else
                fileN = [fish_dir '\' exp_name num2str(counter) '\' filebase h '.tif'];
            end
            if exist(fileN) == 2;
                ims = tiffread2(fileN);
                DAPI_ims(:,:,i) = uint16(ims.data);
                ImSt = i;
            else;
            end;
        end;
    else
    end;
    
    if Ych(9) == 1;                                                    % load TRANS
        filebase = 'img_000000000_7-TRANS_';                            % compare to file name
        i = 1;
        for i = 1:50;
            h =sprintf('%03d',i-1);
            %fileN = [fish_dir '/' exp_name num2str(counter) '/' filebase h '.tif'];
            if Ywin == 0;;
                fileN = [fish_dir '/' exp_name num2str(counter) '/' filebase h '.tif'];
            else
                fileN = [fish_dir '\' exp_name num2str(counter) '\' filebase h '.tif'];
            end
            if exist(fileN) == 2;
                ims = tiffread2(fileN);
                TRANS_ims(:,:,i) = uint16(ims.data);
                ImSt = i;
            else;
            end;
        end;
    else
    end;  
else
    return;
end;
    
clear stack;
DAPI_ims = DAPI_ims(round(size(DAPI_ims,1)/4):size(DAPI_ims,1)-round(size(DAPI_ims,1)/4),round(size(DAPI_ims,2)/4):size(DAPI_ims,2)-round(size(DAPI_ims,2)/4),:);
TRANS_ims = TRANS_ims(round(size(TRANS_ims,1)/4):size(TRANS_ims,1)-round(size(TRANS_ims,1)/4),round(size(TRANS_ims,2)/4):size(TRANS_ims,2)-round(size(TRANS_ims,2)/4),:);
    %% Determine  DAPI threshold
    %counter
    fileTH = ['nuclei_TH_' im_prefixes{im_num} '.mat']
    if false;%exist([outfile_prefix_seg fileTH],'file'); % test if the file exist allready     
       load([outfile_prefix_seg fileTH])   % load the existing file
    else                                                                    % if the file does not exist, process DAPI images
       if Seg_thres == 1
           [dapi_threshold,dapi_label] = B1_autosegment_nuclei8_thres(DAPI_ims,Yim,ManTh,min_nucleus_size,max_nucleus_size);  
            DAPI_ims_added = zeros(size(dapi_label));
            DAPI_ims_added = uint8(DAPI_ims_added);
            dapi_label_low1 = dapi_label;
       else
           [dapi_threshold,dapi_label,counter2,max_nuclei,DAPI_ims_added,dapi_label_low1] = B1_autosegment_nuclei8_adding_indiv(DAPI_ims,Yim,ManTh,counter2,max_nuclei,min_nucleus_size,max_nucleus_size,Ywin,dxy1);
       end
       save([outfile_prefix_seg fileTH],'dapi_threshold','dapi_label','DAPI_ims_added','-mat');     % Save the segmented nuclei
    end;
   % save(strcat(outfile_prefix_seg, 'max_nuclei', exp_date{task_id}), 'max_nuclei')

    if Yth~= 1;
       
    %% Segment DAPI image
    fileDAPI = ['nuclei_' im_prefixes{im_num} '.mat']
    if false;%exist([outfile_prefix_seg fileDAPI],'file');
       load([outfile_prefix_seg fileDAPI]);   % load the existing file
    else;
       if yeast_seg
         [Label_low,Label_mid,Label_hi,nuclei,Nuc_int,Nuc_vol,Maj_axis,Min_axis] = B1_autosegment_nuclei5_yeast(DAPI_ims,dapi_label,Yim,dapi_threshold,min_nucleus_size,max_nucleus_size,dxy1); % Segment DAPI images in 3D with three different thresholds
       save([outfile_prefix_seg fileDAPI],'Label_low','Label_mid','Label_hi','nuclei','Nuc_int','Nuc_vol','Maj_axis','Min_axis','-mat','-v7.3');     % Save the segmented nuclei      
       else
        [Label_low,Label_mid,Label_hi,nuclei,Nuc_int,Nuc_vol,Maj_axis,Min_axis] = B1_autosegment_nuclei5(DAPI_ims,dapi_label,Yim,dapi_threshold,min_nucleus_size,max_nucleus_size,dxy1,dapi_label_low1); % Segment DAPI images in 3D with three different thresholds
       save([outfile_prefix_seg fileDAPI],'Label_low','Label_mid','Label_hi','nuclei','Nuc_int','Nuc_vol','Maj_axis','Min_axis','-mat','-v7.3');     % Save the segmented nuclei      
       end
     end;

    %% Segment TRANS image
    fileTRANS = ['Lab_' im_prefixes{im_num} '.mat'];
    if false; %exist([outfile_prefix_seg fileTRANS],'file');
        load([outfile_prefix_seg fileTRANS],'-mat');
    else;
        if yeast_seg
            %nuclei = max(Label_low,[],3) commented 4/7/16 BK
        [cells,plane,CellInfo,trans_plane,CellImage]= B2_autosegment_cells_new_yeast(TRANS_ims,nuclei,segment_mode,Yim,min_cell_size,max_cell_size,dapi_label); % segment the trans image    
        save([outfile_prefix_seg fileTRANS],'cells','plane','segment_mode','CellInfo','trans_plane','CellImage','-mat','-v7.3'); % save the segmented cells
        else  
        [cells,plane,CellInfo,trans_plane,CellImage]= B2_autosegment_cells_new(TRANS_ims,nuclei,segment_mode,Yim,min_cell_size,max_cell_size,dapi_label,img_stacks); % segment the trans image    
        save([outfile_prefix_seg fileTRANS],'cells','plane','segment_mode','CellInfo','trans_plane','CellImage','-mat','-v7.3'); % save the segmented cells
        end
    end;
    else;
    end;
% catch
% fail_count = fail_count+1;
% if not(Diffstack) & round(ImSt) ~=  ImSt   %Checks if the number of channels was incorrect (for multi-stack)
%     ImSt
%     diff_channel_count = diff_channel_count+1;
% end
% end
end
if diff_channel_count == size(images1,2) 
    msgbox('All image analyses failed due to incorrect channel number');
elseif fail_count == im_num
    msgbox('All image analyses failed');
end
end

