%% This file is to determine the threshold to count the number of mRNA spots.
%% For each image and at each time point a small area of the image is choosen to process a stack of images in 3D.
%% Changing the threshold up () or down () let you find the best threshold range.
%% For each experiment define which dyes have been used.
%% Also ajust the threshold range if needed.
%% The choosen thresholds are save in a file and can be used later for batch processing of images.

function Run_RNA_detection_predefined_variables_v2GUI(thA,ths,outfile_prefix_RNA,Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,node_num,max_int_thres,max_int_spot_th,images1,im_prefixes,img_stacks)
%% Input
% RNA_thresholds = zeros(6,65);
% counter6 = 1;                     %counter for RNA_thresholds
% if not(Ywin)
%     counter6 = (node_num)*sum(Ych(1:7))+1;
% end
ch = Ych(1);
 chi = find(Ych == 1);
% ch = size(chi,2);
max_int_thres = 1;
max_int_spot_th = 0;
images1
im_prefixes
diary 'thresh text2'
diary on
for im_num = 1:size(images1,2) 

%     try
S = im_num;
msgbox(im_prefixes{im_num});
% for S = task_ids
% for pairnum = pairnums
     S 
%     exp_name = [exp_date{S} '_' exp_names{S} '_' strain '_' genepair{pairnum} '_img'];
    counter = 1;  
%     if Ywin
%         img_num = positions_touse{S}                                      %gives loop all images if doing on windows (not in parallel)
%     end
    %% When switching exp date in loop
% if Ywin == 0;;
%     fish_dir = strcat('/gpfs23/scratch/neuertg/BenK/',exp_date{S}); %[ B ]
% else
%     fish_dir = strcat('C:\Users\keslerbk\Desktop\Microscopy\',exp_date{S}); %[A]
% end
%%
%     for counter=img_num
        counter
        h =sprintf('%01d',counter); 
        
        %% Load trans and dapi images
% %         fileDAPI = [outfile_prefix_seg 'Lab_' exp_name num2str(counter) '.mat']
%         fileDAPI = ['Lab_' im_prefixes{im_num} '.mat']
%         load(fileDAPI,'cells','trans_plane','-mat');

        %% Load micromanger images
        if true; %strcmp(file_type,'multiple');
    %    figure(1); clf; imshow(stack(1,7).data,[]);
%             if Ywin
%                 filename = [fish_dir '\' exp_name '_' num2str(counter) '\' exp_name '_' num2str(counter) file_end]
%             else
%                 filename = [fish_dir '/' exp_name '_' num2str(counter) '/' exp_name '_' num2str(counter) file_end]
%             end
%             [stack, img_read] = tiffread2(filename);
images1{im_num}
            [stack, img_read] = tiffread2(images1{im_num});
            ImSt = img_read/ch;
            im_size = size(stack(1,1).data);
         if false;%Ych(3); %Ych(1) == 1;                                                     % CY7
                ij = Ych(3);%ij = find(chi==1);
                TRANSim = [ij:ch:img_read];
                im_size = size(stack(1,TRANSim(1)).data);                        %BK 12/10/15
                TRANS_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    TRANS_ims(:,:,i) = stack(1,TRANSim(i)).data;
                end;
                cells = max(TRANS_ims(:,:,img_stacks(1):img_stacks(2)),[],3);
                trans_plane = NaN(im_size(1),im_size(2),ImSt);
            else
                TRANS_ims = NaN;
                cells = NaN;
                trans_plane = NaN;
            end;


        if false; %Ych(1) == 1;                                                     % CY7
                ij = find(chi==1);
                CY7im = [ij:ch:img_read];
                im_size = size(stack(1,CY7im(1)).data);                        %BK 12/10/15
                CY7_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    CY7_ims(:,:,i) = stack(1,CY7im(i)).data;
                end;
            else
                CY7_ims = NaN;
            end;

            if false; %Ych(2) == 1;                                                     % AF700
                ij = find(chi==2);
                AF700im = [ij:ch:img_read];
                im_size = size(stack(1,AF700im(1)).data);                        %BK 12/10/15
                AF700_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    AF700_ims(:,:,i) = stack(1,AF700im(i)).data;
                end;
            else
                AF700_ims = NaN;
            end;

            if Ych(4);                                                     % CY5
                ij = Ych(4) %ij = find(chi==3);
                CY5im = [ij:ch:img_read];
                im_size = size(stack(1,CY5im(1)).data);                        %BK 12/10/15
                CY5_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    CY5_ims(:,:,i) = stack(1,CY5im(i)).data;
                end;
                %CY5_ims = uint16(CY5_ims);
            else
                CY5_ims = NaN;
            end;

%             if Ych(4) == 1;                                                     % AF594
            if Ych(5)                                                     % AF594
                ij = Ych(5);% find(chi==4);
                AF594im = [ij:ch:img_read];
                im_size = size(stack(1,AF594im(1)).data);                        %BK 12/10/15
                AF594_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    AF594_ims(:,:,i) = stack(1,AF594im(i)).data;
                end;
                %AF594_ims = uint16(AF594_ims);
            else
                AF594_ims = NaN;
            end;


            if Ych(6);                                                     % TMR
                ij = Ych(6);%ij = find(chi==5);
                TMRim = [ij:ch:img_read];
                im_size = size(stack(1,TMRim(1)).data);                        %BK 12/10/15
                TMR_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    TMR_ims(:,:,i) = stack(1,TMRim(i)).data;
                end;
                %TMR_ims = uint16(TMR_ims);
            else
                TMR_ims = NaN;
            end;


            if false; %Ych(6) == 1;                                                     % YFP
                ij = find(chi==6);
                YFPim = [ij:ch:img_read];
                im_size = size(stack(1,YFPim(1)).data);                        %BK 12/10/15
                YFP_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    YFP_ims(:,:,i) = stack(1,YFPim(i)).data;
                end;
            else
                YFP_ims = NaN;
            end;

            if false; %Ych(7) == 1;                                                     % GFP / AF488
                ij = find(chi==7);
                GFPim = [ij:ch:img_read];
                im_size = size(stack(1,GFPim(1)).data);                        %BK 12/10/15
                GFP_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    GFP_ims(:,:,i) = stack(1,GFPim(i)).data;
                end;
            else
                GFP_ims = NaN;
            end;

        
        elseif strcmp(file_type,'single');
            S
            if Ych(1) == 1;                                                     % CY7
                filebase = 'img_000000000_1-CY7_';
                i = 1;
                for i = 1:50;
                    h =sprintf('%03d',i-1); 
                    fileN = [fish_dir exp_date{S} '\' exp_name num2str(counter) '\' filebase h '.tif'];
                    if exist(fileN) == 2;
                        ims = tiffread2(fileN);
                        CY7_ims(:,:,i) = ims.data;
                        ImSt = i;
                    else;
                    end;
                end;
            else
                CY7_ims = NaN;
            end;

            if Ych(2) == 1;                                                     % AF700
                filebase = 'img_000000000_2-AF700_';
                i = 1;
                for i = 1:50;
                    h =sprintf('%03d',i-1); 
                    fileN = [fish_dir exp_date{S} '\' exp_name num2str(counter) '\' filebase h '.tif'];
                    if exist(fileN) == 2;
                        ims = tiffread2(fileN);
                        AF700_ims(:,:,i) = ims.data;
                        ImSt = i;
                    else;
                    end;
                end;
            else
                AF700_ims = NaN;
            end;

            if Ych(3) == 1;                                                     % CY5
                filebase = 'img_000000000_3-CY5_';
                i = 1;
                for i = 1:50;
                    h =sprintf('%03d',i-1); 
                    fileN = [fish_dir exp_date{S} '\' exp_name num2str(counter) '\' filebase h '.tif'];
                    if exist(fileN) == 2;
                        ims = tiffread2(fileN);
                        CY5_ims(:,:,i) = ims.data;
                        ImSt = i;
                    else;
                    end;
                end;
            else
                CY5_ims = NaN
            end;

            if Ych(4) == 1;                                                     % AF594
                filebase = 'img_000000000_4-AF594_';
                i = 1;
                for i = 1:50;
                    h =sprintf('%03d',i-1); 
                    fileN = [fish_dir exp_date{S} '\' exp_name num2str(counter) '\' filebase h '.tif'];
                    if exist(fileN) == 2;
                        ims = tiffread2(fileN);
                        AF594_ims(:,:,i) = ims.data;
                        ImSt = i;
                    else;
                    end;
                end;
            else
                AF594_ims = NaN;
            end;

            if Ych(5) == 1;                                                     % TMR
                filebase = 'img_000000000_5-TMR_';
                i = 1;
                for i = 1:50;
                    h =sprintf('%03d',i-1); 
                    fileN = [fish_dir exp_date{S} '\' exp_name num2str(counter) '\' filebase h '.tif'];
                    if exist(fileN) == 2;
                        ims = tiffread2(fileN);
                        TMR_ims(:,:,i) = ims.data;
                        ImSt = i;
                    else;
                    end;
                end;
            else
                TMR_ims = NaN;
            end;

            if Ych(6) == 1;                                                     % YFP
                filebase = 'img_000000000_6-YFP_';
                i = 1;
                for i = 1:50;
                    h =sprintf('%03d',i-1); 
                    fileN = [fish_dir exp_date{S} '\' exp_name num2str(counter) '\' filebase h '.tif'];
                    if exist(fileN) == 2;
                        ims = tiffread2(fileN);
                        YFP_ims(:,:,i) = ims.data;
                        ImSt = i;
                    else;
                        
                    end;
                end;
            else
                YFP_ims = NaN;
            end;

            if Ych(7) == 1;                                                     % GFP / AF488
                filebase = 'img_000000000_6-GFP_';
                i = 1;
                for i = 1:50;
                    h =sprintf('%03d',i-1); 
                    fileN = [fish_dir exp_date{S} '\' exp_name num2str(counter) '\' filebase h '.tif'];
                    if exist(fileN) == 2;
                        ims = tiffread2(fileN);
                        GFP_ims(:,:,i) = ims.data;
                        ImSt = i;
                    else;
                    end;
                end;
            else
                GFP_ims = NaN;
            end;
        else
            'test'
            return;
        end;

        %load([outfile_prefix_RNA filename],'PixRNA','PixRNA1d')
        %%%%Graphing Code%%%%
        %figure(88); imshow(PixRNA(1:2048,1:2048,8,1,2),[])
%         genepairs = {'SCR-TMR-SCR-CY5-SCR-AF594','no-probes','STL1-TMR-GPD1-CY5-GPP2-AF594'};
%         colors = {'r','b','g'};
%         styles = {'--',':','-'}; 
%         figure(89);clf
%         for m = 1:3
%             for k = 1:3               
%                 genepair = genepairs{k};
%                 filename = ['PixRNA_' exp_date '_' osmo '_' strain '_' genepair '.mat']
%                 load([outfile_prefix_RNA filename])
%                 if k == 3
%                     temp12 = PixRNA1d(1:2048*2048,1,1,m);
%                 else
%                     temp12 = PixRNA1d(1:2048*2048,1,1,m);
%                 end
%                 minp = min(temp12);
%                 maxp = max(temp12);
%                 binnum12 = 1000;
%                 xval = 50:2:800;
%                 xvalp = 51:2:799;
%                 yvalp = [0];
%                 %         xval = minp:(maxp-minp)/binnum12:maxp;
%                 %         xvalp = minp+(maxp-minp)/binnum12/2:ceil((maxp-minp)/binnum12):maxp-(maxp-minp)/binnum12/2;
%                 for j = 1:size(xval,2)-1
%                     yvalp(1,j) = size(find(temp12 <= xval(j+1)),1)-size(find(temp12 <= xval(j)),1);
%                 end
%                 %       yvalp(1,j) = yvalp(1,j) + size(find(temp12 == xval(j)),1);          %Adds the last values
%                 yvalp = yvalp/sum(yvalp);
%                 plot(xvalp,yvalp,'Color',colors{m},'Linestyle',styles{k})
%                 hold on
%             end
%         end
%         title('Comparison of Scrambled Probes, No Probes, and Gene Probes')
%         xlabel('Pixel Intensity');ylabel('Probability within a Cell');
%         legend('scrambled probes CY5','no probes CY5','GPD1 CY5','scrambled probes AF594','no probes AF594','GPP2 AF594','scrambled probes TMR','no probes TMR','STL1 TMR')
       %%%%%%%%% 
    im_size = size(stack(1,1).data);
    whos
    clear stack
    for channels1 = 4:6
        if Ych(4)& channels1 == 4
            thresh_img = CY5_ims;
        elseif Ych(5)& channels1 == 5
            thresh_img = AF594_ims;
        elseif Ych(6)& channels1 == 6
            thresh_img = TMR_ims;
        end
                   img_name = im_prefixes{im_num};
           tif_path = fish_dir;
           rna_ch = channels1;
           light_ch = ch;
           total_ch = ch;
           out_dir = outfile_prefix_seg;
           t_min = thA(1);
           t_max = thA(end);
           cellseg_path = outfile_prefix_seg;
           ztrim = 5;
           ctrl_path = '';
           ctrl_ch = [];
           ctrl_chcount = [];
           ttune_winsize = 0;
           ttune_wscorethresh = -1;
           overwrite_output = 1;
                      guimode = 1;
%            preloaded_images = thresh_img;
           addpath('./core');
           preloaded_images = ImageChannelSet;
           preloaded_images.dat_rna_sample = thresh_img;
    rna_spot_run = Adapter_RNASpots(img_name, tif_path, rna_ch, light_ch, total_ch,...
    out_dir, t_min, t_max, ztrim, cellseg_path, ctrl_path, ctrl_ch, ctrl_chcount, ttune_winsize, ttune_wscorethresh,...
    overwrite_output, guimode, preloaded_images)
   files = char(strcat('mRNA_TH_',im_prefixes{im_num},'channel,',num2str(channels1),'.mat')); %strcat(exp_date,'_',ch,'_th',ths,'_p',ps,'_im',h)
%     save([outfile_prefix_RNA char(files)],'thCY7Man','thAF700Man','thCY5Man','thAF594Man','thTMRMan','thYFPMan','thGFPMan'); % save the analyzed files
    save([outfile_prefix_RNA char(files)],'rna_spot_run')
%         Main_RNASpotDetectv2(thresh_img, ImSt,fish_dir, im_prefixes{im_num},outfile_prefix_RNA, ch,thA(1), thA(end), true,true);
    
    end
%     [thCY7,thAF700,thCY5,thAF594,thTMR,thYFP,thGFP] = AB_FindThreshold_TMR_AF594_CY5_B(im_size,counter,cells,trans_plane,S,CY7_ims,AF700_ims,CY5_ims,AF594_ims,TMR_ims,YFP_ims,GFP_ims,thA,ths,Ych,max_int_thres,max_int_spot_th);
    % collect all the picked thresholds for CY7, AF700, CY5, AF594, TMR, YFP,GFP
%     thCY7Man(1:3,S,counter) = thCY7;
%     thAF700Man(1:3,S,counter) = thAF700;
%     thCY5Man(1:3,S,counter) = thCY5;
%     thAF594Man(1:3,S,counter) = thAF594;
%     thTMRMan(1:3,S,counter) = thTMR; 
%     thYFPMan(1:3,S,counter) = thYFP; 
%     thGFPMan(1:3,S,counter) = thGFP; 
    clear TMR_ims CY5_ims AF594_ims
     whos

%     files = char(strcat('mRNA_TH_',exp_name,num2str(counter),'.mat')); %strcat(exp_date,'_',ch,'_th',ths,'_p',ps,'_im',h)
 
    %     ths(3) = thCY5;
%     ths(1) = thTMR;
%     ths(2) = thAF594;
%     ths(3) = thCY5;
%     end
%     end;
end
diary off
% filename = ['RNA_thresholds_' exp_date '_' osmo '_' strain '_' genepair '.mat']
% if not(Ywin) & exist([filename,'RNA_thresholds'])
%     temp_RNA_thresholds = RNA_thresholds(1:7,(node_num)*sum(Ych(1:7))+1:(node_num)*sum(Ych(1:7))); 
%     load([outfile_prefix_RNA filename],'RNA_thresholds');
%     RNA_thresholds(1:7,(node_num)*sum(Ych(1:7))+1:(node_num)*sum(Ych(1:7))) = temp_RNA_thresholds;
%     save([outfile_prefix_RNA filename],'RNA_thresholds');
% else
%     save([outfile_prefix_RNA filename],'RNA_thresholds')
% end


%% Saves threshold values for each image and each sample separately
% thTMRman, thAF594Man, thCY5Man     
            
            
            
 
