%% Specify Image
for mm = 1:7
pairnum = 4;
task_id = mm;
counter = 6;
Cloud = 0;
Trans_show = 1;
min_RNA_size = 6;
max_RNA_size = 100;
size_filter = 0;

Ych = [0 1 1 0 1 0 0 0 1];                                                  % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
file_type = 'multiple' ;                                                    % 'single' or ' multiple'   % Defines if images have been saved as single files or if all files have been saved in a single file
exp_date = {'20180619','20180620','20180621','20180622','20180623','20180626','20180628','20180629'};                                                        % Date of the experiment and folder name on the hard disk
strain = 'F1-2-1';
osmo = '0.2MStep';
genepair = {'Xist-CY5-Maoa-TMR','Xist-CY5-Pdk3-TMR','Xist-CY5-RepA-TMR','Xist-CY5-Tsix-TMR'};
file_end = '_MMStack.ome.tif';

% Time points
exp_names = {'0d' '6hr' '12hr','1d','2d','5d','9d','9d'};% '6min' '8min' '10min' '12min' '14min',...
   % '15min' '18min' '20min' '22min' '25min' '40min' '50min' '120min' '122min' '129min'};

positions_touse = { ...
 [1:11];... %  0min 4                                                        %Change each number to the range of image numbers each timepoint has
 [1:11];... %  1min 4
 [1:11];... %  2min 4
  [1:11];... %  4min 4 
 % [1:11];... %  6min 6
 % [1:6];... %  8min 4
 % [1:9];... %  10min
 % [1:6];... %  12min
 % [1:4];... %  14min
 % [1:5];... %  16min
 % [1:7];... %  18min
 % [1:6];... %  20min 4
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

chi = find(Ych == 1);                                                       % Determines how many channels have been imaged
ch = size(chi,2);                                                           % Determines how many channels have been imaged

fish_dir = strcat('D:\',exp_date{task_id}); %[A]
exp_name = [exp_date{task_id} '_' exp_names{task_id} '_' strain  '_' genepair{pairnum} '_img'];
 filename = [fish_dir '\' exp_name '_' num2str(counter) '\' exp_name '_' num2str(counter) file_end];
   [stack, img_read] = tiffread2(filename);
        ImSt = img_read/ch;
 
        if Ych(2) == 1;                                                      % load DAPI    CHANGED TO 7 TEMPORARILY BECAUSE OF SOMETHING STRANGE                             
            ij = find(chi==2);                                                  %CHANGED TO 7 TEMPORARILY BECAUSE OF SOMETHING STRANGE
            DAPIim = [ij:ch:img_read];
            im_size = size(stack(1,DAPIim(1)).data);                        %BK 5/15/15 
            DAPI_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/15/15 resets the DAPI image so previous ones do not get incorporated
            for i = 1:ImSt;
                DAPI_ims(:,:,i) = stack(1,DAPIim(i)).data;
            end;
        else
        end;
                if Ych(9) == 1;                                                     % load TRANS
            ij = find(chi==9);
            TRANSim = [ij:ch:img_read];
            im_size = size(stack(1,TRANSim(1)).data);                        %BK 5/15/15
            TRANS_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 5/15/15 resets the TRANS image so previous ones do not get incorporated
            for i = 1:ImSt;
                TRANS_ims(:,:,i) = stack(1,TRANSim(i)).data;
            end;
        else
        end;
         if Ych(3) == 1;                                                     % CY5
                ij = find(chi==3);
                CY5im = [ij:ch:img_read];
                im_size = size(stack(1,CY5im(1)).data);                        %BK 12/10/15
                CY5_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    CY5_ims(:,:,i) = stack(1,CY5im(i)).data;
                end;
            else
                CY5_ims = NaN;
            end;
              if Ych(5) == 1;                                                     % TMR
                ij = find(chi==5);
                TMRim = [ij:ch:img_read];
                im_size = size(stack(1,TMRim(1)).data);                        %BK 12/10/15
                TMR_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 12/10/15 resets the TRANS image so previous ones do not get incorporated
                for i = 1:ImSt;
                    TMR_ims(:,:,i) = stack(1,TMRim(i)).data;
                end;
            else
                TMR_ims = NaN;
            end;
            exp_name
            
load(['D:\Dropbox (VU Basic Sciences)\Vanderbilt Computer\2018-07-12 mESC Segmentation (20170607 timecourse) Output\Lab_' exp_name num2str(counter)],'cells')
load(['D:\Dropbox (VU Basic Sciences)\Vanderbilt Computer\2018-07-12 mESC Segmentation (20170607 timecourse) Output\nuclei_' exp_name num2str(counter)])
if Cloud
load(['D:\Dropbox (VU Basic Sciences)\Vanderbilt Computer\2018-07-12 mESC Segmentation (20170607 timecourse) Output\mRNA_20180623_F1-2-1_Pdk3-TMR-th145_63_Xist-CY5-th65_2d_im' num2str(counter) '.mat'])
end
%% Processing

TMR_max = max(TMR_ims,[],3);
CY5_max = max(CY5_ims,[],3);    figure(5); clf; imshow(CY5_ims(:,:,30),[100 1000])
figure(6); clf; imshow(cells,[])
figure(5); clf
DAPI_max = max(DAPI_ims,[],3);
Nuc_border = Label_low(:,:,30);
figure(2); imshow(Nuc_border,[]);
NucBorder1 = bwmorph(Nuc_border,'remove')*100000;
figure(2); imshow(NucBorder1,[]);
if Cloud
CloudBorder1 = bwmorph(clouds1,'remove')*100000;
CloudBorder2 = imadjust(CloudBorder1,[0 1]);
end
          % figure(1); clf; imshow(CellRNAorg3(:,:,blah),[1 max(CellRNAorg2(:))*.5]); title(num2str(blah))
 %       figure(101);  clf; imshow(CellRNAorg4(:,:,blah),[]); title(num2str(blah))
%           figure(102); clf;; imshow(clouds(:,:,blah),[]); title(num2str(blah))
%            figure(103); clf; imshow(TMR3Dorg(:,:,blah),[]); title(num2str(blah))
           %%for a specific cell
%             CloudBorder1 = bwmorph(CellRNAorg4A(:,:,blah),'remove')*100000;
%             CloudBorder2 = imadjust(CloudBorder1,[0 0.5]);
            NucBorder2 = imadjust(NucBorder1,[0 0.5]);
            CellBorder1 = bwmorph(cells,'remove')*100000;
           % CellBorder2 = imadjust(CellBorder1,[0 0.1]);
                  CellBorder2 = imadjust(CellBorder1,[0 1]);
            RNA1 = CY5_max; %immultiply(CellRNAorg2(:,:,blah),CellRNAorg2(:,:,blah)>mm4-40);
            RNA2 = TMR_max;
            RNAs1 = RNA1/max(RNA1(:));
            RNAs2 = RNA2/max(RNA2(:));
            if size_filter
            RNAbw = RNAs1> prctile(RNAs1(:),95);
             RNA_normal = bwareaopen(RNAbw, min_RNA_size);                        % remove DAPI signal that are too small
            RNAs1 = immultiply(RNAs1,RNA_normal);
              RNAbw = RNAs2> prctile(RNAs2(:),95);
             RNA_normal = bwareaopen(RNAbw, min_RNA_size);                        % remove DAPI signal that are too small
            RNA_huge = bwareaopen(RNAbw, max_RNA_size);
            RNA_bw2 =  RNA_normal-RNA_huge;
            RNAs2 = immultiply(RNAs2,RNA_bw2);
            end
            %             dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
%             dapi_bw2 = dapi_normal - dapi_huge;                       
            RNAs1_2 = imadjust(RNAs1,[prctile(RNAs1(:),98) prctile(RNAs1(:),99.5)], [0 1],2);
            RNAs2_2 = imadjust(RNAs2,[prctile(RNAs2(:),98) prctile(RNAs2(:),99)], [0 1],2);
            figure(3); clf; imshow(RNAs1_2,[])
            figure(4); clf; imshow(RNAs2_2,[])
            %RNAs2 = imadjust(RNA2,[0 0.7]);
            Nuc1 = DAPI_max/max(DAPI_max(:));
            Nuc2 = imadjust(Nuc1,[prctile(Nuc1(:),50) prctile(Nuc1(:),80)],[0 1]);
            Trans1 = TRANS_ims(:,:,50);
            Trans2 = Trans1/max(Trans1(:));
            Trans3 = imadjust(Trans2,[prctile(Trans2(:),85) prctile(Trans2(:),99)],[0 1]);
            R = RNAs2_2+CellBorder2;
            B = NucBorder2+CellBorder2+Nuc2*.6;
            G = RNAs1_2+CellBorder2+Nuc2*.3;
            if Trans_show 
            R = R+Trans3;
            B = B+ Trans3;
            G = G+Trans3*.4; 
            end
            if Cloud
            R = R+CloudBorder2;
            B = B+ CloudBorder2;
            G = G+CloudBorder2;
            end
            RGB = cat(3,R,G,B);
           % RGB2 = imresize(RGB,3);
            figure(30); subplot(2,4,mm); imshow(RGB,[]);
            hold on
         pause(.5)
         counter
end
     %save(['RGB_image_stck_Cloud_' exppath '_Task' num2str(Task_Number) '_cell_' num2str(j)],'RGB_all');  