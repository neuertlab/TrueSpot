%% This file is to determine the threshold to count the number of mRNA spots.
%% For each image and at each time point a small area of the image is choosen to process a stack of images in 3D.
%% Changing the threshold up () or down () let you find the best threshold range.
%% For each experiment define which dyes have been used.
%% Also ajust the threshold range if needed.
%% The choosen thresholds are save in a file and can be used later for batch processing of images.
clear
%% Input
RNA_thresholds = zeros(6,65);
counter6 = 1;
for x = [1]
if x
    thA = [1,10,1000]; % Thresholds for RNA spot intensities
else
    thA = [1,2,200]; % Thresholds for RNA spot intensities
end
ths  = [50 50 21 50 101 50 50 50 50]; % Thresholds for RNA spot intensities % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS; defines which fluorecnet dyes have been imaged; 
file_type = 'multiple';                                                     % 'single' or 'multiple'
S = 1;
Im = 1;
maxSample = 15; % maxium number of sample in the experiment
maxIm = 4;
Ych = [0 0 1 0 1 0 0 1 1]  %third is 1                                                 % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS; defines which fluorecnet dyes have been imaged; 
% Ych = [0 0 1 1 1 0 0 1 1] % YCY7 AF700 CY5 AF594 TMR YFP GFP\AF488 DAPI TRANS
chi = find(Ych == 1);
ch = size(chi,2);

globalpath = 'C:\Users\keslerbk\Documents\Guoliang''Data\'; % main directory for microscopy images
exppath = '20130914'; % Day of experiments
strain = 'WT';
osmo = '0.02Mramp-stop10min';
genepair = 'STL1-TMR-1000_CTT1-CY5-2000';

exp_names = {'0min' '1min' '2min' '4min' '6min' '8min' '10min' '15min',...
    '20min' '25min' '30min' '35min' '40min' '45min' '50min' '55min' '60min'}; %

filepath1 = 'C:\Users\keslerbk\Documents\'; % segmented cells are stored here
positions_touse = { ...
 [1:4];... %  00min
 [1:4];... %  01min
 [1:4];... %  02min
 [1:4];... %  04min
 [1:4];... %  06min
 [1:4];... %  08min
 [1:4];... %  10min
 [1:4];... %  15min
 [1:4];... %  20min
 [1:4];... %  25min
 [1:4];... %  30min
 [1:4];... %  35min
 [2:4];... %  40min
 [1:4];... %  45min
 [1:3];... %  50min
 [1:4];... %  55min
 [1:4];... %  60min
};

for S = 1:17%1:maxSample; % goes through all the samples
    S 
    exp_name = [exppath '_' osmo '_' strain '_' exp_names{S} '_' genepair '_img'];
    counter = 1;    
    for counter=positions_touse{S}
        counter
        h =sprintf('%01d',counter); 
        
        %% Load trans and dapi images
        fileDAPI = [filepath1 exppath '\SegFiles\Lab_' exp_name num2str(counter) '.mat']
        load(fileDAPI,'cells','trans_plane','-mat');

        %% Load micromanger images
        if strcmp(file_type,'multiple');
    %    figure(1); clf; imshow(stack(1,7).data,[]);

            filename = [globalpath exppath '\' exp_name '_' num2str(counter) '\' exp_name '_MMImages.ome.tif']
            [stack, img_read] = tiffread2(filename);
            ImSt = img_read/ch;

            if Ych(1) == 1;                                                     % CY7
                ij = find(chi==1);
                CY7im = [ij:ch:img_read];
                im_size = size(stack(1,CY7im(1)).data);                        %BK 5/21/15 
                CY7_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/21/15 resets the CY7 image so previous ones do not get incorporated
                for i = 1:ImSt;
                    CY7_ims(:,:,i) = stack(1,CY7im(i)).data;
                end;
            else
                CY7_ims = NaN;
            end;

            if Ych(2) == 1;                                                     % AF700
                ij = find(chi==2);
                AF700im = [ij:ch:img_read];
                im_size = size(stack(1,AF700im(1)).data);                        %BK 5/21/15 
                AF700_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/21/15 resets the AF700 image so previous ones do not get incorporated
                for i = 1:ImSt;
                    AF700_ims(:,:,i) = stack(1,AF700im(i)).data;
                end;
            else
                AF700_ims = NaN;
            end;

            if Ych(3) == 1;                                                     % CY5
                ij = find(chi==3);
                CY5im = [ij:ch:img_read];
                im_size = size(stack(1,CY5im(1)).data);                        %BK 5/21/15 
                CY5_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/21/15 resets the CY5 image so previous ones do not get incorporated
                for i = 1:ImSt;
                    CY5_ims(:,:,i) = stack(1,CY5im(i)).data;
                end;
            else
                CY5_ims = NaN;
            end;

            if Ych(4) == 1;                                                     % AF594
                ij = find(chi==4);
                AF594im = [ij:ch:img_read];
                im_size = size(stack(1,AF594im(1)).data);                        %BK 5/21/15 
                AF594_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/21/15 resets the AF594 image so previous ones do not get incorporated
                for i = 1:ImSt;
                    AF594_ims(:,:,i) = stack(1,AF594im(i)).data;
                end;
            else
                AF594_ims = NaN;
            end;


            if Ych(5) == 1;                                                     % TMR
                ij = find(chi==5);
                TMRim = [ij:ch:img_read];
                im_size = size(stack(1,TMRim(1)).data);                        %BK 5/21/15 
                TMR_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/21/15 resets the TMR image so previous ones do not get incorporated
                for i = 1:ImSt;
                    TMR_ims(:,:,i) = stack(1,TMRim(i)).data;
                end;
            else
                TMR_ims = NaN;
            end;


            if Ych(6) == 1;                                                     % YFP
                ij = find(chi==6);
                YFPim = [ij:ch:img_read];
                im_size = size(stack(1,YFPim(1)).data);                        %BK 5/21/15 
                YFP_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/21/15 resets the YFP image so previous ones do not get incorporated
                for i = 1:ImSt;
                    YFP_ims(:,:,i) = stack(1,YFPim(i)).data;
                end;
            else
                YFP_ims = NaN;
            end;

            if Ych(7) == 1;                                                     % GFP / AF488
                ij = find(chi==7);
                GFPim = [ij:ch:img_read];
                im_size = size(stack(1,GFPim(1)).data);                        %BK 5/21/15 
                GFP_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/21/15 resets the GFP image so previous ones do not get incorporated
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
                    fileN = [globalpath exppath '\' exp_name num2str(counter) '\' filebase h '.tif'];
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
                    fileN = [globalpath exppath '\' exp_name num2str(counter) '\' filebase h '.tif'];
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
                    fileN = [globalpath exppath '\' exp_name num2str(counter) '\' filebase h '.tif'];
                    if exist(fileN) == 2;
                        ims = tiffread2(fileN);
                        CY5_ims(:,:,i) = ims.data;
                        ImSt = i;
                    else;
                    end;
                end;
            else
                CY5_ims = NaN;
            end;

            if Ych(4) == 1;                                                     % AF594
                filebase = 'img_000000000_4-AF594_';
                i = 1;
                for i = 1:50;
                    h =sprintf('%03d',i-1); 
                    fileN = [globalpath exppath '\' exp_name num2str(counter) '\' filebase h '.tif'];
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
                    fileN = [globalpath exppath '\' exp_name num2str(counter) '\' filebase h '.tif'];
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
                    fileN = [globalpath exppath '\' exp_name num2str(counter) '\' filebase h '.tif'];
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
                    fileN = [globalpath exppath '\' exp_name num2str(counter) '\' filebase h '.tif'];
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

    [thCY7,thAF700,thCY5,thAF594,thTMR,thYFP,thGFP,RNA_thresholds,counter6] = AB_FindThreshold_TMR_AF594_CY5_B_BK_auto(counter,cells,trans_plane,S,CY7_ims,AF700_ims,CY5_ims,AF594_ims,TMR_ims,YFP_ims,GFP_ims,thA,ths,Ych,RNA_thresholds,counter6);
    % collect all the picked thresholds for CY7, AF700, CY5, AF594, TMR, YFP,GFP
    thCY7Man(1:3,S,counter) = thCY7;
    thAF700Man(1:3,S,counter) = thAF700;
    thCY5Man(1:3,S,counter) = thCY5;
    thAF594Man(1:3,S,counter) = thAF594;
    thTMRMan(1:3,S,counter) = thTMR; 
    thYFPMan(1:3,S,counter) = thYFP; 
    thGFPMan(1:3,S,counter) = thGFP; 
    
    files = char(strcat('mRNA_TH_',exp_name,num2str(counter),'.mat')); %strcat(exppath,'_',ch,'_th',ths,'_p',ps,'_im',h)
    save(char(files),'thCY7Man','thAF700Man','thCY5Man','thAF594Man','thTMRMan','thYFPMan','thGFPMan'); % save the analyzed files
%     ths(3) = thCY5;
%     ths(1) = thTMR;
%     ths(2) = thAF594;
%     ths(3) = thCY5;
filename = ['RNA_thresholds_' exppath]
save(filename,'RNA_thresholds')
    end;
end;
end;

%% Saves threshold valuse for each image and each sample separately
% thTMRman, thAF594Man, thCY5Man     
            
            
            
 