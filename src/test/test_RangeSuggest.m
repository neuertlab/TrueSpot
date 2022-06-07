%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
addpath('./core');

% ========================== Image Channels ==========================

img_paths = cell(48,1);
i = 1;

%----- mESC Set 1
img_paths{i,1} = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d']; i = i+1; %1
img_paths{i,1} = [ImgDir '\data\preprocess\feb2018\Xist_CY5\Xist\Xist-CY5_IMG1_all_3d']; i = i+1; %2

%----- mESC Set 2
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\1Day\Tsix\mESC_1d_Tsix_all_3d']; i = i+1; %3
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\1Day\Xist\mESC_1d_Xist_all_3d']; i = i+1; %4
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\2Day\Tsix\mESC_2d_Tsix_all_3d']; i = i+1; %5
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\2Day\Xist\mESC_2d_Xist_all_3d']; i = i+1; %6
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\3Day\Tsix\mESC_3d_Tsix_all_3d']; i = i+1; %7
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\3Day\Xist\mESC_3d_Xist_all_3d']; i = i+1; %8

%----- yeast RNA
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E1R2\STL1\E1R2-STL1-TMR_all_3d']; i = i+1; %9
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E1R2\CTT1\E1R2-CTT1-CY5_all_3d']; i = i+1; %10
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R1\STL1\E2R1-STL1-TMR_all_3d']; i = i+1; %11
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R1\CTT1\E2R1-CTT1-CY5_all_3d']; i = i+1; %12
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\STL1\E2R2I3-STL1-TMR_all_3d']; i = i+1; %13
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\CTT1\E2R2I3-CTT1-CY5_all_3d']; i = i+1; %14
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\STL1\E2R2I5-STL1-TMR_all_3d']; i = i+1; %15
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\CTT1\E2R2I5-CTT1-CY5_all_3d']; i = i+1; %16

%----- yeast proteins
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d']; i = i+1; %17
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d']; i = i+1; %18
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\4M1m_img3\Msb2_04M_1m_img3_GFP_all_3d']; i = i+1; %19
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\4M5m_img2\Msb2_04M_5m_img2_GFP_all_3d']; i = i+1; %20
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\NoSalt_img1\Msb2_NoSalt_img1_GFP_all_3d']; i = i+1; %21

img_paths{i,1} = [ImgDir '\data\preprocess\opy2\2M1m_img3\Opy2_02M_1m_img3_GFP_all_3d']; i = i+1; %22
img_paths{i,1} = [ImgDir '\data\preprocess\opy2\2M5m_img2\Opy2_02M_5m_img2_GFP_all_3d']; i = i+1; %23
img_paths{i,1} = [ImgDir '\data\preprocess\opy2\4M1m_img1\Opy2_04M_1m_img1_GFP_all_3d']; i = i+1; %24
img_paths{i,1} = [ImgDir '\data\preprocess\opy2\4M5m_img3\Opy2_04M_5m_img3_GFP_all_3d']; i = i+1; %25
img_paths{i,1} = [ImgDir '\data\preprocess\opy2\NoSalt_img2\Opy2_NoSalt_img2_GFP_all_3d']; i = i+1; %26

%----- Xist/Tsix + histones
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I4\Ch2\Histone_D0_img4_ch2_all_3d']; i = i+1; %27
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I4\Ch3\Histone_D0_img4_ch3_all_3d']; i = i+1; %28
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I4\Ch4\Histone_D0_img4_ch4_all_3d']; i = i+1; %29
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I6\Ch2\Histone_D0_img6_ch2_all_3d']; i = i+1; %30
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I6\Ch3\Histone_D0_img6_ch3_all_3d']; i = i+1; %31
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I6\Ch4\Histone_D0_img6_ch4_all_3d']; i = i+1; %32
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I3\Ch2\Histone_D2_img3_ch2_all_3d']; i = i+1; %33
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I3\Ch3\Histone_D2_img3_ch3_all_3d']; i = i+1; %34
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I3\Ch4\Histone_D2_img3_ch4_all_3d']; i = i+1; %35
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I8\Ch2\Histone_D2_img8_ch2_all_3d']; i = i+1; %36
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I8\Ch3\Histone_D2_img8_ch3_all_3d']; i = i+1; %37
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I8\Ch4\Histone_D2_img8_ch4_all_3d']; i = i+1; %38

% ========================== Test One ==========================
% 
% testidx = 17;
% 
% mypath = img_paths{testidx,1};
% path_spotsrun = [mypath '_rnaspotsrun.mat'];
% fprintf("Examining %s...\n", mypath);
% 
% if isfile(path_spotsrun)
%     spotsrun = RNASpotsRun.loadFrom(path_spotsrun);
%     [spotsrun, img_filter] = spotsrun.loadFilteredImage();
%     
%     %Print some info
%     fprintf("Probe: %s\n", spotsrun.type_probe);
%     fprintf("Target Type: %s\n", spotsrun.type_targetmol);
%     fprintf("Species: %s\n", spotsrun.type_species);
%     
%     %Get maximum intensity.
%     fimg_max_val = max(img_filter, [], 'all', 'omitnan');
%     fprintf("Global maximum intensity value: %d\n", fimg_max_val);
%     
%     %Gen "histogram" (to be displayed as line graphs - one log scaled) w/
%     %   bin size of 1
%     ivals = [0:1:fimg_max_val];
%     valbin = histcounts(img_filter, fimg_max_val+1);
%     logvalbin = log10(valbin);
%     
%     fig_norm = figure(1);
%     plot(ivals,valbin);
%     
%     fig_log = figure(2);
%     plot(ivals,logvalbin);
%     
%     img_dbl = double(img_filter);
%     
%     percraw = prctile(img_dbl, [90 99], 'all');
%     perc90 = percraw(1);
%     perc99 = percraw(2);
%     
%     %Remove lowest 75%
%     img_dbl(img_dbl < perc99) = NaN;
%     perc99 = prctile(img_dbl, 99, 'all');
%     perc99 = round(perc99);
%     
%     fprintf("Auto range: %d - %d\n", perc90, perc99);
%     
% end

% ========================== Test All ==========================

path_count = i-1;

for j = 1:path_count
    mypath = img_paths{j,1};
    path_spotsrun = [mypath '_rnaspotsrun.mat'];
    fprintf("Examining %s...\n", mypath);
    
    if isfile(path_spotsrun)
        spotsrun = RNASpotsRun.loadFrom(path_spotsrun);
        [spotsrun, img_filter] = spotsrun.loadFilteredImage();
    
        %Print some info
        fprintf("Probe: %s\n", spotsrun.type_probe);
        fprintf("Target Type: %s\n", spotsrun.type_targetmol);
        fprintf("Species: %s\n", spotsrun.type_species);
    
        %Get maximum intensity.
        fimg_max_val = max(img_filter, [], 'all', 'omitnan');
        valbin = histcounts(img_filter, fimg_max_val+1);
        fprintf("Global maximum intensity value: %d\n", fimg_max_val);
        
        img_dbl = double(img_filter);
        img_dbl(img_dbl < 1) = NaN;
        img_std = std(img_dbl,0,'all', 'omitnan');
        fprintf("Standard Deviation: %f\n", img_std);
        
        percraw = prctile(img_dbl, [90 99], 'all');
        perc90 = percraw(1);
        perc99 = percraw(2);
        fprintf("90th percentile: %d\n", perc90);
        fprintf("99th percentile: %d\n", perc99);
        
        img_dbl(img_dbl < perc90) = NaN;
        mean90 = mean(img_dbl,'all', 'omitnan');
        stdev90 = std(img_dbl,0,'all', 'omitnan');
        med90 = median(img_dbl,'all', 'omitnan');
        
        fprintf("Mean >90 perc: %f\n", mean90);
        fprintf("Stdev >90 perc: %f\n", stdev90);
        fprintf("Median >90 perc: %d\n", med90);
        
        img_dbl(img_dbl < perc99) = NaN;
        mean99 = mean(img_dbl,'all', 'omitnan');
        stdev99 = std(img_dbl,0,'all', 'omitnan');
        med99 = median(img_dbl,'all', 'omitnan');
        
        fprintf("Mean >99 perc: %f\n", mean99);
        fprintf("Stdev >99 perc: %f\n", stdev99);
        fprintf("Median >99 perc: %d\n", med99);
        
%         img_dbl(img_dbl < perc99) = NaN;
%         perc99 = prctile(img_dbl, 99, 'all');
%         perc99 = round(perc99);
%         
%         fprintf("Auto range: %d - %d\n", perc90, perc99);
    end
end
