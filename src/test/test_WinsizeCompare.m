%
%%

ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Image Channels ==========================
%De-comment one at a time

%----- mESC Set 1
%save_stem_rna = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2018\Xist_CY5\Xist\Xist-CY5_IMG1_all_3d'];

%----- mESC Set 2
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Tsix\mESC_1d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Xist\mESC_1d_Xist_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\2Day\Tsix\mESC_2d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\2Day\Xist\mESC_2d_Xist_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\3Day\Tsix\mESC_3d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\3Day\Xist\mESC_3d_Xist_all_3d'];

%----- yeast RNA
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E1R2\Ch1\E1R2-CH1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E1R2\Ch2\E1R2-CH2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R1\Ch1\E2R1-CH1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R1\Ch2\E2R1-CH2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\Ch1\E2R2-IM3-CH1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\Ch2\E2R2-IM3-CH2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\Ch1\E2R2-IM5-CH1_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\Ch2\E2R2-IM5-CH2_all_3d'];

%----- yeast proteins
%save_stem_rna = [ImgDir '\data\preprocess\msb2\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M1m_img3\Msb2_04M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\4M5m_img2\Msb2_04M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\msb2\NoSalt_img1\Msb2_NoSalt_img1_GFP_all_3d'];

%save_stem_rna = [ImgDir '\data\preprocess\opy2\2M1m_img3\Opy2_02M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\2M5m_img2\Opy2_02M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\4M1m_img1\Opy2_04M_1m_img1_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\4M5m_img3\Opy2_04M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\opy2\NoSalt_img2\Opy2_NoSalt_img2_GFP_all_3d'];

%----- Xist/Tsix + histones
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch2\Histone_D0_img4_ch2_all_3d'];
save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch3\Histone_D0_img4_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I4\Ch4\Histone_D0_img4_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch2\Histone_D0_img6_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch3\Histone_D0_img6_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0_I6\Ch4\Histone_D0_img6_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I3\Ch2\Histone_D2_img3_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I3\Ch3\Histone_D2_img3_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I3\Ch4\Histone_D2_img3_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I8\Ch2\Histone_D2_img8_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I8\Ch3\Histone_D2_img8_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2_I8\Ch4\Histone_D2_img8_ch4_all_3d'];

% ========================== Load Run Info ==========================

addpath('./core');
spotsrun = RNASpotsRun.loadFrom(save_stem_rna);

% ========================== Gen Graphs? ==========================

%Load spot count tables
if spotsrun.ztrim > spotsrun.ztrim_auto
    load([spotsrun.out_stem '_ztrim' num2str(spotsrun.ztrim)], 'trimmed_coords');
    [~, spots_sample] = spotsrun.loadSpotsTable();
    T = size(trimmed_coords,1);
    for i = 1:T
        spots_sample(i,2) = size(trimmed_coords{i},1);
    end
    
    load([spotsrun.ctrl_stem '_ztrim' num2str(spotsrun.ztrim)], 'trimmed_coords');
    [~, spots_control] = spotsrun.loadControlSpotsTable();
    T = size(trimmed_coords,1);
    for i = 1:T
        spots_control(i,2) = size(trimmed_coords{i},1);
    end
else
    [~, spots_sample] = spotsrun.loadSpotsTable();
    [~, spots_control] = spotsrun.loadControlSpotsTable();
end

%Try to load fscore table
fscores = [];
fscore_path = [save_stem_rna '_fscores.csv'];
Tmax = 0;
if isfile(fscore_path)
    fscores = csvread(fscore_path);
    Tmax = size(fscores,1);
end

%Determine window sizes to try
win_sizes = [5, 10, 15, 20, 25];
wincount = size(win_sizes,2);
winscore_plots = cell(wincount,1);
picked_threshes = NaN(wincount,2); %Col 1 is intensity thresh, col 2 is winscore thresh

%MAD factors to try
mad_factors = [-2.0, -1.5, -1.0, -0.75,-0.5, -0.25, 0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.0, 1.5, 2.0];
mfcount = size(mad_factors,2);
th_mtx = zeros(mfcount, wincount);
fs_mtx = zeros(mfcount, wincount);

%Generate winscore plots.
for i = 1:wincount
    for j = 1:mfcount
        [thresh, win_scores, score_thresh, ~] = RNA_Threshold_Common.estimateThreshold(spots_sample, spots_control, win_sizes(1,i), 0.0, mad_factors(1,j));
        if mad_factors(1,j) == 0.0
            picked_threshes(i,1) = thresh;
            picked_threshes(i,2) = score_thresh;
            winscore_plots{i,1} = win_scores;
        end
        th_mtx(j,i) = thresh;
        if ~isempty(fscores) && thresh > 0 && thresh <= Tmax
            fs_mtx(j,i) = fscores(thresh,1);
        end
    end
end
