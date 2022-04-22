%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Image Channels ==========================

img_paths = cell(48,1);
i = 1;

%----- mESC Set 1
img_paths{i,1} = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\feb2018\Xist_CY5\Xist\Xist-CY5_IMG1_all_3d']; i = i+1;

%----- mESC Set 2
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\1Day\Tsix\mESC_1d_Tsix_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\1Day\Xist\mESC_1d_Xist_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\2Day\Tsix\mESC_2d_Tsix_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\2Day\Xist\mESC_2d_Xist_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\3Day\Tsix\mESC_3d_Tsix_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\feb2019\3Day\Xist\mESC_3d_Xist_all_3d']; i = i+1;

%----- yeast RNA
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E1R2\STL1\E1R2-STL1-TMR_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E1R2\CTT1\E1R2-CTT1-CY5_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R1\STL1\E2R1-STL1-TMR_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R1\CTT1\E2R1-CTT1-CY5_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\STL1\E2R2I3-STL1-TMR_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img3\CTT1\E2R2I3-CTT1-CY5_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\STL1\E2R2I5-STL1-TMR_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\CTT1\E2R2I5-CTT1-CY5_all_3d']; i = i+1;

%----- yeast proteins
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\4M1m_img3\Msb2_04M_1m_img3_GFP_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\4M5m_img2\Msb2_04M_5m_img2_GFP_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\msb2\NoSalt_img1\Msb2_NoSalt_img1_GFP_all_3d']; i = i+1;

img_paths{i,1} = [ImgDir '\data\preprocess\opy2\2M1m_img3\Opy2_02M_1m_img3_GFP_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\opy2\2M5m_img2\Opy2_02M_5m_img2_GFP_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\opy2\4M1m_img1\Opy2_04M_1m_img1_GFP_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\opy2\4M5m_img3\Opy2_04M_5m_img3_GFP_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\opy2\NoSalt_img2\Opy2_NoSalt_img2_GFP_all_3d']; i = i+1;

%----- Xist/Tsix + histones
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I4\Ch2\Histone_D0_img4_ch2_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I4\Ch3\Histone_D0_img4_ch3_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I4\Ch4\Histone_D0_img4_ch4_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I6\Ch2\Histone_D0_img6_ch2_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I6\Ch3\Histone_D0_img6_ch3_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D0_I6\Ch4\Histone_D0_img6_ch4_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I3\Ch2\Histone_D2_img3_ch2_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I3\Ch3\Histone_D2_img3_ch3_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I3\Ch4\Histone_D2_img3_ch4_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I8\Ch2\Histone_D2_img8_ch2_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I8\Ch3\Histone_D2_img8_ch3_all_3d']; i = i+1;
img_paths{i,1} = [ImgDir '\data\preprocess\histones\D2_I8\Ch4\Histone_D2_img8_ch4_all_3d']; i = i+1;

% ========================== Cycle ==========================

%addpath('..');
%addpath('../core');
addpath('./core');
path_count = i-1;

for j = 1:path_count
    mypath = img_paths{j,1};
    path_spotsrun = [mypath '_rnaspotsrun.mat'];
    fprintf("Now processing %s... (%d of %d)\n", mypath, j, path_count);
    if isfile(path_spotsrun)
        spotsrun = RNASpotsRun.loadFrom(path_spotsrun);
        spotsrun.ttune_winsize = 15;
        RNA_Pipeline_Core(spotsrun, 1, []);
    else
        fprintf("WARNING: Could not find run data for %s. Skipping...\n", mypath);
    end
end


