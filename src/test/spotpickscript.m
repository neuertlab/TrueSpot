%
%%

%%  !! UPDATE TO YOUR BASE DIR
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
save_stem_rna = [ImgDir '\data\preprocess\YeastFISH\E2R2\Img5\Ch2\E2R2-IM5-CH2_all_3d'];

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
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0I4\Ch2\Histone_D0_img4_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0I4\Ch3\Histone_D0_img4_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0I4\Ch4\Histone_D0_img4_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0I6\Ch2\Histone_D0_img6_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0I6\Ch3\Histone_D0_img6_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D0I6\Ch4\Histone_D0_img6_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2I3\Ch2\Histone_D2_img3_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2I3\Ch3\Histone_D2_img3_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2I3\Ch4\Histone_D2_img3_ch4_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2I8\Ch2\Histone_D2_img8_ch2_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2I8\Ch3\Histone_D2_img8_ch3_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\histones\D2I8\Ch4\Histone_D2_img8_ch4_all_3d'];

% ========================== Load Run Info ==========================

addpath('./core');
spotsrun = RNASpotsRun.loadFrom(save_stem_rna);

% ========================== Generation ==========================
%(This is for my use, you don't need to worry about this section if spot
%picking.)

%-> AnnoObj already exists
%selector = RNA_Threshold_SpotSelector.openSelector(save_stem_rna);

%-> New AnnoObj
selector = RNA_Threshold_SpotSelector;
selector = selector.initializeNew(save_stem_rna, spotsrun.intensity_threshold - spotsrun.t_min + 1);
%selector = selector.launchRefSelectGUI();
selector.selmcoords = zeros(4,1);
selector.selmcoords(1,1) = uint16(702);
selector.selmcoords(2,1) = uint16(1726);
selector.selmcoords(3,1) = uint16(198);
selector.selmcoords(4,1) = uint16(1222);

selector = selector.saveMe();

% ========================== Spawn AnnoObj ==========================

%selector = RNA_Threshold_SpotSelector.openSelectorSetPaths(save_stem_rna);
%selector = selector.launchRefSelectGUI();

