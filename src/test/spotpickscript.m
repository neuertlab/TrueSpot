%
%%

%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc\data\preprocess';
% ========================== Image Channels ==========================
%De-comment one at a time

%----- mESC Set 1
%save_stem_rna = [ImgDir '\mESC4d\Tsix-AF594\Tsix-AF594_IMG1_all_3d'];
%save_stem_rna = [ImgDir '\mESC4d\Xist-CY5\Xist-CY5_IMG1_all_3d'];

%----- mESC Set 2 [NOT UPDATED]
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Tsix\mESC_1d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\1Day\Xist\mESC_1d_Xist_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\2Day\Tsix\mESC_2d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\2Day\Xist\mESC_2d_Xist_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\3Day\Tsix\mESC_3d_Tsix_all_3d'];
%save_stem_rna = [ImgDir '\data\preprocess\feb2019\3Day\Xist\mESC_3d_Xist_all_3d'];

%----- yeast RNA
%save_stem_rna = [ImgDir '\yeastrna\E1R2I4\STL1\E1R2-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\yeastrna\E1R2I4\CTT1\E1R2-CTT1-CY5_all_3d'];
%save_stem_rna = [ImgDir '\yeastrna\E2R1I2\STL1\E2R1-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\yeastrna\E2R1I2\CTT1\E2R1-CTT1-CY5_all_3d'];
%save_stem_rna = [ImgDir '\yeastrna\E2R2I3\STL1\E2R2I3-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\yeastrna\E2R2I3\CTT1\E2R2I3-CTT1-CY5_all_3d'];
%save_stem_rna = [ImgDir '\yeastrna\E2R2I5\STL1\E2R2I5-STL1-TMR_all_3d'];
%save_stem_rna = [ImgDir '\yeastrna\E2R2I5\CTT1\E2R2I5-CTT1-CY5_all_3d'];

%----- yeast proteins
%save_stem_rna = [ImgDir '\yeast_protein\Msb2\2M1m_img2\Msb2_02M_1m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\yeast_protein\Msb2\2M5m_img3\Msb2_02M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\yeast_protein\Msb2\4M1m_img3\Msb2_04M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\yeast_protein\Msb2\4M5m_img2\Msb2_04M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\yeast_protein\Msb2\NoSalt_img1\Msb2_NoSalt_img1_GFP_all_3d'];

%save_stem_rna = [ImgDir '\yeast_protein\Opy2\2M1m_img3\Opy2_02M_1m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\yeast_protein\Opy2\2M5m_img2\Opy2_02M_5m_img2_GFP_all_3d'];
%save_stem_rna = [ImgDir '\yeast_protein\Opy2\4M1m_img1\Opy2_04M_1m_img1_GFP_all_3d'];
%save_stem_rna = [ImgDir '\yeast_protein\Opy2\4M5m_img3\Opy2_04M_5m_img3_GFP_all_3d'];
%save_stem_rna = [ImgDir '\yeast_protein\Opy2\NoSalt_img2\Opy2_NoSalt_img2_GFP_all_3d'];

%----- Xist/Tsix + histones (FEB 2020)
%save_stem_rna = [ImgDir '\histones_feb2020\D0I1\Xist-CY5\20200218_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_1_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I1\Tsix-TMR\20200218_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_1_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I1\H3K4me2-AF488\20200218_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_1_MMStack-H3K4me2-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I2\Xist-CY5\20200214_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_2_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I2\Tsix-TMR\20200214_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_2_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I2\H3K4me2-AF488\20200214_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_2_MMStack-H3K4me2-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I4\Xist-CY5\20200207_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_4_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I4\Tsix-TMR\20200207_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_4_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I4\H3K4me2-AF488\20200207_0d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_4_MMStack-H3K4me2-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I9A\Xist-CY5\20200207_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_9_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I9A\Tsix-TMR\20200207_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_9_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I9A\H3K36me3-AF488\20200207_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_9_MMStack-H3K36me3-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I9B\Xist-CY5\20200221_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_9_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I9B\Tsix-TMR\20200221_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_9_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I9B\H3K36me3-AF488\20200221_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_9_MMStack-H3K36me3-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I10\Xist-CY5\20200214_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_10_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I10\Tsix-TMR\20200214_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_10_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D0I10\H3K36me3-AF488\20200214_0d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_10_MMStack-H3K36me3-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I3A\Xist-CY5\20200205_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_3_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I3A\Tsix-TMR\20200205_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_3_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I3A\H3K36me3-AF488\20200205_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_3_MMStack-H3K36me3-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I3B\Xist-CY5\20200221_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_3_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I3B\Tsix-TMR\20200221_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_3_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I3B\H3K36me3-AF488\20200221_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_3_MMStack-H3K36me3-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I5A\Xist-CY5\20200214_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_5_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I5A\Tsix-TMR\20200214_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_5_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I5A\H3K36me3-AF488\20200214_2d_F1-2-1_H3K36me3-AF488-XistInt-CY5-Tsix5Int-TMR_img_5_MMStack-H3K36me3-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I5B\Xist-CY5\20200218_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_5_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I5B\Tsix-TMR\20200218_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_5_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I5B\H3K4me2-AF488\20200218_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_5_MMStack-H3K4me2-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I6\Xist-CY5\20200205_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_6_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I6\Tsix-TMR\20200205_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_6_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I6\H3K4me2-AF488\20200205_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_6_MMStack-H3K4me2-AF488_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I7\Xist-CY5\20200214_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_7_MMStack-Xist-CY5_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I7\Tsix-TMR\20200214_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_7_MMStack-Tsix-TMR_all_3d'];
%save_stem_rna = [ImgDir '\histones_feb2020\D2I7\H3K4me2-AF488\20200214_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_7_MMStack-H3K4me2-AF488_all_3d'];

%----- Xist/Tsix + histones (NOV 2020)
%save_stem_rna = [ImgDir '\histones_nov2020\D0I4\Xist\Histone_D0_img4_ch2_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D0I4\Tsix\Histone_D0_img4_ch3_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D0I4\H3K36me3\Histone_D0_img4_ch4_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D0I6\Xist\Histone_D0_img6_ch2_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D0I6\Tsix\Histone_D0_img6_ch3_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D0I6\H3K4me2\Histone_D0_img6_ch4_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D2I3\Xist\Histone_D2_img3_ch2_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D2I3\Tsix\Histone_D2_img3_ch3_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D2I3\H3K4me2\Histone_D2_img3_ch4_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D2I8\Xist\Histone_D2_img8_ch2_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D2I8\Tsix\Histone_D2_img8_ch3_all_3d'];
%save_stem_rna = [ImgDir '\histones_nov2020\D2I8\H3K36me3\Histone_D2_img8_ch4_all_3d'];

% ========================== Load Run Info ==========================

addpath('./core');
spotsrun = RNASpotsRun.loadFrom(save_stem_rna);
spotsrun.out_stem = save_stem_rna;
spotsrun = spotsrun.saveMe();

% ========================== Generation ==========================
%(This is for my use, you don't need to worry about this section if spot
%picking.)

%-> AnnoObj already exists
%selector = RNA_Threshold_SpotSelector.openSelector(save_stem_rna);

%-> New AnnoObj
% selector = RNA_Threshold_SpotSelector;
% selector = selector.initializeNew(save_stem_rna, spotsrun.intensity_threshold - spotsrun.t_min + 1);
% selector.ztrim = 7;
% selector.selmcoords = zeros(4,1);
% selector.selmcoords(1,1) = uint16(504);
% selector.selmcoords(2,1) = uint16(1528);
% selector.selmcoords(3,1) = uint16(462);
% selector.selmcoords(4,1) = uint16(1486);
% selector = selector.launchRefSelectGUI();

%selector = selector.saveMe();

% ========================== Spawn AnnoObj ==========================

selector = RNA_Threshold_SpotSelector.openSelector(save_stem_rna, true);
selector = selector.launchRefSelectGUI();

