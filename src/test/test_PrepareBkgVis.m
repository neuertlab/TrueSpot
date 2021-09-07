%
%%
%base_dir = 'D:\usr\bghos\labdat\imgproc';
base_dir = 'C:\Users\Blythe\labdata\imgproc';

%Tsix
%SpotDat_Path = [base_dir '\data\preprocess\feb2018\Tsix_AF594\Tsix\all_3d\Tsix-AF594_IMG1_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\mESC_4d\Tsix-AF594_img_1_bkgmask.mat'];

%Xist
%SpotDat_Path = [base_dir '\data\preprocess\feb2018\Xist_CY5\Xist'];
%Bkg_Path = [base_dir '\data\bkgmask\mESC_4d\Xist-CY5_img_1_bkgmask.mat'];

%Other mESC images
%SpotDat_Path = [base_dir '\data\preprocess\feb2019\1Day\Tsix\all_3d\mESC_1d_Tsix_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\feb2019\1Day\Xist\all_3d\mESC_1d_Xist_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\mESC\1d_F1-2-1_Xist-CY5-Tsix-TMR_img_9_bkgmask.mat'];
%SpotDat_Path = [base_dir '\data\preprocess\feb2019\2Day\Tsix\all_3d\mESC_2d_Tsix_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\feb2019\2Day\Xist\all_3d\mESC_2d_Xist_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\mESC\2d_F1-2-1_Xist-CY5-Tsix-TMR_img_9_bkgmask.mat'];
%SpotDat_Path = [base_dir '\data\preprocess\feb2019\3Day\Tsix\all_3d\mESC_3d_Tsix_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\feb2019\3Day\Xist\all_3d\mESC_3d_Xist_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\mESC\3d_F1-2-1_Xist-CY5-Tsix-TMR_img_9_bkgmask.mat'];

%Yeast
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E1R2\Ch1\all_3d\E1R2-CH1_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E1R2\Ch2\all_3d\E1R2-CH2_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\yeast\Exp1_rep2_10min_im4_bkgmask.mat'];
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E2R1\Ch1\all_3d\E2R1-CH1_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E2R1\Ch2\all_3d\E2R1-CH2_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\yeast\Exp2_rep1_10min_im2_bkgmask.mat'];
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E2R2\Img3\Ch1\all_3d\E2R2-IM3-CH1_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E2R2\Img3\Ch2\all_3d\E2R2-IM3-CH2_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\yeast\Exp2_rep2_10min_im3_bkgmask.mat'];
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E2R2\Img5\Ch1\all_3d\E2R2-IM5-CH1_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\YeastFISH\E2R2\Img5\Ch2\all_3d\E2R2-IM5-CH2_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\yeast\Exp2_rep2_10min_im5_bkgmask.mat'];

%Histones
%Tif_path = [base_dir '\img\histones\20201117_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K4me2-AF488_img_6_MMStack.ome.tif'];
%Cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201117_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K4me2-AF488_img_6.mat'];
%Bkg_Path = [base_dir '\data\bkgmask\histones\0d_F1-2-1_H3K4me2_img_6_bkgmask'];
%Tif_path = [base_dir '\img\histones\20201118_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_4_MMStack.ome.tif'];
%Cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201118_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_4.mat'];
%Bkg_Path = [base_dir '\data\bkgmask\histones\0d_F1-2-1_H3K36me3_img_4_bkgmask'];
%Tif_path = [base_dir '\img\histones\20201120_2d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K4me2-AF488_img_3_MMStack.ome.tif'];
%Cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201120_2d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K4me2-AF488_img_3.mat'];
%Bkg_Path = [base_dir '\data\bkgmask\histones\2d_F1-2-1_H3K4me2_img_3_bkgmask'];
%Tif_path = [base_dir '\img\histones\20201120_2d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_8_MMStack.ome.tif'];
%Cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201120_2d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_8.mat'];
%Bkg_Path = [base_dir '\data\bkgmask\histones\2d_F1-2-1_H3K36me3_img_8_bkgmask'];
%Main_BackgroundMask(Tif_path, Cellseg_path, Bkg_Path, '5', '5', true);
%return;
%printf('But it did not return...\n');

%SpotDat_Path = [base_dir '\data\preprocess\histones\D0_I6\Ch2\all_3d\Histone_D0_img6_ch2_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\histones\D0_I6\Ch3\all_3d\Histone_D0_img6_ch3_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\histones\D0_I6\Ch4\all_3d\Histone_D0_img6_ch4_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\histones\0d_F1-2-1_H3K4me2_img_6_bkgmask.mat'];

%SpotDat_Path = [base_dir '\data\preprocess\histones\D0_I4\Ch2\all_3d\Histone_D0_img4_ch2_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\histones\D0_I4\Ch3\all_3d\Histone_D0_img4_ch3_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\histones\D0_I4\Ch4\all_3d\Histone_D0_img4_ch4_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\histones\0d_F1-2-1_H3K36me3_img_4_bkgmask.mat'];

%SpotDat_Path = [base_dir '\data\preprocess\histones\D2_I3\Ch2\all_3d\Histone_D2_img3_ch2_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\histones\D2_I3\Ch3\all_3d\Histone_D2_img3_ch3_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\histones\D2_I3\Ch4\all_3d\Histone_D2_img3_ch4_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\histones\2d_F1-2-1_H3K4me2_img_3_bkgmask.mat'];

%SpotDat_Path = [base_dir '\data\preprocess\histones\D2_I8\Ch2\all_3d\Histone_D2_img8_ch2_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\histones\D2_I8\Ch3\all_3d\Histone_D2_img8_ch3_all_3d'];
%SpotDat_Path = [base_dir '\data\preprocess\histones\D2_I8\Ch4\all_3d\Histone_D2_img8_ch4_all_3d'];
%Bkg_Path = [base_dir '\data\bkgmask\histones\2d_F1-2-1_H3K36me3_img_8_bkgmask.mat'];

%Membrane Proteins

Bkg_CFilter_Path = [SpotDat_Path '_bkgmasked'];

%Do mask (if not already done)


%Generating bkg masked spot set (if not already done so)
addpath('./core');
load(Bkg_Path, 'background_mask');
tbl_path_RNA = [SpotDat_Path '_coordTable'];
load(tbl_path_RNA, 'coord_table');
[masked_spot_table, masked_coord_table] = RNA_Threshold_Common.mask_spots(background_mask, coord_table);
coord_table = masked_coord_table;
spot_table = masked_spot_table;
save([Bkg_CFilter_Path '_coordTable'], 'coord_table');
save([Bkg_CFilter_Path '_spotTable'], 'spot_table');
RNA_Threshold_Common.scale_mask_spot_table(Bkg_Path, Bkg_CFilter_Path);

%Generating masked image structs (if not already done so)
load([SpotDat_Path '_imgviewstructs'], 'my_images');
my_images(1).image = immultiply(my_images(1).image, background_mask);
my_images(2).image = immultiply(my_images(2).image, background_mask);

%Save
save([Bkg_CFilter_Path '_imgviewstructs'], 'my_images');