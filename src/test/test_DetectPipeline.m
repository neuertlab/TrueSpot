%

%%
%base_dir = 'D:\usr\bghos\labdat\imgproc';
%base_dir = 'C:\Users\Blythe\labdata\imgproc';
base_dir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

tif_path = [base_dir '\img\yeast\Exp2_rep2_10min_im5.tif'];
cellseg_path = [base_dir '\data\cell_seg\yeast\Lab_Exp2_rep2_10min_im5.mat'];
light_channel = 4;
channel_count = 4;

img_name = 'E2R2-IM5-CH2';
out_dir = [base_dir '\data\preprocess\YeastFISH\E2R2\Img5\Ch2'];
rna_channel = 2;

%ctrl_path = [base_dir '\img\mESC_4d\20180202_4d_mESC_SCR-CY5_img_1_MMStack.ome.tif'];
%ctrl_ch = 2;
%ctrl_ch_count = 5;

%img_name = 'Histone_D2_img8_ch4';
%out_dir = [base_dir '\data\preprocess\histones\D2_I8\Ch4'];
%rna_channel = 4;

%tif_path = [base_dir '\img\histones\20201118_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_4_MMStack.ome.tif'];
%cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201118_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_4.mat'];
%light_channel = 5;
%channel_count = 5;

%img_name = 'Histone_D0_img4_ch4';
%out_dir = [base_dir '\data\preprocess\histones\D0_I4\Ch4'];
%rna_channel = 4;

%tif_path = [base_dir '\img\histones\20201117_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K4me2-AF488_img_6_MMStack.ome.tif'];
%cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201117_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K4me2-AF488_img_6.mat'];
%light_channel = 5;
%channel_count = 5;

%img_name = 'Histone_D0_img6_ch2';
%out_dir = [base_dir '\data\preprocess\histones\D0_I6\Ch2'];
%rna_channel = 2;

%tif_path = [base_dir '\img\histones\20201120_2d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K4me2-AF488_img_3_MMStack.ome.tif'];
%cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201120_2d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K4me2-AF488_img_3.mat'];
%light_channel = 5;
%channel_count = 5;

%img_name = 'Histone_D2_img3_ch4';
%out_dir = [base_dir '\data\preprocess\histones\D2_I3\Ch4'];
%rna_channel = 4;

%tif_path = [base_dir '\img\histones\20201120_2d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_8_MMStack.ome.tif'];
%cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201120_2d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_8.mat'];
%light_channel = 5;
%channel_count = 5;

%img_name = 'Histone_D2_img8_ch4';
%out_dir = [base_dir '\data\preprocess\histones\D2_I8\Ch4'];
%rna_channel = 4;

ctrl_path = [];
ctrl_ch = 0;
ctrl_ch_count = 0;

t_min = 1; 
t_max = 300; 
z_trim = 3;

ttune_winsize = 10;
ttune_winthresh = 0.9;

addpath('..');
Main_RNASpots(img_name, tif_path, rna_channel, light_channel, channel_count, out_dir, t_min, t_max, z_trim, cellseg_path, ctrl_path, ctrl_ch, ctrl_ch_count, ttune_winsize, ttune_winthresh, false);