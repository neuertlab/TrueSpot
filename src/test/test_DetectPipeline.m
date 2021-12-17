%

%%
%base_dir = 'D:\usr\bghos\labdat\imgproc';
%base_dir = 'C:\Users\Blythe\labdata\imgproc';
base_dir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

tif_path = [base_dir '\img\histones\20201118_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_4_MMStack.ome.tif'];
img_name = 'Histone_D0_img4_ch3';
out_dir = [base_dir '\data\preprocess\histones\D0_I4\Ch3'];
cellseg_path = [base_dir '\data\cell_seg\histones\Lab_20201118_0d_F1-2-1_XistInt-CY5-Tsix5Int-TMR-H3K36me3-AF488_img_4.mat'];
rna_channel = 3;
light_channel = 5;
channel_count = 5;

ctrl_path = [];
ctrl_ch = 0;
ctrl_ch_count = 0;

t_min = 10; 
t_max = 400; 
z_trim = 5;

addpath('..');
Main_RNASpots(img_name, tif_path, rna_channel, light_channel, channel_count, out_dir, t_min, t_max, z_trim, cellseg_path, ctrl_path, ctrl_ch, ctrl_ch_count, -1, -1, false);