% Example/test script for running RNA spot detection
%

%%
%Parameters
%base_dir = 'D:\usr\bghos\labdat\imgproc';
base_dir = 'C:\Users\Blythe\labdata\imgproc';

%tif_path = [base_dir '\img\mESC_4d\20180202_4d_mESC_Tsix-AF594_img_1_MMStack.ome.tif'];
%img_name = 'Tsix-AF594';
%out_dir = [base_dir '\data\preprocess\feb2018\Tsix_AF594\Tsix'];
%rna_channel = 3;
%channel_count = 5;
%t_min = 1; %Min threshold to try
%t_max = 200; %Max threshold to try

tif_path = [base_dir '\img\mESC_4d\20180202_4d_mESC_SCR-AF594_img_1_MMStack.ome.tif'];
img_name = 'SCR-AF594';
out_dir = [base_dir '\data\preprocess\feb2018\Tsix_AF594\SCR'];
rna_channel = 3;
channel_count = 5;
t_min = 1; %Min threshold to try
t_max = 200; %Max threshold to try

%tif_path = [base_dir
%'\img\mESC_4d\20180202_4d_mESC_No-probe_img_2_MMStack.ome.tif'];
%                                                                                                                    
%img_name = 'NoProbe-AF594';
%out_dir = [base_dir '\data\preprocess\feb2018\NoProbe\AF594'];
%rna_channel = 3;
%channel_count = 5;
%t_min = 1; %Min threshold to try
%t_max = 200; %Max threshold to try

%tif_path = [base_dir '\img\mESC_4d\20180202_4d_mESC_No-probe_img_2_MMStack.ome.tif'];
%img_name = 'NoProbe-CY5';
%out_dir = [base_dir '\data\preprocess\feb2018\NoProbe\CY5'];
%rna_channel = 2;
%channel_count = 5;
%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%tif_path = [base_dir '\img\mESC_4d\20180202_4d_mESC_SCR-CY5_img_1_MMStack.ome.tif'];
%img_name = 'SCR-CY5';
%out_dir = [base_dir '\data\preprocess\feb2018\Xist_CY5\SCR'];
%rna_channel = 2;
%channel_count = 5;
%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%tif_path = [base_dir '\img\yeast\Exp1_rep2_10min_im4.tif'];
%img_name = 'E1R2-CH2';
%out_dir = [base_dir '\data\preprocess\YeastFISH\E1R2\Ch2'];
%rna_channel = 2;
%channel_count = 4;
%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%tif_path = [base_dir '\img\yeast\Exp2_rep1_10min_im2.tif'];
%img_name = 'E2R1-CH2';
%out_dir = [base_dir '\data\preprocess\YeastFISH\E2R1\Ch2'];
%rna_channel = 2;
%channel_count = 4;
%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%tif_path = [base_dir '\img\yeast\Exp2_rep2_10min_im3.tif'];
%img_name = 'E2R2-IM3-CH2';
%out_dir = [base_dir '\data\preprocess\YeastFISH\E2R2\Img3\Ch2'];
%rna_channel = 2;
%channel_count = 4;
%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%tif_path = [base_dir '\img\yeast\Exp2_rep2_10min_im5.tif'];
%img_name = 'E2R2-IM5-CH2';
%out_dir = [base_dir '\data\preprocess\YeastFISH\E2R2\Img5\Ch2'];
%rna_channel = 2;
%channel_count = 4;
%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%tif_path = [base_dir '\img\mESC\20190218_1d_F1-2-1_Xist-CY5-Tsix-TMR_img_9_MMStack.ome.tif'];
%img_name = 'mESC_1d_Xist';
%img_name = 'mESC_1d_Tsix';
%out_dir = [base_dir '\data\preprocess\feb2019\1Day\Xist'];
%out_dir = [base_dir '\data\preprocess\feb2019\1Day\Tsix'];
%rna_channel = 2;
%rna_channel = 3;
%channel_count = 4;
%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%tif_path = [base_dir '\img\mESC\20190220_2d_F1-2-1_Xist-CY5-Tsix-TMR_img_9_MMStack.ome.tif'];
%img_name = 'mESC_2d_Xist';
%img_name = 'mESC_2d_Tsix';
%out_dir = [base_dir '\data\preprocess\feb2019\2Day\Xist'];
%out_dir = [base_dir '\data\preprocess\feb2019\2Day\Tsix'];
%rna_channel = 2;
%rna_channel = 3;
%channel_count = 4;
%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%tif_path = [base_dir '\img\mESC\20190221_3d_F1-2-1_Xist-CY5-Tsix-TMR_img_9_MMStack.ome.tif'];
%img_name = 'mESC_3d_Xist';
%img_name = 'mESC_3d_Tsix';
%out_dir = [base_dir '\data\preprocess\feb2019\3Day\Xist'];
%out_dir = [base_dir '\data\preprocess\feb2019\3Day\Tsix'];
%rna_channel = 2;
%rna_channel = 3;
%channel_count = 4;

%tif_path = [base_dir '\img\msb2\20180306_Msb2-GFP_02M_NaCl_1min_2_MMStack.ome.tif'];
%img_name = 'Msb2_02M_1m_img2_GFP';
%out_dir = [base_dir '\data\preprocess\msb2\2M1m_img2'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\msb2\20180306_Msb2-GFP_02M_NaCl_5min_3_MMStack.ome.tif'];
%img_name = 'Msb2_02M_5m_img3_GFP';
%out_dir = [base_dir '\data\preprocess\msb2\2M5m_img3'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\msb2\20180306_Msb2-GFP_04M_NaCl_1min_3_MMStack.ome.tif'];
%img_name = 'Msb2_04M_1m_img3_GFP';
%out_dir = [base_dir '\data\preprocess\msb2\4M1m_img3'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\msb2\20180306_Msb2-GFP_04M_NaCl_5min_2_MMStack.ome.tif'];
%img_name = 'Msb2_04M_5m_img2_GFP';
%out_dir = [base_dir '\data\preprocess\msb2\4M5m_img2'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\msb2\20180306_Msb2-GFP_no_salt_1_MMStack.ome.tif'];
%img_name = 'Msb2_NoSalt_img1_GFP';
%out_dir = [base_dir '\data\preprocess\msb2\NoSalt_img1'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\opy2\20180306_Opy2-GFP_02M_NaCl_1min_3_MMStack.ome.tif'];
%img_name = 'Opy2_02M_1m_img3_GFP';
%out_dir = [base_dir '\data\preprocess\opy2\2M1m_img3'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\opy2\20180306_Opy2-GFP_02M_NaCl_5min_2_MMStack.ome.tif'];
%img_name = 'Opy2_02M_5m_img2_GFP';
%out_dir = [base_dir '\data\preprocess\opy2\2M5m_img2'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\opy2\20180306_Opy2-GFP_04M_NaCl_1min_1_MMStack.ome.tif'];
%img_name = 'Opy2_04M_1m_img1_GFP';
%out_dir = [base_dir '\data\preprocess\opy2\4M1m_img1'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\opy2\20180306_Opy2-GFP_04M_NaCl_5min_3_MMStack.ome.tif'];
%img_name = 'Opy2_04M_5m_img3_GFP';
%out_dir = [base_dir '\data\preprocess\opy2\4M5m_img3'];
%rna_channel = 1;
%channel_count = 2;

%tif_path = [base_dir '\img\opy2\20180306_Opy2-GFP_no_salt_2_MMStack.ome.tif'];
%img_name = 'Opy2_NoSalt_img2_GFP';
%out_dir = [base_dir '\data\preprocess\opy2\NoSalt_img2'];
%rna_channel = 1;
%channel_count = 2;

%t_min = 1; %Min threshold to try
%t_max = 300; %Max threshold to try

%And call the main!
addpath('..');
Main_RNASpotDetect(img_name, tif_path, out_dir, rna_channel, channel_count, t_min, t_max, true);