%
%%

base_dir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

out_stem = [base_dir '\data\preprocess\histones\D0_I4\Ch3\Histone_D0_img4_ch3_all_3d'];

spotsrun = RNASpotsRun.loadFrom(out_stem);
%[~, coord_table] = spotsrun.loadCoordinateTable();
%[~, background_mask] = spotsrun.loadBackgroundMask();
[~, spot_table] = spotsrun.loadSpotsTable();
[~, th_list] = spotsrun.loadThresholdTable();