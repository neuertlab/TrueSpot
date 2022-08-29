%
%%

addpath('./core');
% ========================== Input Paths ==========================

%% !! Replace with your own tif paths and info!
ch_count = 5;
ch_smpl = 2;
tif_path = 'C:\Users\hospelb\labdata\imgproc\img\histones_feb\20200205_2d_F1-2-1_H3K4me2-AF488-XistInt-CY5-Tsix5Int-TMR_img_6_MMStack.ome.tif';

% ========================== Output Paths ==========================

%% !! Replace with your own desired paths!
base_out_dir = 'D:\Users\hospelb\labdata\imgproc\imgproc\data\refsets';

out_stem = [base_out_dir filesep 'D2I6\Xist'];

% ========================== Generate New ==========================

selector = RNA_Threshold_SpotSelector.createEmptyRefSelector(tif_path, ch_count, ch_smpl, out_stem);
selector.launchRefSelectGUI();

% ========================== Open Saved ==========================

% selector = RNA_Threshold_SpotSelector.openSelector(out_stem);
% selector.launchRefSelectGUI();