%
%%

addpath('./core');
DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Constants ==========================

csv_dir = [DataDir '\data\rsfish\yeast_tc\E2R2\CH1\sctc_E2R2_0m_I1_STL1'];
out_stem = [DataDir '\data\rsfish\yeast_tc\E2R2\CH1\sctc_E2R2_0m_I1_STL1\RSFISH_E2R2_0m_I1_STL1'];

% ========================== Run ==========================

Main_RSFish2Mat(csv_dir, out_stem);