%
%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

outdir_root = [ImgDir '\data\preprocess\sim'];

%----- Test
data_file = [ImgDir '\img\sim\mESC_RNA_100x_1.mat'];

% ========================== Run ==========================

addpath('./core');

[~, fname, ~] = fileparts(data_file);

Main_RNASpots('-imgname', fname, '-matimg', data_file, '-outdir', [outdir_root filesep fname], '-thmax', 500, '-ztrim', 5, ...
    '-voxelsize', '(300, 65, 65)', '-debug');