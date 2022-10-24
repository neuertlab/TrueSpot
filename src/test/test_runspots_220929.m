%
%%  !! UPDATE TO YOUR BASE DIR
%BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
BaseDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Parameters ==========================

imgname = 'ROI004_XY1657822935_Z00_T0_C1';
tifpath = [BaseDir '\img\munsky_lab_test\ROI004_XY1657822935_Z00_T0_C1.tif'];
outdir = [BaseDir '\data\preprocess\munsky_lab\test\GFP'];
sample_ch = 1;
total_ch = 1;
light_ch = 0;
ztrim = 3;

% ========================== Call Main ==========================

addpath('./core');

rna_spot_run = Main_RNASpots('-imgname', imgname, '-tif', tifpath, '-outdir', outdir,...
    '-chsamp', sample_ch, '-chtrans', light_ch, '-chtotal', total_ch, '-ztrim', ztrim);