%
%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

%----- Test
data_file = [ImgDir '\img\sim\yeast_proteinGFP_100x_1.mat'];

% ========================== Run ==========================

addpath('./core');
load(data_file, 'imgdat');

%Extra Gaussian?
%imgdat = RNA_Threshold_Common.applyGaussianFilter(imgdat, 7, 2, 2);
%imgdat = RNA_Threshold_Common.applyGaussianFilter(imgdat, 5, 3, 2);

fig_handle = MatImages.viewMaxProjection(imgdat, true, 505);