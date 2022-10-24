%
%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

img_paths = cell(48,1);
i = 1;

%----- Test
img_paths{i,1} = [ImgDir '\img\sim\mESC_RNA_100x_1']; i = i+1;

% ========================== Cycle ==========================

addpath('./core');
path_count = i-1;

for j = 1:path_count
    mypath = img_paths{j,1};
    load([mypath '.mat'], 'imgdat');
    Mtx2Tiff(imgdat, [mypath '.tif']);
end