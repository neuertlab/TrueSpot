%

%Updated to use https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack

%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

img_paths = cell(48,1);
i = 1;

%----- Test
img_paths{i,1} = [ImgDir '\img\sim\mESC_histone_100x_1']; i = i+1;
img_paths{i,1} = [ImgDir '\img\sim\mESC_histone_100x_2']; i = i+1;
img_paths{i,1} = [ImgDir '\img\sim\mESC_histone_100x_3']; i = i+1;
img_paths{i,1} = [ImgDir '\img\sim\mESC_RNA_20x_1']; i = i+1;
img_paths{i,1} = [ImgDir '\img\sim\mESC_RNA_20x_2']; i = i+1;
img_paths{i,1} = [ImgDir '\img\sim\mESC_RNA_20x_3']; i = i+1;
img_paths{i,1} = [ImgDir '\img\sim\mESC_RNA_100x_1']; i = i+1;

% ========================== Read from dir ==========================

dir_path = [ImgDir '\img\histones_febc\mat'];
dir_contents = dir(dir_path);

% ========================== Cycle ==========================

addpath('./core');
%path_count = i-1;
path_count = size(dir_contents,1);

for j = 1:path_count
    %mypath = img_paths{j,1};
    fname = dir_contents(j).name;
    if dir_contents(j).isdir; continue; end
    if ~endsWith(fname, '.mat'); continue; end
    if ~contains(fname, 'H3K4me2'); continue; end
    [~, fname, ~] = fileparts(fname);
    mypath = [dir_path filesep fname];
    
    fprintf("MAT Found: %s\n", fname);
    %load([mypath '.mat'], 'imgdat');
    
    %Ben's images
    load([mypath '.mat'], 'CY5_ims', 'DAPI_ims', 'TMR_ims', 'TRANS_ims', 'YFP_ims');
%     X = size(TRANS_ims,2);
%     Y = size(TRANS_ims,1);
%     Z = size(TRANS_ims,3);
%     imgdat = uint16(zeros(Y,X,Z,3));
%     imgdat(:,:,:,1) = DAPI_ims(:,:,:); clear DAPI_ims;
%     imgdat(:,:,:,2) = CY5_ims(:,:,:); clear CY5_ims;
%     imgdat(:,:,:,3) = TMR_ims(:,:,:); clear TMR_ims;
%     imgdat(:,:,:,4) = YFP_ims(:,:,:); clear YFP_ims;
%     imgdat(:,:,:,5) = TRANS_ims(:,:,:); clear TRANS_ims;
    
    %Mtx2Tiff(imgdat, [mypath '.tif']);
    %options.big = true;
    %saveastiff(imgdat, [mypath '.tif']);
    saveastiff(DAPI_ims, [mypath '_C1.tif']); clear DAPI_ims;
    saveastiff(CY5_ims, [mypath '_C2.tif']); clear CY5_ims;
    saveastiff(TMR_ims, [mypath '_C3.tif']); clear TMR_ims;
    saveastiff(YFP_ims, [mypath '_C4.tif']); clear YFP_ims;
    saveastiff(TRANS_ims, [mypath '_C5.tif']); clear TRANS_ims;
end