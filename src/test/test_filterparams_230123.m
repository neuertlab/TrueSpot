%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Settings ==========================
addpath('./core');

ImgName = 'rsfish_sim_bgs_30A';
gaussian_rad = 7;

% ========================== Read Table ==========================

InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

rec_row = 0;
rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    if strcmp(iname, ImgName)
        rec_row = r;
        break;
    end
end

if rec_row < 1
    fprintf('Image with name %s could not be found!\n', ImgName);
    return;
end

% ========================== Run ==========================

tif_path = [ImgDir replace(getTableValue(image_table, rec_row, 'IMAGEPATH'), '/', filesep)];
rna_ch = image_table{rec_row, 'CHANNEL'};
total_ch = image_table{rec_row, 'CH_TOTAL'};

[tif_channels, ~] = LoadTif(tif_path, total_ch, [rna_ch], 1);
my_channel = tif_channels{rna_ch,1};
clear tif_channels;

my_channel = uint16(my_channel);
imin = min(my_channel, [], 'all');
imax = max(my_channel, [], 'all');

%Linear rescale??

%Max proj draw
max_proj = double(max(my_channel,[],3));
Lmin = min(max_proj(:));
Lmax = median(max_proj(:)) + round(10 * std(max_proj(:)));

figure(1);
imshow(max_proj, [Lmin Lmax]);

%Apply filter
IMG_filtered = RNA_Threshold_Common.applyGaussianFilter(my_channel, gaussian_rad, 2);
IMG_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(IMG_filtered);
IMG_filtered = RNA_Threshold_Common.blackoutBorders(IMG_filtered, gaussian_rad+1, 0);

%Rescale the filtered image if not enough range
ifmin = min(IMG_filtered, [], 'all');
ifmax = max(IMG_filtered, [], 'all');
ifrng = ifmax - ifmin;
if ifrng < 25
    %Linear rescale to 0-255
    IMG_filtered = ((IMG_filtered - ifmin) .* 255) ./ ifrng;
end

IMG_filtered = uint16(IMG_filtered);

ifmin = min(IMG_filtered, [], 'all');
ifmax = max(IMG_filtered, [], 'all');
max_proj_f = double(max(IMG_filtered,[],3));
LminF = median(max_proj_f(:)) - round(0 * std(max_proj_f(:)));
LmaxF = median(max_proj_f(:)) + round(10 * std(max_proj_f(:)));

figure(2);
imshow(max_proj_f, [LminF LmaxF]);

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end