%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Constants ==========================

SingleImgName = 'rsfish_sim_bgl_30A';
%OverrideName = 'Tsix-AF594_IMG1';
OverrideName = SingleImgName;

InputTablePath = [DataDir filesep 'test_images.csv'];

percentile_values = [50 75 80 85 90 95 99];
pcheck_count = size(percentile_values,2);
perc_index_99 = 7;
gaussian_rad = 7;

RAW_RESCALE_TRIGGER = 256;
RAW_RESCALE_TARGET = 512;

% ========================== Load csv Table ==========================

imgtbl = testutil_opentable(InputTablePath);

% ========================== Find Record ==========================
addpath('./core');

rec_row = 0;
rec_count = size(imgtbl,1);
for r = 1:rec_count
    myname = getTableValue(imgtbl, r, 'IMGNAME');
    if strcmp(myname, SingleImgName); rec_row = r; break; end
end

if rec_row < 1
    fprintf("Couldn't find image: %s!\n", SingleImgName);
    return;
end

% ========================== Do the thing ==========================

%Get tif path

tifpath_raw = getTableValue(imgtbl, r, 'IMAGEPATH');
tifpath = [ImgDir replace(tifpath_raw, '/', filesep)];
if endsWith(tifpath, '.mat')
    [tifdir, tifname, ~] = fileparts(tifpath);
    tifpath = [tifdir filesep 'tif' filesep tifname '.tif'];
end

if ~isfile(tifpath)
    fprintf('Tif file %s not found. Skipping...', tifpath);
    return;
end

%Load tif
ch_total = getTableValue(imgtbl, r, 'CH_TOTAL');
ch_sample = getTableValue(imgtbl, r, 'CHANNEL');

[channels, ~] = LoadTif(tifpath, ch_total, [ch_sample], 1);
sample_image = channels{ch_sample,1};
clear channels;

%Get stats on raw image
rawimg_stats = struct('imin', 0);
rawimg_stats.imin = min(sample_image, [], 'all');
rawimg_stats.imax = max(sample_image, [], 'all');
rawimg_stats.irange = rawimg_stats.imax - rawimg_stats.imin;
rawimg_stats.ipercentiles = NaN(pcheck_count,2);
rawimg_stats.ipercentiles(:,1) = percentile_values;
rawimg_stats.ipercentiles(:,2) = prctile(sample_image,percentile_values,'all');

top1 = find(sample_image >= rawimg_stats.ipercentiles(perc_index_99,2));
rawimg_stats.ipercentiles_top = NaN(pcheck_count,2);
rawimg_stats.ipercentiles_top(:,1) = percentile_values;
rawimg_stats.ipercentiles_top(:,2) = prctile(sample_image(top1),percentile_values,'all');
    
%Rescale raw
if rawimg_stats.irange < RAW_RESCALE_TRIGGER
    fprintf("WARNING: Raw image has low dynamic range. Rescale triggered!\n");
    sample_image = ((sample_image - rawimg_stats.imin) .* RAW_RESCALE_TARGET) ./ rawimg_stats.irange;
end

%Get stats on rescaled raw image
rsraw_stats = struct('imin', 0);
rsraw_stats.imin = min(sample_image, [], 'all');
rsraw_stats.imax = max(sample_image, [], 'all');
rsraw_stats.irange = rsraw_stats.imax - rsraw_stats.imin;
rsraw_stats.ipercentiles = NaN(pcheck_count,2);
rsraw_stats.ipercentiles(:,1) = percentile_values;
rsraw_stats.ipercentiles(:,2) = prctile(sample_image,percentile_values,'all');

top1 = find(sample_image >= rsraw_stats.ipercentiles(perc_index_99,2));
rsraw_stats.ipercentiles_top = NaN(pcheck_count,2);
rsraw_stats.ipercentiles_top(:,1) = percentile_values;
rsraw_stats.ipercentiles_top(:,2) = prctile(sample_image(top1),percentile_values,'all');

%Display raw image

%Apply filter
dead_pix_path = 'deadpix.mat';
RNA_Threshold_Common.saveDeadPixels(sample_image, dead_pix_path, true);

IMG3D = uint16(sample_image);
clear sample_image;
IMG3D = RNA_Threshold_Common.cleanDeadPixels(IMG3D, dead_pix_path, true);
delete(dead_pix_path);

IMG_filtered = RNA_Threshold_Common.applyGaussianFilter(IMG3D, gaussian_rad, 2);
IMG_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(IMG_filtered);
IMG_filtered = RNA_Threshold_Common.blackoutBorders(IMG_filtered, gaussian_rad+1, 0);
IMG_filtered = uint16(IMG_filtered);

%Get stats on filtered image
filtimg_stats = struct('imin', 0);
filtimg_stats.imin = min(IMG_filtered, [], 'all');
filtimg_stats.imax = max(IMG_filtered, [], 'all');
filtimg_stats.irange = filtimg_stats.imax - filtimg_stats.imin;
filtimg_stats.ipercentiles = NaN(pcheck_count,2);
filtimg_stats.ipercentiles(:,1) = percentile_values;
filtimg_stats.ipercentiles(:,2) = prctile(IMG_filtered,percentile_values,'all');

top1 = find(IMG_filtered >= filtimg_stats.ipercentiles(perc_index_99,2));
filtimg_stats.ipercentiles_top = NaN(pcheck_count,2);
filtimg_stats.ipercentiles_top(:,1) = percentile_values;
filtimg_stats.ipercentiles_top(:,2) = prctile(IMG_filtered(top1),percentile_values,'all');

%Display filtered image
max_proj = double(max(IMG3D,[],3));
max_proj_f = double(max(IMG_filtered,[],3));
Lmin = min(max_proj(:));
Lmax = median(max_proj(:)) + round(10 * std(max_proj(:)));
LminF = median(max_proj_f(:)) - round(0 * std(max_proj_f(:)));
LmaxF = median(max_proj_f(:)) + round(10 * std(max_proj_f(:)));

figure(1);
imshow(max_proj, [Lmin Lmax]);

figure(2);
imshow(max_proj_f, [LminF LmaxF]);

% ========================== Helper funcs ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end