%
%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

img_paths = cell(48,1);
i = 1;

%----- Test
%img_paths{i,1} = [ImgDir '\img\sim\mESC_RNA_TMRLike_100x_3']; i = i+1;
img_paths{i,1} = [ImgDir '\img\sim\mESC_RNA_LEHB_100x_3']; i = i+1;

% ========================== Cycle ==========================

addpath('./core');
path_count = i-1;

for j = 1:path_count
    mypath = img_paths{j,1};
    %Find highest number slice file
    z = 0;
    filetblname = [mypath '_imgz_' sprintf('%04d', z) '.csv'];
    while isfile(filetblname)
        z = z + 1;
        filetblname = [mypath '_imgz_' sprintf('%04d', z) '.csv'];
    end
    Z = z;

    %Open first file to determine X and Y
    z = 1;
    filetblname = [mypath '_imgz_' sprintf('%04d', (z-1)) '.csv'];
    raw_coord_table = double(csvread(filetblname));
    Y = size(raw_coord_table, 1);
    X = size(raw_coord_table, 2);
   
    %Read in
    img_dat = NaN(Y,X,Z);
    img_dat(:,:,z) = raw_coord_table(:,:);
    for z = 2:Z
        filetblname = [mypath '_imgz_' sprintf('%04d', (z-1)) '.csv'];
        raw_coord_table = double(csvread(filetblname));
        img_dat(:,:,z) = raw_coord_table(:,:);
    end

    %Convert to uint16
    imgdat16 = uint16(img_dat);
    clear img_dat;

    %Load truthset
    filetblname = [mypath '_key.csv'];
    key_data = double(csvread(filetblname));
    spot_count = size(key_data,1);
    key_out(spot_count) = struct('x', 0, 'y', 0, 'z', 0, 'sigma_xy', 0.0, 'sigma_z', 0.0, 'amplitude', 0.0);
    for s = 1:spot_count
        key_out(s).x = uint16(key_data(s,3) + 1);
        key_out(s).y = uint16(key_data(s,2) + 1);
        key_out(s).z = uint16(key_data(s,1) + 1);
        key_out(s).sigma_xy = key_data(s,5);
        key_out(s).sigma_z = key_data(s,4);
        key_out(s).amplitude = key_data(s,6);
    end

    %Save
    imgdat = imgdat16;
    key = key_out;
    save([mypath '.mat'], 'imgdat', 'key');

    %(Optional) Apply blur and save again
%     xy_rad = 7;
%     z_rad = 2;
%     g_amt = 4;
% 
%     min_o = min(imgdat, [], 'all');
%     max_o = max(imgdat, [], 'all');
%     ao = double(min_o);
%     do = double(max_o - min_o);
% 
%     imgdat = RNA_Threshold_Common.applyGaussianFilter(imgdat, xy_rad, g_amt, z_rad, false);
%     imgdat = RNAUtils.medianifyBorder(imgdat, [xy_rad xy_rad z_rad]);
% 
%     min_f = min(imgdat, [], 'all');
%     max_f = max(imgdat, [], 'all');
%     af = double(min_f);
%     df = double(max_f - min_f);
%     props = (double(imgdat) - af) ./ df;
%     imgdat = round((props .* do) + ao);
% 
%     save([mypath '_blur.mat'], 'imgdat', 'key');

end