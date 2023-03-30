%
%% BASE DIR

%ImgProcBaseDir = 'D:\usr\bghos\labdat\imgproc';
ImgProcBaseDir = 'C:\Users\hospelb\labdata\imgproc';

% ========================== Constants ==========================

addpath('./core');

UseDir = [ImgProcBaseDir filesep 'img' filesep 'simytc'];
table_path = [UseDir filesep 'siminfo.csv'];

TifDir = [UseDir filesep 'tif'];
KeyDir = [UseDir filesep 'keybackup'];
TrashDir = [UseDir filesep 'trash'];

VOX_XY = 65;
VOX_Z = 200;
POINT_XY = 95;
POINT_Z = 210;
CLUSTER_COUNT = 3;
CLUSTER_SIZE = 6;
RANDOM_SIG = 0.15;

% ========================== Read Image Table ==========================

ImageTable = readtable(table_path,'Delimiter',',','ReadVariableNames',true,'Format',...
    '%s%d%d%f%d%f');

% ========================== Import Images ==========================

if ~isfolder(TifDir)
    mkdir(TifDir);
end
if ~isfolder(KeyDir)
    mkdir(KeyDir);
end
if ~isfolder(TrashDir)
    mkdir(TrashDir);
end

imgcount = size(ImageTable,1);

for r = 1:imgcount
    clear key;
    clear imgdat;
    clear simparam;

    simid = getTableValue(ImageTable, r, 1);
    fprintf('Processing %s (%d of %d)...\n', simid, r, imgcount);

    z = 1;
    simstem = [UseDir filesep simid];
    filetblname = [simstem '_imgz_' sprintf('%04d', z) '.csv'];
    while isfile(filetblname)
        z = z + 1;
        filetblname = [simstem '_imgz_' sprintf('%04d', z) '.csv'];
    end
    Z = z;

    %Open first file to determine X and Y
    z = 1;
    filetblname = [simstem '_imgz_' sprintf('%04d', (z-1)) '.csv'];
    raw_slice = double(csvread(filetblname));
    Y = size(raw_slice, 1);
    X = size(raw_slice, 2);

    %Read in
    img_dat = NaN(Y,X,Z);
    img_dat(:,:,z) = raw_slice(:,:);
    movefile(filetblname, [TrashDir filesep simid '_imgz_' sprintf('%04d', (z-1)) '.csv']);
    for z = 2:Z
        filetblname = [simstem '_imgz_' sprintf('%04d', (z-1)) '.csv'];
        raw_slice = double(csvread(filetblname));
        img_dat(:,:,z) = raw_slice(:,:);
        movefile(filetblname, [TrashDir filesep simid '_imgz_' sprintf('%04d', (z-1)) '.csv']);
    end
    imgdat_raw = uint16(img_dat); %For saving prefiltered image
    
    %Apply filter to sliiiiightly darken the edges of the image in xy
    mindark = 0.01;
    maxdark = 0.10;
    rdark = rand();
    dark_factor = mindark + (rdark * (maxdark-mindark));
    dark_noise_var = 0.05;
    
    ctr_x = (X+1)/2;
    ctr_y = (Y+1)/2;
    [xdist, ydist] = meshgrid(1:X,1:Y);
    xdist = abs(xdist - ctr_x);
    ydist = abs(ydist - ctr_y);
    xydist = sqrt((double(xdist).^2) + (double(ydist).^2));
    maxdist = sqrt((double(ctr_x).^2) + (double(ctr_y).^2));
    xydist = xydist ./ maxdist;
    xydark = xydist .* dark_factor;
    xydark = 1.0 - xydark;
    
    for z = 1:Z
        somenoise = 1.0 - (rand(Y,X) .* dark_noise_var);
        darkmtx = xydark .* somenoise;
        img_dat(:,:,z) = immultiply(img_dat(:,:,z), darkmtx(:,:));
    end
    
    %Apply filter to slightly blur edge z slices
    ctr_z = ((Z+1) / 2.1); %Juuuust below center
    zdist = [1:Z];
    zdist = abs(zdist - ctr_z);
    gauss_rad = 5;
    %Amount will range from 0 to between 1.1 and 2.5
    blur_scale = (rand() * 1.4) + 1.1;
    %For slices that get blur, expand in xy so that borders can be blurred
    %too
    %Slices with a value < 1.1 get no blur applied.
    zdistmax = max(zdist, [], 'all');
    zblur = (zdist ./ zdistmax) .* blur_scale;
    for z = 1:Z
        zblur_amt = zblur(z);
        %if zblur_amt < 1.1
        %    continue;
        %end
        
        gadd = gauss_rad + 1;
        GX = X + (gadd * 2);
        GY = Y + (gadd * 2);
        gslice = zeros(GY,GX);
        
        gx_min = gadd + 1;
        gx_max = GX - gadd;
        gy_min = gadd + 1;
        gy_max = GY - gadd;
        gslice(gy_min:gy_max, gx_min:gx_max) = img_dat(:,:,z);

        %Fill in borders
        imgedge = img_dat(:,:,z);
        imgedge(3:(Y-3),3:(X-3)) = NaN;
        slice_mean = nanmean(img_dat(:,:,z), 'all');
        edge_mean = nanmean(imgedge(:,:), 'all');
        clear imgedge;
        randamp = 50;
        fillnoise = rand(GY,GX) .* randamp;
        fillnoise = (fillnoise .* 2) - randamp;
        fillmtx = fillnoise + edge_mean;
        fillmtx(gy_min:gy_max, gx_min:gx_max) = 0;
        gslice = gslice + fillmtx;
        
        %Apply blur
        gslice = RNA_Threshold_Common.applyGaussianFilter(gslice, gauss_rad, zblur_amt);
        
        img_dat(:,:,z) = gslice(gy_min:gy_max, gx_min:gx_max);
    end
    

    %Convert to uint16
    imgdat16 = uint16(img_dat);
    clear img_dat;

    %Read key and copy raw file to backup
    filetblname = [simstem '_key.csv'];
    key_data = double(csvread(filetblname));
    spot_count = size(key_data,1);
    clear key_out;
    key_out(spot_count) = struct('x', 0, 'y', 0, 'z', 0, 'sigma_xy', 0.0, 'sigma_z', 0.0, 'amplitude', 0.0);
    for s = 1:spot_count
        key_out(s).x = uint16(key_data(s,3) + 1);
        key_out(s).y = uint16(key_data(s,2) + 1);
        key_out(s).z = uint16(key_data(s,1) + 1);
        key_out(s).sigma_xy = key_data(s,5);
        key_out(s).sigma_z = key_data(s,4);
        key_out(s).amplitude = key_data(s,6);
    end
    movefile(filetblname, [KeyDir filesep simid '_key.csv']);

    %Make a struct for the other sim parameters
    simparam = struct('voxel_x', VOX_XY);
    simparam.voxel_y = VOX_XY;
    simparam.voxel_z = VOX_Z;
    simparam.point_x = POINT_XY;
    simparam.point_y = POINT_XY;
    simparam.point_z = POINT_Z;
    simparam.sigma_var = RANDOM_SIG;
    simparam.cluster_count = CLUSTER_COUNT;
    simparam.cluster_size = CLUSTER_SIZE;
    simparam.xy_dark_factor = dark_factor;
    simparam.xy_dark_noise = dark_noise_var;
    simparam.zblur_max = blur_scale;

    simparam.spots_seed = getTableValue(ImageTable, r, 'SpotsSeed');
    simparam.spots_actual = spot_count;
    simparam.amplitude_mean = getTableValue(ImageTable, r, 'Ampl');
    simparam.amplitude_var = getTableValue(ImageTable, r, 'AmplVar');
    simparam.bg_level = getTableValue(ImageTable, r, 'BG');
    simparam.bg_var = getTableValue(ImageTable, r, 'BGVar');

    %Save
    imgdat = imgdat16;
    key = key_out;
    save([UseDir filesep 'simvarmass_' simid '.mat'], 'imgdat', 'key', 'simparam', 'imgdat_raw');
    
    tifpath = [TifDir filesep 'simvarmass_' simid '.tif'];
    if isfile(tifpath)
        delete(tifpath);
    end
    saveastiff(imgdat, tifpath);

end

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end