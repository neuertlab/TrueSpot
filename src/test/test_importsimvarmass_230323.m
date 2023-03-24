%
%% BASE DIR

ImgProcBaseDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Constants ==========================

addpath('./core');

UseDir = [ImgProcBaseDir filesep 'img' filesep 'simvarmass'];
table_path = [UseDir filesep 'siminfo.csv'];

TifDir = [UseDir filesep 'tif'];
KeyDir = [UseDir filesep 'keybackup'];
TrashDir = [UseDir filesep 'trash'];

VOX_XY = 65;
VOX_Z = 200;
POINT_XY = 95;
POINT_Z = 210;
BG_LEVEL = 500;
CLUSTER_COUNT = 5;
CLUSTER_SIZE = 3;
RANDOM_SIG = 0.15;

% ========================== Read Image Table ==========================

ImageTable = readtable(table_path,'Delimiter',',','ReadVariableNames',true,'Format',...
    '%s%d%d%f%f');

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
    simparam.bg_level = BG_LEVEL;
    simparam.sigma_var = RANDOM_SIG;
    simparam.cluster_count = CLUSTER_COUNT;
    simparam.cluster_size = CLUSTER_SIZE;

    simparam.spots_seed = getTableValue(ImageTable, r, 'SpotsSeed');
    simparam.spots_actual = spot_count;
    simparam.amplitude_mean = getTableValue(ImageTable, r, 'Ampl');
    simparam.amplitude_var = getTableValue(ImageTable, r, 'AmplVar');
    simparam.bg_var = getTableValue(ImageTable, r, 'BGVar');

    %Save
    imgdat = imgdat16;
    key = key_out;
    save([UseDir filesep 'simvarmass_' simid '.mat'], 'imgdat', 'key', 'simparam');
    
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
