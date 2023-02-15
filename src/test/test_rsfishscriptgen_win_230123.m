%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc';

JarPath = [DataDir '\tools\rsfish_scriptgen.jar'];

% ========================== Constants ==========================

ExtractedImgDir = 'C:\Users\hospelb.VUDS\Desktop\splitimg';
DataOutputDir = [DataDir '\imgproc\data\rsfish'];

THREADS = 2;

DO_IMG_SPLIT = false;
KEEP_CSVS = false;
USE_WSL = false;

% ========================== Load csv Table ==========================
InputTablePath = [DataDir filesep 'imgproc\test_images.csv'];
image_table = testutil_opentable(InputTablePath);

%ImageName='scrna_E2R2I5_CTT1';
GroupPrefix = 'sctc_E2R2_';
GroupSuffix = 'STL1';
GroupDirSfx = 'yeast_tc/E2R2/CH1';

NewTifDir = ['/img/' GroupDirSfx];
OutputDir = [DataOutputDir filesep replace(GroupDirSfx, '/', filesep)];

% ========================== Do things ==========================
addpath('./core');

imin = 65535;
imax = 0;
rec_count = size(image_table,1);
incl_count = 0;
voxel_size = [];
point_size = [];

newtifdir = [ExtractedImgDir replace(NewTifDir, '/', filesep)];
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if ~startsWith(iname, GroupPrefix); continue; end
    if ~isempty(GroupSuffix)
        if ~endsWith(iname, GroupSuffix); continue; end
    end
    
    %General steps:
    %   1. Isolate TIF channels if needed. Rearrange input paths.
    %   2. Generate wrapper script. This calls the Java module that writes
    %       the ijm script for this batch, then ImageJ itself.
    %       Then it runs a MATLAB script that bundles the output and
    %       deletes the temp csvs.
    %   3. Generate an sbatch command for submitting the main script
    
    %Determine if needs splitting, and if so what the new tif path is.
    ogtifpath = getTableValue(image_table, r, 'IMAGEPATH');
    chcount = getTableValue(image_table, r, 'CH_TOTAL');
    
    ogtifpath = [ImgDir replace(ogtifpath, '/', filesep)];
    newtifpath = [newtifdir filesep iname '.tif'];
    fprintf('Input TIF added: %s\n', newtifpath);
    if chcount > 1
        trgch = getTableValue(image_table, r, 'CHANNEL');
        
        %Open (need to check imin and imax)
        [channels, ~] = LoadTif(ogtifpath, chcount, [trgch], 1);
        my_image = channels{1,trgch};
        chmin = min(my_image, [], 'all');
        chmax = max(my_image, [], 'all');
        if chmin < imin; imin = chmin; end
        if chmax > imax; imax = chmax; end
        
        if DO_IMG_SPLIT
            tifops.overwrite = true;
            saveastiff(my_image, newtifpath, tifops);
        end
        clear my_image;
    else
        %Copy to new dir and rename if DO_IMG_SPLIT
        [channels, ~] = LoadTif(ogtifpath, 1, [1], 1);
        my_image = channels{1,1};
        chmin = min(my_image, [], 'all');
        chmax = max(my_image, [], 'all');
        if chmin < imin; imin = chmin; end
        if chmax > imax; imax = chmax; end
        
        if DO_IMG_SPLIT
            tifops.overwrite = true;
            saveastiff(my_image, newtifpath, tifops);
        end
        clear my_image;
    end
    incl_count = incl_count + 1;
    
    if isempty(voxel_size)
        voxel_size.x = getTableValue(image_table, r, 'VOXEL_X');
        voxel_size.y = getTableValue(image_table, r, 'VOXEL_Y');
        voxel_size.z = getTableValue(image_table, r, 'VOXEL_Z');
    end
    
    if isempty(point_size)
        point_size.x = getTableValue(image_table, r, 'POINT_X');
        point_size.y = getTableValue(image_table, r, 'POINT_Y');
        point_size.z = getTableValue(image_table, r, 'POINT_Z');
    end
    
end

%Generate command
if USE_WSL
    JarPath_Unix = unixify_winpath(JarPath);
    InputDir_Unix = unixify_winpath(newtifdir);
    OutputDir_Unix = unixify_winpath(OutputDir);

    fprintf('java -jar "%s"', JarPath_Unix);
    fprintf(' --input "%s"', InputDir_Unix);
    fprintf(' --output "%s"', OutputDir_Unix);
else
    fprintf('java -jar "%s"', JarPath);
    fprintf(' --input "%s"', newtifdir);
    fprintf(' --output "%s"', OutputDir);
end

%fprintf(' --minTh %d', TH_MIN);
%fprintf(' --maxTh %d', TH_MAX);
fprintf(' --voxsz "(%d,%d,%d)"', voxel_size.x, voxel_size.y, voxel_size.z);
fprintf(' --spotsz "(%d,%d,%d)"', point_size.x, point_size.y, point_size.z);
fprintf(' --threads %d', THREADS);
fprintf(' --imMin %d', imin);
fprintf(' --imMax %d', imax);
%fprintf(' --thdog %.3f', DOG_TH);

fprintf('\n');

% ========================== Helper funcs ==========================

function unix_path = unixify_winpath(winpath)
    unix_path = replace(winpath, 'C:', '/mnt/c');
    unix_path = replace(unix_path, 'D:', '/mnt/d');
    unix_path = replace(unix_path, '\', '/');
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end