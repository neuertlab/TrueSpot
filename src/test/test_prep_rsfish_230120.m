%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

ScriptDir = 'C:\Users\hospelb.VUDS\Desktop\slurm';
%ScriptDir = 'C:\Users\bghos\Desktop\slurm';

ExtractedImgDir = 'C:\Users\hospelb.VUDS\Desktop\splitimg';

ClusterImageJPath = '~/imgproc/Fiji.app/ImageJ-linux64';
ClusterWorkDir = '/scratch/hospelb/imgproc';
ClusterSlurmDir = '/scratch/hospelb/imgproc/slurm/script';
ClusterScriptsDir = '/scratch/hospelb/scripts';

MatlabImportFunc = 'Main_RSFish2Mat';
JavaJarName = 'rsfish_scriptgen.jar';

MODULE_NAME = 'MATLAB/2018b';
MATLAB_DIR = [ClusterWorkDir '/matlab'];

% ========================== Constants ==========================

DETECT_THREADS = 8;
RAM_PER_CORE = 2;
HRS_PER_IMAGE = 1;

TH_MIN = 10;
TH_MAX = 1000;
DOG_TH = 0.007;

DO_IMG_SPLIT = true;
KEEP_CSVS = false;

% ========================== Load csv Table ==========================
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

NewTifDir = '/img/yeast_tc/E2R2/CH1';

%ImageName='scrna_E2R2I5_CTT1';
GroupPrefix = 'sctc_E2R2_';
GroupSuffix = 'STL1';

% ========================== Do things ==========================
addpath('./core');

imin = 65535;
imax = 0;
rec_count = size(image_table,1);
incl_count = 0;
voxel_size = [];
point_size = [];
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
    fieldstr = getTableValue(image_table, r, 'CH_TOTAL');
    chcount = str2num(fieldstr);
    
    ogtifpath = [ImgDir replace(ogtifpath, '/', filesep)];
    newtifpath = [ExtractedImgDir replace(NewTifDir, '/', filesep) filesep iname '.tif'];
    fprintf('Input TIF added: %s\n', newtifpath);
    if chcount > 1
        fieldstr = getTableValue(image_table, r, 'CHANNEL');
        trgch = str2num(fieldstr);
        
        %Open (need to check imin and imax)
        [channels, ~] = LoadTif(ogtifpath, chcount, [trgch], 1);
        my_image = channels{1,trgch};
        chmin = min(my_image, [], 'all');
        chmax = max(my_image, [], 'all');
        if chmin < imin; imin = chmin; end
        if chmax < imax; imax = chmax; end
        
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
        if chmax < imax; imax = chmax; end
        
        if DO_IMG_SPLIT
            tifops.overwrite = true;
            saveastiff(my_image, newtifpath, tifops);
        end
        clear my_image;
    end
    incl_count = incl_count + 1;
    
    if isempty(voxel_size)
        voxel_size.x = str2num(getTableValue(image_table, r, 'VOXEL_X'));
        voxel_size.y = str2num(getTableValue(image_table, r, 'VOXEL_Y'));
        voxel_size.z = str2num(getTableValue(image_table, r, 'VOXEL_Z'));
    end
    
    if isempty(point_size)
        point_size.x = str2num(getTableValue(image_table, r, 'POINT_X'));
        point_size.y = str2num(getTableValue(image_table, r, 'POINT_Y'));
        point_size.z = str2num(getTableValue(image_table, r, 'POINT_Z'));
    end
    
end

%Write script
script_path = [ScriptDir filesep 'rsfish_batch.sh'];
script_file = fopen(script_path, 'w');

fprintf(script_file, '#!/bin/bash\n\n');

%Java script generator
fprintf(script_file, 'JAVA_SCRGEN_PATH=%s\n', [ClusterScriptsDir '/' JavaJarName]);
fprintf(script_file, 'module load Java/13.0.2\n');
fprintf(script_file, 'java -jar ${JAVA_SCRGEN_PATH}');
fprintf(script_file, ' --input "%s"', [ClusterWorkDir NewTifDir]);

ds = datestr(datetime, 'yyyymmddHHMM');
run_out_dir = [ClusterWorkDir '/data/rsfish/' ds];
fprintf(script_file, ' --output "%s"', run_out_dir);

fprintf(script_file, ' --minTh %d', TH_MIN);
fprintf(script_file, ' --maxTh %d', TH_MAX);

fprintf(script_file, ' --voxsz "(%d,%d,%d)"', voxel_size.x, voxel_size.y, voxel_size.z);
fprintf(script_file, ' --spotsz "(%d,%d,%d)"', point_size.x, point_size.y, point_size.z);
fprintf(script_file, ' --threads %d', DETECT_THREADS);
fprintf(script_file, ' --imMin %d', imin);
fprintf(script_file, ' --imMax %d', imax);
fprintf(script_file, ' --thdog %f', DOG_TH);
fprintf(script_file, '\n\n');

fprintf(script_file, 'IMAGEJ_PATH=%s\n',ClusterImageJPath);
fprintf(script_file, 'IJM_PATH=%s\n',[run_out_dir '/rsfish.ijm']);
fprintf(script_file, 'chmod 770 ${IJM_PATH}\n');
fprintf(script_file, '${IMAGEJ_PATH} --headless --run ${IJM_PATH} %s\n\n', [run_out_dir '/rsfish.log']);

fprintf(script_file, 'module load %s\n', MODULE_NAME);
fprintf(script_file, 'cd %s\n', MATLAB_DIR);
fprintf(script_file, 'matlab -nodisplay -nosplash -logfile "%s_mat.log" -r "cd %s; ', [run_out_dir 'rsfish'], MATLAB_DIR);
fprintf(script_file, '%s(', MatlabImportFunc);
fprintf(script_file, '''%s''', run_out_dir);

if KEEP_CSVS
    fprintf(script_file, ', ''-keepcsvs''');
end
fprintf(script_file, '\n');

fclose(script_file);

%Slurm command
cluster_script_path = [ClusterSlurmDir '/rsfish_batch.sh'];
fprintf('chmod 770 %s\n', cluster_script_path);
fprintf('mkdir -p %s\n', run_out_dir);
fprintf('sbatch');
fprintf(' --job-name="%s"', ['RSFish_' ds]);
fprintf(' --cpus-per-task=%d', DETECT_THREADS);

hours_req = ceil(incl_count * HRS_PER_IMAGE);
fprintf(' --time=%d:00:00', hours_req);
fprintf(' --mem=%dg', (DETECT_THREADS * RAM_PER_CORE));
fprintf(' --error="%s"', [run_out_dir '/rsfish_sbatch.err']);
fprintf(' --output="%s"', [run_out_dir '/rsfish_sbatch.out']);
fprintf(' %s\n', cluster_script_path);

% ========================== Helper funcs ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end