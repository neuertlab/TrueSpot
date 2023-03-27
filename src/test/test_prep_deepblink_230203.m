%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

ScriptDir = 'C:\Users\hospelb.VUDS\Desktop\slurm';
%ScriptDir = 'C:\Users\bghos\Desktop\slurm';

ExtractedImgDir = 'C:\Users\hospelb.VUDS\Desktop\splitimg';

ClusterWorkDir = '/scratch/hospelb/imgproc';
ClusterSlurmDir = '/scratch/hospelb/imgproc/slurm/script';
ClusterScriptsDir = '/scratch/hospelb/scripts';
ClusterPyenvDir = '/scratch/hospelb/pyvenv';
ClusterModelsDir = '/scratch/hospelb/imgproc/data/deepblink_training/og_models';

PyenvModule = 'deepblink';
MatlabImportFunc = 'Main_DeepBlink2Mat';

MODULE_NAME = 'MATLAB/2018b';
MATLAB_DIR = [ClusterWorkDir '/matlab'];

% ========================== Constants ==========================

DETECT_THREADS = 4;
RAM_PER_CORE = 16;
HRS_PER_IMAGE = 4;

MIN_PROB = 0.01;

DO_IMG_SPLIT = false;
OVERWRITE = false;

OutDirTail = '/mESC4d';
ModelName = 'smfish';

NewTifDir = ['/img' OutDirTail];

% ========================== Load csv Table ==========================
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

%ImageName='scrna_E2R2I5_CTT1';
GroupPrefix = 'mESC4d_';
GroupSuffix = [];

% ========================== Do things ==========================
addpath('./core');

master_script_path = [ScriptDir filesep 'deep_runall.sh'];
master_script = fopen(master_script_path, 'w');

fprintf(master_script, '#!/bin/bash\n\n');
fprintf(master_script, 'SCRIPTDIR=%s\n', ClusterSlurmDir);

rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if ~startsWith(iname, GroupPrefix); continue; end
    if ~isempty(GroupSuffix)
        if ~endsWith(iname, GroupSuffix); continue; end
    end
    
    ogtifpath = getTableValue(image_table, r, 'IMAGEPATH');
    chcount = getTableValue(image_table, r, 'CH_TOTAL');
    
    ogtifpath = [ImgDir replace(ogtifpath, '/', filesep)];
    [~,ogtifname,~] = fileparts(ogtifpath);
    newtifpath = [ExtractedImgDir replace(NewTifDir, '/', filesep) filesep iname '.tif'];
    
    if chcount > 1
        if DO_IMG_SPLIT
            fprintf('Input TIF added: %s\n', newtifpath);
            trgch = getTableValue(image_table, r, 'CHANNEL');
            [channels, ~] = LoadTif(ogtifpath, chcount, [trgch], 1);
            my_image = channels{trgch,1};
            my_image = uint16(my_image);
            tifops.overwrite = true;
            tifops.color = false;
            saveastiff(my_image, newtifpath, tifops);
            clear my_image;
        end
        cluster_tif_path = [ClusterWorkDir NewTifDir '/' iname '.tif'];
    else
        short_tif_path = getTableValue(image_table, r, 'IMAGEPATH');
        if endsWith(short_tif_path, '.mat')
            [tif_pref, filename, ~] = fileparts(short_tif_path);
            cluster_tif_path = [ClusterWorkDir tif_pref '/' filename '.tif'];
        else
            cluster_tif_path = [ClusterWorkDir short_tif_path];
        end
    end
    dboutdir = [ClusterWorkDir '/data/deepblink' OutDirTail '/' iname];
    
    my_script_path = [ScriptDir filesep iname '_deepblink.sh'];
    my_script = fopen(my_script_path, 'w');
    
    fprintf(my_script, '#!/bin/bash\n\n');
    fprintf(my_script, 'PYVENV_DIR=%s\n', ClusterPyenvDir);
    fprintf(my_script, 'OUTPUT_DIR=%s\n', dboutdir);
    fprintf(my_script, 'if [ ! -d ${OUTPUT_DIR} ]; then\n');
    fprintf(my_script, '\tmkdir -p "${OUTPUT_DIR}"\n');
    fprintf(my_script, 'fi\n\n');
    
    fprintf(my_script, 'module load GCCcore/.10.2.0\n');
    fprintf(my_script, 'module load Python/3.8.6\n');
    fprintf(my_script, 'source ${PYVENV_DIR}/%s/bin/activate\n', PyenvModule);
    fprintf(my_script, 'deepblink check "%s"\n', cluster_tif_path);
    fprintf(my_script, 'echo $(date +''%%Y/%%m/%%d %%H:%%M:%%S:%%3N'')\n');
    fprintf(my_script, 'deepblink predict -i "%s" -o "${OUTPUT_DIR}" -m "%s" -p %.2f -ps 1\n', ...
        cluster_tif_path, [ClusterModelsDir '/deepblink_' ModelName '.h5'], MIN_PROB);
    fprintf(my_script, 'echo $(date +''%%Y/%%m/%%d %%H:%%M:%%S:%%3N'')\n');
    fprintf(my_script, 'deactivate\n\n');
    
    dboutstem = [dboutdir '/DeepBlink_' iname];
    fprintf(my_script, 'module load %s\n', MODULE_NAME);
    fprintf(my_script, 'cd %s\n', MATLAB_DIR);
	fprintf(my_script, 'matlab -nodisplay -nosplash -logfile "%s_mat.log" -r "cd %s; ', [dboutdir '/dbimport'], MATLAB_DIR);
	fprintf(my_script, '%s(', MatlabImportFunc);
    fprintf(my_script, '''%s''', [dboutdir '/' ogtifname '.csv']);
    fprintf(my_script, ', ''%s''); quit;"\n', dboutstem);
    fprintf(my_script, 'module unload\n');
    
    fclose(my_script);
    
    %Add lines to master script
    fprintf(master_script, 'if [ -s "%s" ]; then\n', cluster_tif_path);
    fprintf(master_script, '\tif [ ! -d "%s" ]; then\n', dboutdir);
    fprintf(master_script, '\t\tmkdir -p "%s"\n', dboutdir);
    fprintf(master_script, '\tfi\n');
    fprintf(master_script, '\tchmod 750 ${SCRIPTDIR}/%s\n', [iname '_deepblink.sh']);
    
    if ~OVERWRITE
        fprintf(master_script, '\tif [ ! -s "%s" ]; then\n', [dboutstem '_coordTable.mat']);
        fprintf(master_script, '\t\tsbatch');
    else
        fprintf(master_script, '\tsbatch');
    end
    
    fprintf(master_script, ' --job-name="%s"', ['DeepBlink_' iname]);
    if DETECT_THREADS > 1
        fprintf(master_script, ' --cpus-per-task=%d', DETECT_THREADS);
        fprintf(master_script, ' --time=%d:00:00', HRS_PER_IMAGE);
        fprintf(master_script, ' --mem=%dg', (RAM_PER_CORE * DETECT_THREADS));
    else
        fprintf(master_script, ' --cpus-per-task=2');
        fprintf(master_script, ' --time=%d:00:00', HRS_PER_IMAGE);
        fprintf(master_script, ' --mem=%dg', RAM_PER_CORE);
    end
    fprintf(master_script, ' --error="%s"', [dboutdir '/deepblink_slurm.err']);
    fprintf(master_script, ' --output="%s"', [dboutdir '/deepblink_slurm.out']);
    fprintf(master_script, ' "${SCRIPTDIR}/%s"\n', [iname '_deepblink.sh']);
    
    if ~OVERWRITE
        fprintf(master_script, '\telse\n');
        fprintf(script_master, '\t\techo -e "DeepBlink run for %s found! Not resubmitting..."\n', iname);
        fprintf(master_script, '\tfi\n');
    end
    
    fprintf(master_script, 'fi\n');
end

fclose(master_script);

% ========================== Helper funcs ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end