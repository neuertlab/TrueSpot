%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

ScriptDir = 'C:\Users\hospelb.VUDS\Desktop\slurm';
%ScriptDir = 'C:\Users\bghos\Desktop\slurm';

ExtractedImgDir = 'C:\Users\hospelb.VUDS\Desktop\splitimg';

ClusterWorkDir = '/nobackup/p_neuert_lab/hospelb/imgproc';
ClusterSlurmDir = '/nobackup/p_neuert_lab/hospelb/imgproc/slurm/script';
ClusterScriptsDir = '/nobackup/p_neuert_lab/hospelb/scripts';
ClusterRsFishPath = '/home/hospelb/Git/RS-FISH/rs-fish';
ScriptAdjJarPath = [ClusterScriptsDir '/rsfish_scriptadj.jar'];

MatlabImportFunc = 'Main_RSFish2Mat';

MODULE_NAME = 'MATLAB/2018b';
JAVA_MODULE = 'Java/13.0.2';
MATLAB_DIR = [ClusterWorkDir '/matlab'];

% ========================== Constants ==========================
addpath('./core');

DETECT_THREADS = 8;
RAM_PER_CORE = 8;
HRS_PER_IMAGE = 12;

OVERWRITE = false;

USE_ANISOTROPY = false;

TH_ITR = 250;
TH_NTR = 0.1/TH_ITR;
MIN_TH = 0.0036;
%MIN_TH = TH_NTR;
START_I = 1; %Maybe cutting out the lowest thresholds will help the memory problem?
while ((START_I * TH_NTR) < MIN_TH)
    START_I = START_I + 1;
end

START_I = 3; %Override

DO_IMG_SPLIT = false;

% ========================== Load csv Table ==========================
%InputTablePath = [DataDir filesep 'test_images_simytc.csv'];
%InputTablePath = [DataDir filesep 'test_images_simvarmass.csv'];
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

RSDirTail = '/yeast_tc/E2R1/CH1';
NewTifDir = ['/img' RSDirTail];

%ImageName='scrna_E2R2I5_CTT1';
GroupPrefix = 'sctc_E2R1_';
GroupSuffix = 'STL1';

% ========================== Do things ==========================

master_script_path = [ScriptDir filesep 'rs_runall.sh'];
master_script = fopen(master_script_path, 'w');

fprintf(master_script, '#!/bin/bash\n\n');
fprintf(master_script, 'module load %s\n', JAVA_MODULE);
fprintf(master_script, 'SCRIPTDIR=%s\n', ClusterSlurmDir);
fprintf(master_script, 'ADJJARPATH=%s\n\n', ScriptAdjJarPath);

rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if ~startsWith(iname, GroupPrefix); continue; end
    if ~isempty(GroupSuffix)
        if ~endsWith(iname, GroupSuffix); continue; end
    end
    
    %Determine if needs splitting, and if so what the new tif path is.
    ogtifpath = getTableValue(image_table, r, 'IMAGEPATH');
    chcount = getTableValue(image_table, r, 'CH_TOTAL');
    
    ogtifpath = [ImgDir replace(ogtifpath, '/', filesep)];
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
    rsoutdir = [ClusterWorkDir '/data/rsfish' RSDirTail '/' iname];
    
    
    my_script_path = [ScriptDir filesep iname '_rsfish.sh'];
    my_script = fopen(my_script_path, 'w');
    %Have line in the script that makes it quit (success) if output file
    %from rs fish at that threshold is empty (no point in doing higher
    %thresholds)
    
    vox_x = double(getTableValue(image_table, r, 'VOXEL_X'));
    vox_z = double(getTableValue(image_table, r, 'VOXEL_Z'));
    pnt_x = double(getTableValue(image_table, r, 'POINT_X'));
    
    anisotropy = (vox_x/vox_z) * 2.0;
    sigma = (pnt_x/vox_x)/2.25;
    
    fprintf(my_script, '#!/bin/bash\n\n');
    fprintf(my_script, 'RSFISH_PATH=%s\n', ClusterRsFishPath);
    fprintf(my_script, 'OUTPUT_DIR=%s\n', rsoutdir);
    fprintf(my_script, 'if [ ! -d ${OUTPUT_DIR} ]; then\n');
    fprintf(my_script, '\tmkdir -p "${OUTPUT_DIR}"\n');
    fprintf(my_script, 'fi\n');
    fprintf(my_script, 'for i in {%d..%d}\n', START_I, TH_ITR);
    fprintf(my_script, 'do\n');
    fprintf(my_script, '\tOUTPUT_PATH=${OUTPUT_DIR}/thitr_${i}.csv\n');
    fprintf(my_script, '\tTH_DOG=$(echo "scale = 5; %f * $i" | bc)\n', TH_NTR);
    fprintf(my_script, '\techo "itr $i, TH_DOG = ${TH_DOG}"\n');
    fprintf(my_script, '\techo $(date +''%%Y/%%m/%%d %%H:%%M:%%S:%%3N'')\n');
    fprintf(my_script, '\t${RSFISH_PATH}');
    if USE_ANISOTROPY
        fprintf(my_script, ' -a %.2f', anisotropy);
    end
    fprintf(my_script, ' -i "%s"', cluster_tif_path);
    fprintf(my_script, ' -o "${OUTPUT_PATH}"');
    fprintf(my_script, ' -s %.2f', sigma);
    fprintf(my_script, ' -t ${TH_DOG}\n');
    fprintf(my_script, '\techo $(date +''%%Y/%%m/%%d %%H:%%M:%%S:%%3N'')\n');
    
    %If output empty, quit.
    fprintf(my_script, '\tif [ ! -s "${OUTPUT_PATH}" ]; then\n');
    fprintf(my_script, '\t\tbreak\n');
    fprintf(my_script, '\tfi\n');
    
    fprintf(my_script, 'done\n\n');
    
    %Convert to mat
    rsoutstem = [rsoutdir '/RSFISH_' iname];
    fprintf(my_script, 'module load %s\n', MODULE_NAME);
    fprintf(my_script, 'cd %s\n', MATLAB_DIR);
	fprintf(my_script, 'matlab -nodisplay -nosplash -logfile "%s_mat.log" -r "cd %s; ', [rsoutdir '/rsfish'], MATLAB_DIR);
	fprintf(my_script, '%s(', MatlabImportFunc);
    fprintf(my_script, '''%s''', rsoutdir);
    fprintf(my_script, ', ''%s''); quit;"\n', rsoutstem);
    fprintf(my_script, 'module unload\n');
    
    fclose(my_script);
    
    %Add lines to master script
    fprintf(master_script, 'if [ -s "%s" ]; then\n', cluster_tif_path);
    fprintf(master_script, '\tif [ ! -d "%s" ]; then\n', rsoutdir);
    fprintf(master_script, '\t\tmkdir -p "%s"\n', rsoutdir);
    fprintf(master_script, '\tfi\n');
    fprintf(master_script, '\tchmod 770 ${SCRIPTDIR}/%s\n', [iname '_rsfish.sh']);
    
    if ~OVERWRITE
        fprintf(master_script, '\tif [ ! -s "%s" ]; then\n', [rsoutstem '_coordTable.mat']);
        
        %Check if we need to adjust the start point...
        fprintf(master_script, '\t\tjava -jar "${ADJJARPATH}" "%s" "%s"\n', rsoutdir, ['${SCRIPTDIR}/' iname '_rsfish.sh']);
        fprintf(master_script, '\t\tsbatch');
    else
        fprintf(master_script, '\tsbatch');
    end
    
    fprintf(master_script, ' --job-name="%s"', ['RSFISH_' iname]);
    if DETECT_THREADS > 1
        fprintf(master_script, ' --cpus-per-task=%d', DETECT_THREADS);
        fprintf(master_script, ' --time=%d:00:00', HRS_PER_IMAGE);
        fprintf(master_script, ' --mem=%dg', (RAM_PER_CORE * DETECT_THREADS));
    else
        fprintf(master_script, ' --cpus-per-task=2');
        fprintf(master_script, ' --time=%d:00:00', HRS_PER_IMAGE);
        fprintf(master_script, ' --mem=%dg', RAM_PER_CORE);
    end
    fprintf(master_script, ' --error="%s"', [rsoutdir '/rsfish_slurm.err']);
    fprintf(master_script, ' --output="%s"', [rsoutdir '/rsfish_slurm.out']);
    fprintf(master_script, ' "${SCRIPTDIR}/%s"\n', [iname '_rsfish.sh']);
    
    if ~OVERWRITE
        fprintf(master_script, '\telse\n');
        fprintf(master_script, '\t\techo -e "RS-FISH run for %s found! Not resubmitting..."\n', iname);
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
