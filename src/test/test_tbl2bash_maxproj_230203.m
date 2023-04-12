%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

ScriptDir = 'C:\Users\hospelb.VUDS\Desktop\slurm';
%ScriptDir = 'C:\Users\bghos\Desktop\slurm';

ClusterWorkDir = '/nobackup/p_neuert_lab/hospelb/imgproc';
ClusterSlurmDir = '/nobackup/p_neuert_lab/hospelb/imgproc/slurm/script';
ClusterScriptsDir = '/nobackup/p_neuert_lab/hospelb/scripts';

% ========================== Constants ==========================

DETECT_THREADS = 2;
RAM_PER_CORE = 8;
HB_PARALLEL_HR = 2;
HB_SERIAL_HR = 8;

QUANT_SERIAL_HR = 6;
QUANT_PARALLEL_HR = 4;
RAM_PER_CORE_QUANT = 12;
DETECT_THREADS_QUANT = 1;
QUANT_DO_CLOUDS = false;
QUANT_FIXED_TH = 0;

TH_MIN = 10;
TH_MAX = 500;
Z_TRIM = 0;
NO_DPC = false;

RUN_HB = true;
RUN_QUANT = true;

MODULE_NAME = 'MATLAB/2018b';
MATLAB_DIR = [ClusterWorkDir '/matlab'];

% ========================== Load csv Table ==========================
addpath('./core');
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

GroupPrefix = 'sctc_E1R1_';
GroupSuffix = [];

% ========================== Generate Bash Script & Slurm Command ==========================

script_master = fopen([ScriptDir filesep 'runall_maxp.sh'], 'w');
fprintf(script_master, '#!/bin/bash\n\n');
fprintf(script_master, 'SCRIPTDIR=%s\n', ClusterSlurmDir);

rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if ~startsWith(iname, GroupPrefix); continue; end
    if ~isempty(GroupSuffix)
        if ~endsWith(iname, GroupSuffix); continue; end
    end
    
    hb_outstem = getTableValue(image_table, r, 'OUTSTEM');
    [hb_outdir_base,~,~] = fileparts(hb_outstem);
    mpoutdir = [hb_outdir_base '/maxproj'];
    mpoutstem = [mpoutdir '/' iname '_max_proj'];
    
    ipath = getTableValue(image_table, r, 'IMAGEPATH');
    if endsWith(ipath, '.mat')
        ipath = replace(ipath, '.mat', '.tif');
    end
    
    %HB Script
    if RUN_HB
        script_hb = fopen([ScriptDir filesep iname '_mphb.sh'], 'w');
        fprintf(script_hb, '#!/bin/bash\n\n');
        fprintf(script_hb, 'module load %s\n', MODULE_NAME);
        fprintf(script_hb, 'cd %s\n', MATLAB_DIR);
        fprintf(script_hb, 'echo $(date +''%%Y/%%m/%%d %%H:%%M:%%S:%%3N'')\n');
        fprintf(script_hb, 'matlab -nodisplay -nosplash -logfile "%s_mat.log" -r "cd %s; ', [ClusterWorkDir mpoutstem], MATLAB_DIR);
        fprintf(script_hb, 'Main_RNASpots(');

        printMatArg(script_hb, 'imgname', iname, false);
        if endsWith(ipath, '.mat')
            printMatArg(script_hb, 'matimg', [ClusterWorkDir ipath], true);
        elseif endsWith(ipath, '.tif')
            printMatArg(script_hb, 'tif', [ClusterWorkDir ipath], true);
        end

        printMatArg(script_hb, 'outdir', [ClusterWorkDir mpoutdir], true);
        csegdir = getTableValue(image_table, r, 'CELLSEG_DIR');
        csegsfx = getTableValue(image_table, r, 'CELLSEG_SFX');
        csegpath = [];
        if ~strcmp(csegdir, '.')
            csegpath = [ClusterWorkDir csegdir '/Lab_' csegsfx '.mat'];
            printMatArg(script_hb, 'cellseg', csegpath, true);
        end

        printMatArg(script_hb, 'chsamp', num2str(getTableValue(image_table, r, 'CHANNEL')), true);
        printMatArg(script_hb, 'chtrans', num2str(getTableValue(image_table, r, 'CH_LIGHT')), true);
        printMatArg(script_hb, 'chtotal', num2str(getTableValue(image_table, r, 'CH_TOTAL')), true);

        sfield = getTableValue(image_table, r, 'CONTROL_TIF');
        if ~strcmp(sfield, '.')
            printMatArg(script_hb, 'ctrltif', [ClusterWorkDir sfield], true);
            printMatArg(script_hb, 'chctrsamp', num2str(getTableValue(image_table, r, 'CHANNEL')), true);
            printMatArg(script_hb, 'chctrtotal', num2str(getTableValue(image_table, r, 'CH_TOTAL')), true);
        end

        printMatArg(script_hb, 'thmin', num2str(TH_MIN), true);
        printMatArg(script_hb, 'thmax', num2str(TH_MAX), true);
        printMatArg(script_hb, 'ztrim', num2str(Z_TRIM), true);
        printSpecificityArg(script_hb, getTableValue(image_table, r, 'THRESH_SETTING'));

        x_str = num2str(getTableValue(image_table, r, 'VOXEL_X'));
        y_str = num2str(getTableValue(image_table, r, 'VOXEL_Y'));
        z_str = num2str(getTableValue(image_table, r, 'VOXEL_Z'));
        printMatArg(script_hb, 'voxelsize', ['(' x_str ',' y_str ',' z_str ')'], true);
        x_str = num2str(getTableValue(image_table, r, 'POINT_X'));
        y_str = num2str(getTableValue(image_table, r, 'POINT_Y'));
        z_str = num2str(getTableValue(image_table, r, 'POINT_Z'));
        printMatArg(script_hb, 'expspotsize', ['(' x_str ',' y_str ',' z_str ')'], true);

        printMatArg(script_hb, 'probetype', getTableValue(image_table, r, 'PROBE'), true);
        printMatArg(script_hb, 'target', getTableValue(image_table, r, 'TARGET'), true);
        printMatArg(script_hb, 'targettype', getTableValue(image_table, r, 'TARGET_TYPE'), true);
        printMatArg(script_hb, 'species', getTableValue(image_table, r, 'SPECIES'), true);
        printMatArg(script_hb, 'celltype', getTableValue(image_table, r, 'CELLTYPE'), true);

        %Only 3D has explicit multithreading
%         if DETECT_THREADS > 1
%             printMatArg(script_hb, 'threads', num2str(DETECT_THREADS), true);
%         end

        if NO_DPC
            printMatFlagArg(script_hb, 'nodpc', true);
        end
        printMatFlagArg(script_hb, 'verbose', true);
        printMatFlagArg(script_hb, 'maxzproj', true);

        if RUN_QUANT & ~isempty(csegpath)
            fprintf(script_hb, '); Main_QuickCellQuant_FromRun(');
            fprintf(script_hb, '''%s'', ', [ClusterWorkDir mpoutdir '/' iname '_max_proj']);
            fprintf(script_hb, '''%s''', csegpath);
            %fprintf(script_hb, '''%s'', ', [mpoutstem '_quickQuant.mat']);
            if QUANT_FIXED_TH > 0
                useth_idx = QUANT_FIXED_TH - TH_MIN + 1;
                fprintf(script_hb, ', ''%d''', useth_idx);
            end
            %dim_x = getTableValue(image_table, r, 'IDIM_X');
            %dim_y = getTableValue(image_table, r, 'IDIM_Y');
            %dim_z = getTableValue(image_table, r, 'IDIM_Z');
            %fprintf(script_hb, '''(%d,%d,%d)''', dim_x, dim_y, dim_z);
        end
        fprintf(script_hb, '); quit;"\n');
        fprintf(script_hb, 'echo $(date +''%%Y/%%m/%%d %%H:%%M:%%S:%%3N'')\n');
        fprintf(script_hb, 'module unload\n');
        
        fclose(script_hb);
    end
    
    fprintf(script_master, 'if [ -s "%s" ]; then\n', [ClusterWorkDir ipath]);
    
    if RUN_HB
        fprintf(script_master, '\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_mphb.sh']);
        fprintf(script_master, '\tmkdir -p "%s"\n', [ClusterWorkDir mpoutdir]);
        fprintf(script_master, '\tsbatch');
        fprintf(script_master, ' --job-name="%s"', ['SpotDetectMP_' iname]);
        if DETECT_THREADS > 1
            fprintf(script_master, ' --cpus-per-task=%d', DETECT_THREADS);
            fprintf(script_master, ' --time=%d:00:00', HB_PARALLEL_HR);
            fprintf(script_master, ' --mem=%dg', (RAM_PER_CORE * DETECT_THREADS));
        else
            fprintf(script_master, ' --cpus-per-task=2');
            fprintf(script_master, ' --time=%d:00:00', HB_SERIAL_HR);
            fprintf(script_master, ' --mem=%dg', RAM_PER_CORE);
        end
        fprintf(script_master, ' --error="%s"', [ClusterWorkDir mpoutstem '_slurm.err']);
        fprintf(script_master, ' --output="%s"', [ClusterWorkDir mpoutstem '_slurm.out']);
        fprintf(script_master, ' "${SCRIPTDIR}/%s"\n', [iname '_mphb.sh']);
    end
    fprintf(script_master, 'fi\n\n');
end

fclose(script_master);

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

function printMatFlagArg(fhandle, key, leading_comma)
    if leading_comma
        fprintf(fhandle, ', ');
    end
    fprintf(fhandle, '''-%s''', key);
end


function printMatArg(fhandle, key, value, leading_comma)
    if leading_comma
        fprintf(fhandle, ', ');
    end
    fprintf(fhandle, '''-%s'', ''%s''', key, value);
end

function printSpecificityArg(fhandle, t_setting)
    if t_setting == 1
        printMatArg(fhandle, 'sensitivity', '2', true);
    elseif t_setting == 2
        printMatArg(fhandle, 'sensitivity', '1', true);
    elseif t_setting == 3
        printMatArg(fhandle, 'sensitivity', '0', true);
    elseif t_setting == 4
        printMatArg(fhandle, 'specificity', '1', true);
    elseif t_setting == 5
        printMatArg(fhandle, 'specificity', '2', true);
    end
end