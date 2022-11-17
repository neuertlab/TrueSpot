%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

ScriptDir = 'C:\Users\hospelb.VUDS\Desktop\slurm';

ClusterWorkDir = '/scratch/hospelb/imgproc';
ClusterSlurmDir = '/scratch/hospelb/imgproc/slurm/script';
ClusterScriptsDir = '/scratch/hospelb/scripts';
ClusterPyenvDir = '/scratch/hospelb/pyvenv';

% ========================== Constants ==========================

DETECT_THREADS = 4;
RAM_PER_CORE=16;

TH_MIN = 10;
TH_MAX = 500;
Z_TRIM = 3;
BF_SOBJSZ = 50;
BF_NUCSZ = 200;
BF_RESCALE = false;

RUN_HB = true;
RUN_BF = false;

MODULE_NAME = 'MATLAB/2018b';
MATLAB_DIR = [ClusterWorkDir '/matlab'];

% ========================== Load csv Table ==========================
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

%ImageName='scrna_E2R2I5_CTT1';
GroupPrefix='mESC4d_';

% ========================== Find Record ==========================
addpath('./core');

% rec_row = 0;
% rec_count = size(image_table,1);
% for r = 1:rec_count
%     myname = getTableValue(image_table, r, 'IMGNAME');
%     if strcmp(myname, ImageName); rec_row = r; break; end
% end
% 
% if rec_row < 1
%     fprintf("Couldn't find image: %s!\n", ImageName);
%     return;
% end

% ========================== Generate Bash Script & Slurm Command ==========================

script_master = fopen([ScriptDir filesep 'runall.sh'], 'w');
fprintf(script_master, '#!/bin/bash\n\n');
fprintf(script_master, 'SCRIPTDIR=%s\n', ClusterSlurmDir);

rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if ~startsWith(iname, GroupPrefix); continue; end
    hb_outstem = getTableValue(image_table, r, 'OUTSTEM');
    
    %HB Script
    if RUN_HB
        script_hb = fopen([ScriptDir filesep iname '_hb.sh'], 'w');
        fprintf(script_hb, '#!/bin/bash\n\n');
        fprintf(script_hb, 'module load %s\n', MODULE_NAME);
        fprintf(script_hb, 'cd %s\n', MATLAB_DIR);
        fprintf(script_hb, 'matlab -nodisplay -nosplash -logfile "%s_mat.log" -r "cd %s; ', [ClusterWorkDir hb_outstem], MATLAB_DIR);
        fprintf(script_hb, 'Main_RNASpots(');

        printMatArg(script_hb, 'imgname', iname, false);
        ipath = getTableValue(image_table, r, 'IMAGEPATH');
        if endsWith(ipath, '.mat')
            printMatArg(script_hb, 'matimg', [ClusterWorkDir ipath], true);
        elseif endsWith(ipath, '.tif')
            printMatArg(script_hb, 'tif', [ClusterWorkDir ipath], true);
        end

        [hb_outdir, ~, ~] = fileparts(hb_outstem);
        printMatArg(script_hb, 'outdir', [ClusterWorkDir hb_outdir], true);
        csegdir = getTableValue(image_table, r, 'CELLSEG_DIR');
        csegsfx = getTableValue(image_table, r, 'CELLSEG_SFX');
        if ~strcmp(csegdir, '.')
            printMatArg(script_hb, 'cellseg', [ClusterWorkDir csegdir '/Lab_' csegsfx '.mat'], true);
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

        if DETECT_THREADS > 1
            printMatArg(script_hb, 'threads', num2str(DETECT_THREADS), true);
        end

        if startsWith(iname, 'sim_')
            printMatFlagArg(script_hb, 'nodpc', true);
        end
        printMatFlagArg(script_hb, 'debug', true);

        fprintf(script_hb, '); quit;"\n');
        fclose(script_hb);
    end

%BF Script
    bfstem = getTableValue(image_table, r, 'BIGFISH_OUTSTEM');
    [bfoutdir, ~, ~] = fileparts(bfstem);
    if RUN_BF
        script_bf = fopen([ScriptDir filesep iname '_bf.sh'], 'w');
        fprintf(script_bf, '#!/bin/bash\n\n');
        fprintf(script_bf, 'source %s/bigfish/bin/activate\n', ClusterPyenvDir);
        fprintf(script_bf, 'python3 %s/bigfish_wrapper.py "%s" "%s"', ClusterScriptsDir, [ClusterWorkDir ipath], [ClusterWorkDir bfoutdir]);
    
        %- BF Options
        fprintf(script_bf, ' --ch_dapi %d', getTableValue(image_table, r, 'CH_DAPI'));
        fprintf(script_bf, ' --ch_light %d', getTableValue(image_table, r, 'CH_LIGHT'));
        fprintf(script_bf, ' --ch_target %d', getTableValue(image_table, r, 'CHANNEL'));
        fprintf(script_bf, ' --minth 10');
        fprintf(script_bf, ' --maxth 1000');
        
        x = getTableValue(image_table, r, 'VOXEL_X');
        y = getTableValue(image_table, r, 'VOXEL_Y');
        z = getTableValue(image_table, r, 'VOXEL_Z');
        fprintf(script_bf, ' --voxelsz "(%d,%d,%d)"', z,y,x);
        x = getTableValue(image_table, r, 'POINT_X');
        y = getTableValue(image_table, r, 'POINT_Y');
    	z = getTableValue(image_table, r, 'POINT_Z');
        fprintf(script_bf, ' --pointsz "(%d,%d,%d)"', z,y,x);

        fprintf(script_bf, ' --sobjsznuc %d', BF_SOBJSZ);
        fprintf(script_bf, ' --trgsznuc %d', BF_NUCSZ);

        if startsWith(iname, 'sim_')
            fprintf(script_bf, ' --zkeep 1.0');
        else
            fprintf(script_bf, ' --zkeep 0.8');
        end

        if ~BF_RESCALE
            fprintf(script_bf, ' --norescale');
        end
        fprintf(script_bf, '\ndeactivate\n\n');
        fprintf(script_bf, 'module load %s\n', MODULE_NAME);
        fprintf(script_bf, 'cd %s\n', MATLAB_DIR);
        fprintf(script_bf, 'matlab -nodisplay -nosplash -logfile "%s" -r "cd %s; Main_Bigfish2Mat(''%s'',''%s''); quit;"\n',...
            [ClusterWorkDir bfstem '_mat.log'], MATLAB_DIR, [ClusterWorkDir bfoutdir], [ClusterWorkDir bfstem]);
        fprintf(script_bf, 'module unload\n\n');

        fprintf(script_bf, 'if [ -s "%s_coordTable.mat" ]; then\n', [ClusterWorkDir bfstem]);
        fprintf(script_bf, '\techo -e "Coord table file found. Assuming conversion success."\n');
        fprintf(script_bf, '\tcd "%s"\n', [ClusterWorkDir bfoutdir]);
        fprintf(script_bf, '\trm ./spots_*.csv\n');
        fprintf(script_bf, 'fi\n\n');

        fclose(script_bf);
    end

    %When adding master script lines, make sure to check if image is there
    %before submitting job...
    
    fprintf(script_master, 'if [ -s "%s" ]; then\n', [ClusterWorkDir ipath]);
    
    if RUN_HB
        fprintf(script_master, '\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_hb.sh']);
        fprintf(script_master, '\tmkdir -p "%s"\n', [ClusterWorkDir hb_outdir]);
        fprintf(script_master, '\tsbatch');
        fprintf(script_master, ' --job-name="%s"', ['SpotDetect_' iname]);
        if DETECT_THREADS > 1
            fprintf(script_master, ' --cpus-per-task=%d', DETECT_THREADS);
            fprintf(script_master, ' --time=4:00:00');
            fprintf(script_master, ' --mem=%dg', (RAM_PER_CORE * DETECT_THREADS));
        else
            fprintf(script_master, ' --cpus-per-task=2');
            fprintf(script_master, ' --time=8:00:00');
            fprintf(script_master, ' --mem=%dg', RAM_PER_CORE);
        end
        fprintf(script_master, ' --error="%s"', [ClusterWorkDir hb_outstem '_slurm.err']);
        fprintf(script_master, ' --output="%s"', [ClusterWorkDir hb_outstem '_slurm.out']);
        fprintf(script_master, ' "${SCRIPTDIR}/%s"\n', [iname '_hb.sh']);
    end
    if RUN_BF
        fprintf(script_master, '\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_bf.sh']);
        fprintf(script_master, '\tmkdir -p "%s"\n', [ClusterWorkDir bfoutdir]);
        fprintf(script_master, '\tsbatch');
        fprintf(script_master, ' --job-name="%s"', ['Bigfish_' iname]);
        fprintf(script_master, ' --cpus-per-task=4');
        fprintf(script_master, ' --time=8:00:00');
        fprintf(script_master, ' --mem=%dg', RAM_PER_CORE);
        fprintf(script_master, ' --error="%s"', [ClusterWorkDir bfstem '_slurm.err']);
        fprintf(script_master, ' --output="%s"', [ClusterWorkDir bfstem '_slurm.out']);
        fprintf(script_master, ' "${SCRIPTDIR}/%s"\n', [iname '_bf.sh']);
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