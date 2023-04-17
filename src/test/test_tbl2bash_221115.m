%
%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

%DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
DataDir = 'D:\usr\bghos\labdat\imgproc';

%ScriptDir = 'C:\Users\hospelb.VUDS\Desktop\slurm';
ScriptDir = 'C:\Users\bghos\Desktop\slurm';

ClusterWorkDir = '/nobackup/p_neuert_lab/hospelb/imgproc';
ClusterSlurmDir = '/nobackup/p_neuert_lab/hospelb/imgproc/slurm/script';
ClusterScriptsDir = '/nobackup/p_neuert_lab/hospelb/scripts';
ClusterPyenvDir = '/home/hospelb/pyvenv';

% ========================== Constants ==========================

DETECT_THREADS = 8;
RAM_PER_CORE = 8;
RAM_PER_CORE_BF = 64;
HB_PARALLEL_HR = 6;
HB_SERIAL_HR = 10;
BF_SERIAL_HR = 12;

QUANT_SERIAL_HR = 6;
QUANT_PARALLEL_HR = 4;
RAM_PER_CORE_QUANT = 20;
DETECT_THREADS_QUANT = 1;
QUANT_DO_CLOUDS = false;
QUANT_FIXED_TH = 0;

TH_MIN = 0;
TH_MAX = 0;
TH_MIN_BF = 10;
Z_TRIM = 0;
BF_SOBJSZ = 10;
BF_NUCSZ = 256; %200 yeast, 256 mesc
%BF_RESCALE = false;

RUN_HB = true;
RUN_BFNR = false;
RUN_BFRS = true;
RUN_QUANT = true;
OVERWRITE = false;

MODULE_NAME = 'MATLAB/2018b';
MATLAB_DIR = [ClusterWorkDir '/matlab'];
PYVENV_NAME = 'bigfish';

% ========================== Load csv Table ==========================
%InputTablePath = [DataDir filesep 'test_images_simytc.csv'];
%InputTablePath = [DataDir filesep 'test_images_simvarmass.csv'];
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

%ImageName='scrna_E2R2I5_CTT1';
GroupPrefix = 'histonesc_D2';
GroupSuffix = [];
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

if RUN_QUANT
    script_master_quant = fopen([ScriptDir filesep 'runall_q.sh'], 'w');
    fprintf(script_master_quant, '#!/bin/bash\n\n');
    fprintf(script_master_quant, 'SCRIPTDIR=%s\n', ClusterSlurmDir);
end

rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if ~startsWith(iname, GroupPrefix); continue; end
    if ~isempty(GroupSuffix)
        if ~endsWith(iname, GroupSuffix); continue; end
    end
    hb_outstem = getTableValue(image_table, r, 'OUTSTEM');
    
    ipath = getTableValue(image_table, r, 'IMAGEPATH');
    if endsWith(ipath, '.mat')
        ipath = replace(ipath, '.mat', '.tif');
    end
    
    %HB Script
    if RUN_HB
        script_hb = fopen([ScriptDir filesep iname '_hb.sh'], 'w');
        fprintf(script_hb, '#!/bin/bash\n\n');
        fprintf(script_hb, 'module load %s\n', MODULE_NAME);
        fprintf(script_hb, 'cd %s\n', MATLAB_DIR);
        fprintf(script_hb, 'matlab -nodisplay -nosplash -logfile "%s_mat.log" -r "cd %s; ', [ClusterWorkDir hb_outstem], MATLAB_DIR);
        fprintf(script_hb, 'Main_RNASpots(');

        printMatArg(script_hb, 'imgname', iname, false);
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

        if TH_MIN > 0
            printMatArg(script_hb, 'thmin', num2str(TH_MIN), true);
        else
            printMatFlagArg(script_hb, 'autominth', true);
        end
        if TH_MAX > 0
            printMatArg(script_hb, 'thmax', num2str(TH_MAX), true);
        else
            printMatFlagArg(script_hb, 'automaxth', true);
        end
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
    vx = 0; vy = 0; vz = 0;
    px = 0; py = 0; pz = 0;
    if RUN_BFNR
        script_bfnr = fopen([ScriptDir filesep iname '_bfnr.sh'], 'w');
        fprintf(script_bfnr, '#!/bin/bash\n\n');
        fprintf(script_bfnr, 'module load GCC/6.4.0-2.28\n');
        fprintf(script_bfnr, 'module load Intel/2017.4.196\n');
        fprintf(script_bfnr, 'module load Python/3.6.3\n');
        fprintf(script_bfnr, 'source %s/%s/bin/activate\n', ClusterPyenvDir, PYVENV_NAME);
        fprintf(script_bfnr, 'python3 %s/bigfish_wrapper.py "%s" "%s"', ClusterScriptsDir, [ClusterWorkDir ipath], [ClusterWorkDir bfoutdir]);
    
        %- BF Options
        fprintf(script_bfnr, ' --ch_dapi %d', getTableValue(image_table, r, 'CH_DAPI'));
        fprintf(script_bfnr, ' --ch_light %d', getTableValue(image_table, r, 'CH_LIGHT'));
        fprintf(script_bfnr, ' --ch_target %d', getTableValue(image_table, r, 'CHANNEL'));
        fprintf(script_bfnr, ' --minth %d', TH_MIN_BF);
        fprintf(script_bfnr, ' --maxth 1000');
        
        vx = getTableValue(image_table, r, 'VOXEL_X');
        vy = getTableValue(image_table, r, 'VOXEL_Y');
        vz = getTableValue(image_table, r, 'VOXEL_Z');
        fprintf(script_bfnr, ' --voxelsz "(%d,%d,%d)"', vz,vy,vx);
        px = getTableValue(image_table, r, 'POINT_X');
        py = getTableValue(image_table, r, 'POINT_Y');
    	pz = getTableValue(image_table, r, 'POINT_Z');
        
        %Big-FISH has been acting wierd if point size is not >= vox size?
        %I don't know if it is supposed to do that, but I'm sick of it
        %crashing.
        px = max(vx,px);
        py = max(vy,py);
        pz = max(vz,pz);
        fprintf(script_bfnr, ' --pointsz "(%d,%d,%d)"', pz,py,px);

        fprintf(script_bfnr, ' --sobjsznuc %d', BF_SOBJSZ);
        fprintf(script_bfnr, ' --trgsznuc %d', BF_NUCSZ);

        if startsWith(iname, 'sim_') | (Z_TRIM < 1)
            fprintf(script_bfnr, ' --zkeep 1.0');
        else
            fprintf(script_bfnr, ' --zkeep 0.8');
        end

        fprintf(script_bfnr, ' --norescale');
        
        if RUN_QUANT
            fprintf(script_bfnr, ' --gaussfit');
        end

        fprintf(script_bfnr, '\ndeactivate\n\n');
        fprintf(script_bfnr, 'module load %s\n', MODULE_NAME);
        fprintf(script_bfnr, 'cd %s\n', MATLAB_DIR);
        fprintf(script_bfnr, 'matlab -nodisplay -nosplash -logfile "%s" -r "cd %s; Main_Bigfish2Mat(''%s'',''%s''); quit;"\n',...
            [ClusterWorkDir bfstem '_mat.log'], MATLAB_DIR, [ClusterWorkDir bfoutdir], [ClusterWorkDir bfstem]);
        fprintf(script_bfnr, 'module unload\n\n');

        fprintf(script_bfnr, 'if [ -s "%s_coordTable.mat" ]; then\n', [ClusterWorkDir bfstem]);
        fprintf(script_bfnr, '\techo -e "Coord table file found. Assuming conversion success."\n');
        fprintf(script_bfnr, '\tcd "%s"\n', [ClusterWorkDir bfoutdir]);
        fprintf(script_bfnr, '\trm ./spots_*.csv\n');
        fprintf(script_bfnr, 'fi\n\n');

        fclose(script_bfnr);
    end
    
    bfrs_stem = replace(bfstem, '/data/bigfish/', '/data/bigfish/_rescaled/');
    [bfrs_outdir, ~, ~] = fileparts(bfrs_stem);
    if RUN_BFRS
        script_bfrs = fopen([ScriptDir filesep iname '_bfrs.sh'], 'w');
        fprintf(script_bfrs, '#!/bin/bash\n\n');
        fprintf(script_bfrs, 'module load GCC/6.4.0-2.28\n');
        fprintf(script_bfrs, 'module load Intel/2017.4.196\n');
        fprintf(script_bfrs, 'module load Python/3.6.3\n');
        fprintf(script_bfrs, 'source %s/%s/bin/activate\n', ClusterPyenvDir, PYVENV_NAME);
        fprintf(script_bfrs, 'python3 %s/bigfish_wrapper.py "%s" "%s"', ClusterScriptsDir, [ClusterWorkDir ipath], [ClusterWorkDir bfrs_outdir]);
    
        %- BF Options
        fprintf(script_bfrs, ' --ch_dapi %d', getTableValue(image_table, r, 'CH_DAPI'));
        fprintf(script_bfrs, ' --ch_light %d', getTableValue(image_table, r, 'CH_LIGHT'));
        fprintf(script_bfrs, ' --ch_target %d', getTableValue(image_table, r, 'CHANNEL'));
        fprintf(script_bfrs, ' --minth 10');
        fprintf(script_bfrs, ' --maxth 1000');
        
        vx = getTableValue(image_table, r, 'VOXEL_X');
        vy = getTableValue(image_table, r, 'VOXEL_Y');
        vz = getTableValue(image_table, r, 'VOXEL_Z');
        fprintf(script_bfrs, ' --voxelsz "(%d,%d,%d)"', vz,vy,vx);
        px = getTableValue(image_table, r, 'POINT_X');
        py = getTableValue(image_table, r, 'POINT_Y');
    	pz = getTableValue(image_table, r, 'POINT_Z');
        
        px = max(vx,px);
        py = max(vy,py);
        pz = max(vz,pz);
        fprintf(script_bfrs, ' --pointsz "(%d,%d,%d)"', pz,py,px);

        fprintf(script_bfrs, ' --sobjsznuc %d', BF_SOBJSZ);
        fprintf(script_bfrs, ' --trgsznuc %d', BF_NUCSZ);

        if Z_TRIM <= 0
            fprintf(script_bfrs, ' --zkeep 1.0');
        else
            fprintf(script_bfrs, ' --zkeep 0.8');
        end

        if RUN_QUANT
            fprintf(script_bfrs, ' --gaussfit');
        end

        fprintf(script_bfrs, '\ndeactivate\n\n');
        fprintf(script_bfrs, 'module load %s\n', MODULE_NAME);
        fprintf(script_bfrs, 'cd %s\n', MATLAB_DIR);
        fprintf(script_bfrs, 'matlab -nodisplay -nosplash -logfile "%s" -r "cd %s; Main_Bigfish2Mat(''%s'',''%s''); quit;"\n',...
            [ClusterWorkDir bfrs_stem '_mat.log'], MATLAB_DIR, [ClusterWorkDir bfrs_outdir], [ClusterWorkDir bfrs_stem]);
        fprintf(script_bfrs, 'module unload\n\n');

        fprintf(script_bfrs, 'if [ -s "%s_coordTable.mat" ]; then\n', [ClusterWorkDir bfrs_stem]);
        fprintf(script_bfrs, '\techo -e "Coord table file found. Assuming conversion success."\n');
        fprintf(script_bfrs, '\tcd "%s"\n', [ClusterWorkDir bfrs_outdir]);
        fprintf(script_bfrs, '\trm ./spots_*.csv\n');
        fprintf(script_bfrs, 'fi\n\n');

        fclose(script_bfrs);
    end
    
    if RUN_QUANT
        script_qt = fopen([ScriptDir filesep iname '_qt.sh'], 'w');
        fprintf(script_qt, '#!/bin/bash\n\n');
        fprintf(script_qt, 'module load %s\n', MODULE_NAME);
        fprintf(script_qt, 'cd %s\n', MATLAB_DIR);
        fprintf(script_qt, 'matlab -nodisplay -nosplash -logfile "%s_quant_mat.log" -r "cd %s; ', [ClusterWorkDir hb_outstem], MATLAB_DIR);
        fprintf(script_qt, 'Main_RNAQuant(');
        
        printMatArg(script_qt, 'runinfo', [ClusterWorkDir hb_outstem '_rnaspotsrun.mat'], false);
        printMatArg(script_qt, 'tif', [ClusterWorkDir ipath], true);
        
        [hb_outdir, ~, ~] = fileparts(hb_outstem);
        printMatArg(script_qt, 'outdir', [ClusterWorkDir hb_outdir], true);
        
        csegdir = getTableValue(image_table, r, 'CELLSEG_DIR');
        csegsfx = getTableValue(image_table, r, 'CELLSEG_SFX');

        %If there is no cellseg data, be sure to use the nocells arg
        if strcmp(csegdir, '.') | strcmp(csegsfx, '.')
            printMatFlagArg(script_qt, 'nocells', true);
        else
            printMatArg(script_qt, 'cellsegdir', [ClusterWorkDir csegdir], true);
            printMatArg(script_qt, 'cellsegname', csegsfx, true);
        end

        if QUANT_FIXED_TH > 0
            printMatArg(script_qt, 'mthresh', num2str(QUANT_FIXED_TH), true);
        end
        
        if DETECT_THREADS_QUANT > 1
            printMatArg(script_qt, 'workers', num2str(DETECT_THREADS_QUANT), true);
        end
        
        if ~QUANT_DO_CLOUDS
            printMatFlagArg(script_qt, 'noclouds', true);
        end
        
        fprintf(script_qt, '); quit;"\n');
        fclose(script_qt);
    end

    %When adding master script lines, make sure to check if image is there
    %before submitting job...
    
    fprintf(script_master, 'if [ -s "%s" ]; then\n', [ClusterWorkDir ipath]);
    
    if RUN_HB
        if ~OVERWRITE
            fprintf(script_master, '\tif [ ! -s "%s" ]; then\n', [ClusterWorkDir hb_outdir '/' iname '_all_3d_coordTable.mat']);
            fprintf(script_master, '\t\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_hb.sh']);
            fprintf(script_master, '\t\tmkdir -p "%s"\n', [ClusterWorkDir hb_outdir]);
            fprintf(script_master, '\t\tsbatch');
        else
            fprintf(script_master, '\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_hb.sh']);
            fprintf(script_master, '\tmkdir -p "%s"\n', [ClusterWorkDir hb_outdir]);
            fprintf(script_master, '\tsbatch');
        end
        
        fprintf(script_master, ' --job-name="%s"', ['SpotDetect_' iname]);
        if DETECT_THREADS > 1
            fprintf(script_master, ' --cpus-per-task=%d', DETECT_THREADS);
            fprintf(script_master, ' --time=%d:00:00', HB_PARALLEL_HR);
            fprintf(script_master, ' --mem=%dg', (RAM_PER_CORE * DETECT_THREADS));
        else
            fprintf(script_master, ' --cpus-per-task=2');
            fprintf(script_master, ' --time=%d:00:00', HB_SERIAL_HR);
            fprintf(script_master, ' --mem=%dg', RAM_PER_CORE);
        end
        fprintf(script_master, ' --error="%s"', [ClusterWorkDir hb_outstem '_slurm.err']);
        fprintf(script_master, ' --output="%s"', [ClusterWorkDir hb_outstem '_slurm.out']);
        fprintf(script_master, ' "${SCRIPTDIR}/%s"\n', [iname '_hb.sh']);
        
        if ~OVERWRITE
            fprintf(script_master, '\telse\n');
            fprintf(script_master, '\t\techo -e "HB run for %s found! Not resubmitting..."\n', iname);
            fprintf(script_master, '\tfi\n');
        end
    end
    if RUN_BFNR
        if ~OVERWRITE
            fprintf(script_master, '\tif [ ! -s "%s" ]; then\n', [ClusterWorkDir bfstem '_coordTable.mat']);
            fprintf(script_master, '\t\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_bfnr.sh']);
            fprintf(script_master, '\t\tmkdir -p "%s"\n', [ClusterWorkDir bfoutdir]);
            fprintf(script_master, '\t\tsbatch');
        else
            fprintf(script_master, '\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_bfnr.sh']);
            fprintf(script_master, '\tmkdir -p "%s"\n', [ClusterWorkDir bfoutdir]);
            fprintf(script_master, '\tsbatch');
        end
        
        fprintf(script_master, ' --job-name="%s"', ['Bigfish_' iname]);
        fprintf(script_master, ' --cpus-per-task=4');
        fprintf(script_master, ' --time=%d:00:00', BF_SERIAL_HR);
        fprintf(script_master, ' --mem=%dg', RAM_PER_CORE_BF);
        fprintf(script_master, ' --error="%s"', [ClusterWorkDir bfstem '_slurm.err']);
        fprintf(script_master, ' --output="%s"', [ClusterWorkDir bfstem '_slurm.out']);
        fprintf(script_master, ' "${SCRIPTDIR}/%s"\n', [iname '_bfnr.sh']);
        
        if ~OVERWRITE
            fprintf(script_master, '\telse\n');
            fprintf(script_master, '\t\techo -e "BFNR run for %s found! Not resubmitting..."\n', iname);
            fprintf(script_master, '\tfi\n');
        end
    end
    if RUN_BFRS
        if ~OVERWRITE
            fprintf(script_master, '\tif [ ! -s "%s" ]; then\n', [ClusterWorkDir bfrs_stem '_coordTable.mat']);
            fprintf(script_master, '\t\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_bfrs.sh']);
            fprintf(script_master, '\t\tmkdir -p "%s"\n', [ClusterWorkDir bfrs_outdir]);
            fprintf(script_master, '\t\tsbatch');
        else
            fprintf(script_master, '\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_bfrs.sh']);
            fprintf(script_master, '\tmkdir -p "%s"\n', [ClusterWorkDir bfrs_outdir]);
            fprintf(script_master, '\tsbatch');
        end
        
        fprintf(script_master, ' --job-name="%s"', ['BigfishRS_' iname]);
        fprintf(script_master, ' --cpus-per-task=4');
        fprintf(script_master, ' --time=%d:00:00', BF_SERIAL_HR);
        fprintf(script_master, ' --mem=%dg', RAM_PER_CORE_BF);
        fprintf(script_master, ' --error="%s"', [ClusterWorkDir bfrs_stem '_slurm.err']);
        fprintf(script_master, ' --output="%s"', [ClusterWorkDir bfrs_stem '_slurm.out']);
        fprintf(script_master, ' "${SCRIPTDIR}/%s"\n', [iname '_bfrs.sh']);
        
        if ~OVERWRITE
            fprintf(script_master, '\telse\n');
            fprintf(script_master, '\t\techo -e "BFRS run for %s found! Not resubmitting..."\n', iname);
            fprintf(script_master, '\tfi\n');
        end
    end

    fprintf(script_master, 'fi\n\n');
    
    if RUN_QUANT
        fprintf(script_master_quant, 'if [ -s "%s" ]; then\n', [ClusterWorkDir ipath]);
        fprintf(script_master_quant, '\tif [ -s "%s" ]; then\n', [ClusterWorkDir hb_outstem '_rnaspotsrun.mat']);
        fprintf(script_master_quant, '\t\techo -e "Requested run found: %s"\n', hb_outstem);
        
        fprintf(script_master_quant, '\t\tchmod 770 "${SCRIPTDIR}/%s"\n', [iname '_qt.sh']);
        fprintf(script_master_quant, '\t\tsbatch');
        fprintf(script_master_quant, ' --job-name="%s"', ['RNAQuant_' iname]);
        
        if DETECT_THREADS_QUANT > 1
            fprintf(script_master_quant, ' --cpus-per-task=%d', DETECT_THREADS_QUANT);
            fprintf(script_master_quant, ' --time=%d:00:00', QUANT_PARALLEL_HR);
            fprintf(script_master_quant, ' --mem=%dg', (RAM_PER_CORE_QUANT * DETECT_THREADS_QUANT));
        else
            fprintf(script_master_quant, ' --cpus-per-task=2');
            fprintf(script_master_quant, ' --time=%d:00:00', QUANT_SERIAL_HR);
            fprintf(script_master_quant, ' --mem=%dg', RAM_PER_CORE_QUANT);
        end
        fprintf(script_master_quant, ' --error="%s"', [ClusterWorkDir hb_outstem '_quant_slurm.err']);
        fprintf(script_master_quant, ' --output="%s"', [ClusterWorkDir hb_outstem '_quant_slurm.out']);
        fprintf(script_master_quant, ' "${SCRIPTDIR}/%s"\n', [iname '_qt.sh']);
        
        fprintf(script_master_quant, '\tfi\n');
        fprintf(script_master_quant, 'fi\n\n');
    end
    
end
fclose(script_master);

if RUN_QUANT
    fclose(script_master_quant);
end

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