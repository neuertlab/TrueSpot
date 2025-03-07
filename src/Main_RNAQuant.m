%
%%

function Main_RNAQuant(varargin)
addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2025.03.07.00';
VERSION_STRING = 'v1.1.2';

DEFAULT_PRESET_INDEX = 6;
MAX_TH_PRESET_LEVEL = 5;

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.
param_struct.runpath = [];
param_struct.tifpath = [];
param_struct.outdir = [];
param_struct.cellsegdir = [];
param_struct.cellsegname = [];
param_struct.cellsegpath = [];
param_struct.nucsegpath = [];
param_struct.coord_tbl_path = [];
param_struct.rethresh = false;
param_struct.rethreshPreset = 0;
param_struct.noclouds = false;
param_struct.no_refilter = false;
param_struct.dbgcell = 0;
param_struct.no_bkg_subtract = false;
param_struct.use_nuc_mask = 2; %lbl_mid
param_struct.workers = 1;

param_struct.nocells = false; %If no cell seg data (ie. sim image)

%Specs if don't provide run file
param_struct.man_thresh = 0;
param_struct.rna_channel = 1;
param_struct.channel_count = 1;
param_struct.z_adj = 1.0;

param_struct.small_obj_size = 3;
param_struct.gaussian_radius = 0; %Defaults to 7 or spotsrun val if no override.
param_struct.spotzoom_r_xy = 4;
param_struct.spotzoom_r_z = 2;

param_struct.preloaded_image = [];
param_struct.preloaded_cellmask = [];
param_struct.preloaded_nucmask = [];
param_struct.coords = [];

senspe_set = false;

lastkey = [];
for i = 1:nargin
    argval = varargin{i};
    
    if ischar(argval) & startsWith(argval, "-")
        %Key
        if size(argval,2) >= 2
            lastkey = argval(2:end);
        else
            lastkey = [];
        end
        
        if strcmp(lastkey, "rethresh")
            param_struct.rethresh = true;
            if arg_debug; fprintf("Rerun Thresholder: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "noclouds")
            param_struct.noclouds = true;
            if arg_debug; fprintf("Detect Clouds: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "norefilter")
            param_struct.no_refilter = true;
            if arg_debug; fprintf("Refilter: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "nmhi")
            param_struct.use_nuc_mask = 3;
            if arg_debug; fprintf("Using high nuc mask.\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "nmlo")
            param_struct.use_nuc_mask = 1;
            if arg_debug; fprintf("Using low nuc mask.\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "nm2d")
            param_struct.use_nuc_mask = 0;
            if arg_debug; fprintf("Using 2D nuc mask.\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "nocells")
            param_struct.nocells = true;
            if arg_debug; fprintf("Bypassing cell segmentation data.\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "sensitive")
            if ~senspe_set
                if arg_debug; fprintf("Rethreshold tuning Preset: Sensitivity\n"); end
                senspe_set = true;
                param_struct.rethreshPreset = DEFAULT_PRESET_INDEX + 3;
                lastkey = [];
            end
        elseif strcmp(lastkey, "precise")
            if ~senspe_set
                if arg_debug; fprintf("Rethreshold tuning Preset: Precise\n"); end
                senspe_set = true;
                param_struct.rethreshPreset = DEFAULT_PRESET_INDEX - 3;
                lastkey = [];
            end
        elseif strcmp(lastkey, "nosubbkg")
            param_struct.no_bkg_subtract = true;
            if arg_debug; fprintf("Calculate & Subtract background: Off\n"); end
            lastkey = [];
        end
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "runinfo")
            param_struct.runpath = argval;
            if arg_debug; fprintf("Spots Run Info Path Set: %s\n", param_struct.runpath); end
        elseif strcmp(lastkey, "tif")
            param_struct.tifpath = argval;
            if arg_debug; fprintf("TIF Path Set: %s\n", param_struct.tifpath); end
        elseif strcmp(lastkey, "outdir")
            param_struct.outdir = argval;
            if arg_debug; fprintf("Output Directory Set: %s\n", param_struct.outdir); end
        elseif strcmp(lastkey, "cellsegdir")
            param_struct.cellsegdir = argval;
            if arg_debug; fprintf("CellSeg Directory Set: %s\n", param_struct.cellsegdir); end
        elseif strcmp(lastkey, "cellsegname")
            param_struct.cellsegname = argval;
            if arg_debug; fprintf("CellSeg Name Set: %s\n", param_struct.cellsegname); end
        elseif strcmp(lastkey, "cellsegpath")
            param_struct.cellsegpath = argval;
            if arg_debug; fprintf("CellSeg Path Set: %s\n", param_struct.cellsegpath); end
        elseif strcmp(lastkey, "nucsegpath")
            param_struct.nucsegpath = argval;
            if arg_debug; fprintf("NucSeg Path Set: %s\n", param_struct.nucsegpath); end
        elseif strcmp(lastkey, "mthresh")
            param_struct.man_thresh = Force2Num(argval);
            if arg_debug; fprintf("Manual threshold set: %d\n", param_struct.man_thresh); end
        elseif strcmp(lastkey, "coordtable")
            param_struct.coord_tbl_path = argval;
            if arg_debug; fprintf("Manual coord table path set: %s\n", param_struct.coord_tbl_path); end
        elseif strcmp(lastkey, "ch")
            param_struct.rna_channel = Force2Num(argval);
            if arg_debug; fprintf("Signal channel set: %d\n", param_struct.rna_channel); end
        elseif strcmp(lastkey, "chcount")
            param_struct.channel_count = Force2Num(argval);
            if arg_debug; fprintf("Channel count set: %d\n", param_struct.channel_count); end
        elseif strcmp(lastkey, "zadj")
            param_struct.z_adj = Force2Num(argval);
            if arg_debug; fprintf("Z to XY voxel dim ratio set: %f\n", param_struct.z_adj); end
        elseif strcmp(lastkey, "smobjsz")
            param_struct.small_obj_size = Force2Num(argval);
            if arg_debug; fprintf("Small object size set: %d\n", param_struct.small_obj_size); end
        elseif strcmp(lastkey, "gaussrad")
            param_struct.gaussian_radius = Force2Num(argval);
            if arg_debug; fprintf("Gaussian radius set: %d\n", param_struct.gaussian_radius); end
        elseif strcmp(lastkey, "radxy")
            param_struct.spotzoom_r_xy = Force2Num(argval);
            if arg_debug; fprintf("Sampling xy radius set: %d\n", param_struct.spotzoom_r_xy); end
        elseif strcmp(lastkey, "radz")
            param_struct.spotzoom_r_z = Force2Num(argval);
            if arg_debug; fprintf("Sampling z radius set: %d\n", param_struct.spotzoom_r_z); end
        elseif strcmp(lastkey, "workers")
            param_struct.workers = Force2Num(argval);
            if arg_debug; fprintf("Worker threads/processes requested: %d\n", param_struct.workers); end
        elseif strcmp(lastkey, "dbgcell")
            param_struct.dbgcell = Force2Num(argval);
            if arg_debug; fprintf("DEBUG Cell Index: %d\n", param_struct.dbgcell); end
        elseif strcmp(lastkey, "sensitivity")
            specval = Force2Num(argval);
            if specval > MAX_TH_PRESET_LEVEL; specval = MAX_TH_PRESET_LEVEL; end
            param_struct.rethreshPreset = DEFAULT_PRESET_INDEX + specval;
            if arg_debug; fprintf("Rethreshold sensitivity Preset Level Set: %d\n", specval); end
        elseif strcmp(lastkey, "precision")
            specval = Force2Num(argval);
            if specval > MAX_TH_PRESET_LEVEL; specval = MAX_TH_PRESET_LEVEL; end
            param_struct.rethreshPreset = DEFAULT_PRESET_INDEX - specval;
            if arg_debug; fprintf("Rethreshold precision Preset Level Set: %d\n", specval); end
        end
    end
end %End of argin for loop

% ========================== Run ==========================

if isempty(param_struct.outdir)
%     fprintf("Output directory is required for command line interface! Exiting... \n");
%     return;
%Set to run dir, or tif dir if run dir is not provided.
    if ~isempty(param_struct.runpath)
        [param_struct.outdir, ~, ~] = fileparts(param_struct.runpath);
    elseif ~isempty(param_struct.tifpath)
        [param_struct.outdir, ~, ~] = fileparts(param_struct.tifpath);
    else
        fprintf("Output directory is required for command line interface! Exiting... \n");
        return;
    end
    fprintf("WARNING: Output directory was not provided. Set to %s \n", param_struct.outdir);
end

quant_results = RNAQuantPipe(param_struct, false);

% ========================== Save Results ==========================
if isempty(quant_results)
    fprintf("Quantification failed! See log for details. Exiting...\n");
    return;
end

%Clean up unneeded input fields.
quant_results.img_raw = [];
quant_results.cell_mask = [];
quant_results.nuc_mask = [];
quant_results.t_coord_table = [];

runMeta = struct();
runMeta.modifiedDate = datetime;
runMeta.tsQuantBuild = BUILD_STRING;
runMeta.tsQuantVersion = VERSION_STRING;

imgname = guessImageName(param_struct);

outpath = [param_struct.outdir filesep imgname '_quantResults.csv'];
quant_table = RNAQuant.cellData2Table(quant_results.cell_rna_data);
writetable(quant_table, outpath);

outpath = [param_struct.outdir filesep imgname '_quantData.mat'];
quant_results = RNAQuant.results2SavePackage(quant_results);
save(outpath, 'quant_results', 'runMeta');
%save(outpath, 'quant_results', 'runMeta', '-v7.3');

end %end of Main function

function iname = guessImageName(param_struct)
    %Do we have a valid run?
    if ~isempty(param_struct.runpath)
        rnaspots_run = RNASpotsRun.loadFrom(param_struct.runpath);
        if ~isempty(rnaspots_run)
            if ~isempty(rnaspots_run.img_name)
                iname = rnaspots_run.img_name;
                return;
            end
        end
    end
    
    %Just use tif name.
    if ~isempty(param_struct.tifpath)
        [~, iname, ~] = fileparts(param_struct.tifpath);
        return;
    end
    
    iname = 'my_image';
end
