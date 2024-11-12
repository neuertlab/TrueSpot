%%
%%

function rna_spot_run = Main_RNASpots(varargin)
addpath('./core');
addpath('./thirdparty');
addpath('./celldissect');

DEFAULT_PRESET_INDEX = 6;
MAX_TH_PRESET_LEVEL = 5;

% ========================== Process args ==========================
rna_spot_run = RNASpotsRun.initDefault();
rna_spot_run = setThPreset(rna_spot_run, DEFAULT_PRESET_INDEX, 0);

%Highest level set to true is what is used, regardless of other flags.
arg_debug = true; %CONSTANT used for debugging arg parser.
debug_lvl = 0;
senspe_set = false;
matvar = [];
thread_request = 1;
automaxth_flag = false;
autominth_flag = false;
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
        
        %Account for boolean keys...
        if strcmp(lastkey, "ovrw")
            rna_spot_run.options.overwrite_output = true;
            if arg_debug; fprintf("Overwrite Output: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "nodpc")
            rna_spot_run.options.deadpix_detect = false;
            if arg_debug; fprintf("Dead Pixel Cleaning: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "usespc")
            rna_spot_run.th_params.test_data = true;
            if arg_debug; fprintf("Use Raw Spot Count Curve: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "usedfc")
            rna_spot_run.th_params.test_diff = true;
            if arg_debug; fprintf("Use Diff Curve: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "nospline")
            rna_spot_run.th_params.spline_iterations = 0;
            if arg_debug; fprintf("Use Spline: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "fitsegrsq")
            rna_spot_run.th_params.reweight_fit = false;
            if arg_debug; fprintf("Use weighted segmented fit: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "fitwavg")
            rna_spot_run.th_params.reweight_fit = true;
            if arg_debug; fprintf("Use weighted segmented fit: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "verbose")
            if debug_lvl < 1; debug_lvl = 1; end
            if arg_debug; fprintf("Verbose: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "quiet")
            if arg_debug; fprintf("Verbose: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "plotout")
            if debug_lvl < 2; debug_lvl = 2; end
            if arg_debug; fprintf("Output Plots: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "debug")
            if debug_lvl < 2; debug_lvl = 2; end
            if arg_debug; fprintf("Debug Mode: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "debugv")
            if debug_lvl < 3; debug_lvl = 3; end
            if arg_debug; fprintf("Verbose Debug Mode: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "csvzero")
            rna_spot_run.options.csv_zero_based_coords = true;
            if arg_debug; fprintf("csv Output Use 0-Based Coords: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "csvfull")
            rna_spot_run.options.csv_output_level = 1;
            if arg_debug; fprintf("csv Output Level: All Thresholds\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "csvrange")
            rna_spot_run.options.csv_output_level = 2;
            if arg_debug; fprintf("csv Output Level: Auto-threshold range only\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "csvthonly")
            rna_spot_run.options.csv_output_level = 3;
            if arg_debug; fprintf("csv Output Level: Auto selected threshold only\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "maxzproj")
            rna_spot_run.options.use_max_proj = true;
            if arg_debug; fprintf("Use Max Z Projection: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "automaxth")
            rna_spot_run.options.t_max = 0;
            automaxth_flag = true;
            if arg_debug; fprintf("Automatically Determine Scan Max: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "autominth")
            rna_spot_run.options.t_min = 0;
            autominth_flag = true;
            if arg_debug; fprintf("Automatically Determine Scan Min: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "threepiece")
            if arg_debug; fprintf("Three Piece Fit: On\n"); end
            rna_spot_run.th_params.fit_strat = 'three_piece';
            rna_spot_run.th_params.reweight_fit = false;
            rna_spot_run.th_params.fit_to_log = true;
            rna_spot_run.th_params.madth_weight = 0.0;
            rna_spot_run.th_params.fit_weight = 1.0;
            rna_spot_run.th_params.fit_ri_weight = 0.0;
            lastkey = [];
        elseif strcmp(lastkey, "sensitive")
            if ~senspe_set
                if arg_debug; fprintf("Tuning Preset: Sensitivity\n"); end
                senspe_set = true;
                setThPreset(rna_spot_run, DEFAULT_PRESET_INDEX, 3);
                lastkey = [];
            end
        elseif strcmp(lastkey, "precise")
            if ~senspe_set
                if arg_debug; fprintf("Tuning Preset: Precise\n"); end
                senspe_set = true;
                setThPreset(rna_spot_run, DEFAULT_PRESET_INDEX, -3);
                lastkey = [];
            end
        elseif strcmp(lastkey, "noprobe")
            rna_spot_run.meta.noProbe_flag = true;
            if arg_debug; fprintf("Control/No Probe Flag: On\n"); end
            lastkey = [];
        end
        
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "imgname")
            rna_spot_run.img_name = argval;
            if arg_debug; fprintf("Image Name Set: %s\n", rna_spot_run.img_name); end
        elseif strcmp(lastkey, "tif")
            rna_spot_run.paths.img_path = argval;
            if arg_debug; fprintf("TIF Path Set: %s\n", rna_spot_run.paths.img_path); end
        elseif strcmp(lastkey, "matimg")
            rna_spot_run.paths.img_path = argval;
            if arg_debug; fprintf("MAT Path Set: %s\n", rna_spot_run.paths.img_path); end
        elseif strcmp(lastkey, "matvar")
            matvar = argval;
            if arg_debug; fprintf("MAT Variable Name: %s\n", matvar); end
        elseif strcmp(lastkey, "input")
            rna_spot_run.paths.img_path = argval;
            if arg_debug; fprintf("Input Path Set: %s\n", rna_spot_run.paths.img_path); end
        elseif strcmp(lastkey, "outdir")
            rna_spot_run.paths.out_dir = argval;
            if arg_debug; fprintf("Output Directory Set: %s\n", rna_spot_run.paths.out_dir); end
        elseif strcmp(lastkey, "outstem")
            [rna_spot_run.paths.out_dir, rna_spot_run.paths.out_namestem, ~] = fileparts(argval);
            if arg_debug; fprintf("Output Directory Set: %s\n", rna_spot_run.paths.out_dir); end
            if arg_debug; fprintf("Output Filename Stem Set: %s\n", rna_spot_run.paths.out_namestem); end
        elseif strcmp(lastkey, "cellseg")
            rna_spot_run.paths.cellseg_path = argval;
            if arg_debug; fprintf("CellSeg Data Path Set: %s\n", rna_spot_run.paths.cellseg_path); end
        elseif strcmp(lastkey, "csvout")
            rna_spot_run.paths.csv_out_path = argval;
            if arg_debug; fprintf("csv Output Path Set: %s\n", rna_spot_run.paths.csv_out_path); end
        elseif strcmp(lastkey, "runparamout")
            rna_spot_run.paths.params_out_path = argval;
            if arg_debug; fprintf("Run parameter dump path set: %s\n", rna_spot_run.paths.params_out_path); end
        elseif strcmp(lastkey, "chsamp")
            rna_spot_run.channels.rna_ch = round(Force2Num(argval));
            if arg_debug; fprintf("Sample Channel Set: %d\n", rna_spot_run.channels.rna_ch); end
        elseif strcmp(lastkey, "chtrans")
            rna_spot_run.channels.light_ch = round(Force2Num(argval));
            if arg_debug; fprintf("Light/TRANS Channel Set: %d\n", rna_spot_run.channels.light_ch); end
        elseif strcmp(lastkey, "chtotal")
            rna_spot_run.channels.total_ch = round(Force2Num(argval));
            if arg_debug; fprintf("Channel Count Set: %d\n", rna_spot_run.channels.total_ch); end
        elseif strcmp(lastkey, "ctrltif")
            rna_spot_run.paths.ctrl_img_path = argval;
            if arg_debug; fprintf("Control TIF Path Set: %s\n", rna_spot_run.paths.ctrl_img_path); end
        elseif strcmp(lastkey, "ctrlimg")
            rna_spot_run.paths.ctrl_img_path = argval;
            if arg_debug; fprintf("Control Image Path Set: %s\n", rna_spot_run.paths.ctrl_img_path); end
        elseif strcmp(lastkey, "chctrsamp")
            rna_spot_run.channels.ctrl_ch = round(Force2Num(argval));
            if arg_debug; fprintf("Control Channel Set: %d\n", rna_spot_run.channels.ctrl_ch); end
        elseif strcmp(lastkey, "chctrtotal")
            rna_spot_run.channels.ctrl_chcount = round(Force2Num(argval));
            if arg_debug; fprintf("Control Channel Count Set: %d\n", rna_spot_run.channels.ctrl_chcount); end
        elseif strcmp(lastkey, "thmin")
            if ~autominth_flag
                rna_spot_run.options.t_min = round(Force2Num(argval));
                if arg_debug; fprintf("Min Scan Threshold Set: %d\n", rna_spot_run.options.t_min); end
            end
        elseif strcmp(lastkey, "thmax")
            if ~automaxth_flag
                rna_spot_run.options.t_max = round(Force2Num(argval));
                if arg_debug; fprintf("Max Scan Threshold Set: %d\n", rna_spot_run.options.t_max); end
            end
        elseif strcmp(lastkey, "ztrim")
            rna_spot_run.dims.ztrim = round(Force2Num(argval));
            if arg_debug; fprintf("Z Trim Set: %d\n", rna_spot_run.dims.ztrim); end
        elseif strcmp(lastkey, "zmin")
            rna_spot_run.dims.z_min = round(Force2Num(argval));
            if arg_debug; fprintf("Z Min Set: %d\n", rna_spot_run.dims.z_min); end
        elseif strcmp(lastkey, "zmax")
            rna_spot_run.dims.z_max = round(Force2Num(argval));
            if arg_debug; fprintf("Z Max Set: %d\n", rna_spot_run.dims.z_max); end
        elseif strcmp(lastkey, "mfmin")
            rna_spot_run.th_params.mad_factor_min = Force2Num(argval);
            if arg_debug; fprintf("mad Factor Minimum Set: %.3f\n", rna_spot_run.th_params.mad_factor_min); end
        elseif strcmp(lastkey, "mfmax")
            rna_spot_run.th_params.mad_factor_max = Force2Num(argval);
            if arg_debug; fprintf("mad Factor Maximum Set: %.3f\n", rna_spot_run.th_params.mad_factor_max); end
        elseif strcmp(lastkey, "wszmin")
            rna_spot_run.options.winsize_min = round(Force2Num(argval));
            if arg_debug; fprintf("Window Size Minimum Set: %d\n", rna_spot_run.options.winsize_min); end
        elseif strcmp(lastkey, "wszmax")
            rna_spot_run.options.winsize_max = round(Force2Num(argval));
            if arg_debug; fprintf("Window Size Maximum Set: %d\n", rna_spot_run.options.winsize_max); end
        elseif strcmp(lastkey, "wszincr")
            rna_spot_run.options.winsize_incr = round(Force2Num(argval));
            if arg_debug; fprintf("Window Size Increment Set: %d\n", rna_spot_run.options.winsize_incr); end
        elseif strcmp(lastkey, "splitr")
            if rna_spot_run.th_params.spline_iterations > 0
                %Ignore if it was previously set to 0 by something like 'nospline'
                rna_spot_run.th_params.spline_iterations = round(Force2Num(argval));
                if arg_debug; fprintf("Spline Iterations: %d\n", rna_spot_run.th_params.spline_iterations); end
            end
        elseif strcmp(lastkey, "thmw")
            rna_spot_run.th_params.madth_weight = round(Force2Num(argval));
            if arg_debug; fprintf("Med/MAD Thresholding Weight Set: %f\n", rna_spot_run.th_params.madth_weight); end
        elseif strcmp(lastkey, "thfw")
            rna_spot_run.th_params.fit_weight = round(Force2Num(argval));
            if arg_debug; fprintf("Fit X Thresholding Weight Set: %f\n", rna_spot_run.th_params.fit_weight); end
        elseif strcmp(lastkey, "thiw")
            rna_spot_run.th_params.fit_ri_weight = round(Force2Num(argval));
            if arg_debug; fprintf("Fit Intersect Thresholding Weight Set: %f\n", rna_spot_run.th_params.fit_ri_weight); end
        elseif strcmp(lastkey, "fitlog")
            rna_spot_run.th_params.fit_to_log = Force2Bool(argval);
            if arg_debug; fprintf("Fit Piecewise to Log Plot: %d\n", rna_spot_run.th_params.fit_to_log); end
        elseif strcmp(lastkey, "stdfac")
            rna_spot_run.th_params.std_factor = Force2Num(argval);
            if arg_debug; fprintf("StDev Add Factor Set: %f\n", rna_spot_run.th_params.std_factor); end
        elseif strcmp(lastkey, "sensitivity")
            specval = Force2Num(argval);
            if specval > MAX_TH_PRESET_LEVEL; specval = MAX_TH_PRESET_LEVEL; end
            setThPreset(rna_spot_run, DEFAULT_PRESET_INDEX, specval);
            if arg_debug; fprintf("Sensitivity Preset Level Set: %d\n", specval); end
        elseif strcmp(lastkey, "precision")
            specval = Force2Num(argval);
            if specval > MAX_TH_PRESET_LEVEL; specval = MAX_TH_PRESET_LEVEL; end
            setThPreset(rna_spot_run, DEFAULT_PRESET_INDEX, specval * -1);
            if arg_debug; fprintf("Precision Preset Level Set: %d\n", specval); end
        elseif strcmp(lastkey, "voxelsize") | strcmp(lastkey, "pixelsize")
            rna_spot_run.meta.idims_voxel = parseDimsTo(argval, rna_spot_run.meta.idims_voxel);
            if rna_spot_run.meta.idims_voxel.z > 0
                if arg_debug; fprintf("Voxel Size: %d x %d x %d nm\n", ...
                        rna_spot_run.meta.idims_voxel.x, rna_spot_run.meta.idims_voxel.y, rna_spot_run.meta.idims_voxel.z); end
            else
                if arg_debug; fprintf("Pixel Size: %d x %d nm\n", rna_spot_run.meta.idims_voxel.x, rna_spot_run.meta.idims_voxel.y); end
            end
        elseif strcmp(lastkey, "expspotsize")
            rna_spot_run.meta.idims_expspot = parseDimsTo(argval, rna_spot_run.meta.idims_expspot);
            if rna_spot_run.meta.idims_expspot.z > 0
                if arg_debug; fprintf("Estimated Spot Size: %d x %d x %d nm\n", ...
                        rna_spot_run.meta.idims_expspot.x, rna_spot_run.meta.idims_expspot.y, rna_spot_run.meta.idims_expspot.z); end
            else
                if arg_debug; fprintf("Estimated Spot Size: %d x %d nm\n", rna_spot_run.meta.idims_expspot.x, rna_spot_run.meta.idims_expspot.y); end
            end
        elseif strcmp(lastkey, "probetype")
            rna_spot_run.meta.type_probe = argval;
            if arg_debug; fprintf("Probe Set: %s\n", rna_spot_run.meta.type_probe); end
        elseif strcmp(lastkey, "target")
            rna_spot_run.meta.type_target = argval;
            if arg_debug; fprintf("Target Set: %s\n", rna_spot_run.meta.type_target); end
        elseif strcmp(lastkey, "targettype")
            rna_spot_run.meta.type_targetmol = argval;
            if arg_debug; fprintf("Target Type Set: %s\n", rna_spot_run.meta.type_targetmol); end
        elseif strcmp(lastkey, "species")
            rna_spot_run.meta.type_species = argval;
            if arg_debug; fprintf("Species Set: %s\n", rna_spot_run.meta.type_species); end
        elseif strcmp(lastkey, "celltype")
            rna_spot_run.meta.type_cell = argval;
            if arg_debug; fprintf("Cell Type Set: %s\n", rna_spot_run.meta.type_cell); end
        elseif strcmp(lastkey, "threads")
            thread_request = Force2Num(argval);
            if thread_request < 1; thread_request = 1; end
            if arg_debug; fprintf("Requested Threads: %d\n", thread_request); end
        elseif strcmp(lastkey, "gaussrad")
            rna_spot_run.options.dtune_gaussrad = Force2Num(argval);
            if arg_debug; fprintf("XY Gaussian Radius Set: %s\n", rna_spot_run.options.dtune_gaussrad); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

% ========================== Arg Check ==========================

if strcmp(rna_spot_run.th_params.fit_strat, 'three_piece')
    rna_spot_run.th_params.reweight_fit = false;
end

% ========================== Run ==========================

% rna_spot_run = Adapter_RNASpots(img_name, tif_path, rna_ch, light_ch, total_ch,...
%     out_dir, t_min, t_max, ztrim, cellseg_path, ctrl_path, ctrl_ch, ctrl_chcount, ttune_winsize,...
%     ttune_madfactor, overwrite_output, false);

if ~isempty(rna_spot_run.paths.img_path) & endsWith(rna_spot_run.paths.img_path, '.mat')
    rna_spot_run.channels.total_ch = 1;
    rna_spot_run.channels.rna_ch = 1;
    if ~isempty(matvar)
        img_ch_set = MatImages.loadImageChannels(rna_spot_run.paths.img_path, matvar);
    else
        img_ch_set = MatImages.loadImageChannels(rna_spot_run.paths.img_path);
    end
    rna_spot_run = Adapter_RNASpots(rna_spot_run, false, img_ch_set, debug_lvl, thread_request);
else
    rna_spot_run = Adapter_RNASpots(rna_spot_run, false, [], debug_lvl, thread_request);
end

end

function idims = parseDimsTo(dims_arg, trg_struct)
    idims = trg_struct;
	if ischar(dims_arg)
        %Read as string
        %'(x,y,z)'
        repl_str = replace(dims_arg, {'(', ')'}, '');
        split_str = split(repl_str, ',');
        ndims = size(split_str,1);
        idims.x = Force2Num(split_str{1,1});
        if ndims > 1; idims.y = Force2Num(split_str{2,1}); end
        if ndims > 2; idims.z = Force2Num(split_str{3,1}); end
	else
        %Try to read as vector
        if isvector(dims_arg)
            ndims = size(dims_arg,2);
            idims.x = dims_arg(1);
            if ndims > 1; idims.y = dims_arg(2); end
            if ndims > 2; idims.z = dims_arg(3); end
        end
	end
end

function spotsrun = setThPreset(spotsrun, middle, modifier)
    spotsrun = RNAThreshold.applyPreset(spotsrun, middle + modifier);
end
