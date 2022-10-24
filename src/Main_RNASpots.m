%%
%%

function rna_spot_run = Main_RNASpots(varargin)
addpath('./core');

% ========================== Process args ==========================
rna_spot_run = RNASpotsRun.initDefault();
rna_spot_run = setDefaultParams(rna_spot_run);

%Highest level set to true is what is used, regardless of other flags.
arg_debug = true; %CONSTANT used for debugging arg parser.
debug_lvl = 0;
senspe_set = false;
matvar = [];

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
            rna_spot_run.overwrite_output = true;
            if arg_debug; fprintf("Overwrite Output: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "usespc")
            rna_spot_run.ttune_use_rawcurve = true;
            if arg_debug; fprintf("Use Raw Spot Count Curve: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "usedfc")
            rna_spot_run.ttune_use_diffcurve = true;
            if arg_debug; fprintf("Use Diff Curve: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "nospline")
            rna_spot_run.ttune_spline_itr = 0;
            if arg_debug; fprintf("Use Spline: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "fitsegrsq")
            rna_spot_run.ttune_reweight_fit = false;
            if arg_debug; fprintf("Use weighted segmented fit: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "fitwavg")
            rna_spot_run.ttune_reweight_fit = true;
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
        elseif strcmp(lastkey, "threepiece")
            if arg_debug; fprintf("Three Piece Fit: On\n"); end
            rna_spot_run.ttune_fit_strat = 3;
            rna_spot_run.ttune_reweight_fit = false;
            rna_spot_run.ttune_fit_to_log = true;
            rna_spot_run.ttune_thweight_med = 0.0;
            rna_spot_run.ttune_thweight_fit = 1.0;
            rna_spot_run.ttune_thweight_fisect = 0.0;
            lastkey = [];
        elseif strcmp(lastkey, "sensitive")
            if ~senspe_set
                if arg_debug; fprintf("Tuning Preset: Sensitivity\n"); end
                senspe_set = true;
                rna_spot_run = setSensitive1(rna_spot_run);
                lastkey = [];
            end
        elseif strcmp(lastkey, "specific")
            if ~senspe_set
                if arg_debug; fprintf("Tuning Preset: Specificity\n"); end
                senspe_set = true;
                rna_spot_run = setSpecific1(rna_spot_run);
                lastkey = [];
            end
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
            rna_spot_run.tif_path = argval;
            if arg_debug; fprintf("TIF Path Set: %s\n", rna_spot_run.tif_path); end
        elseif strcmp(lastkey, "matimg")
            rna_spot_run.sample_matpath = argval;
            if arg_debug; fprintf("MAT Path Set: %s\n", rna_spot_run.sample_matpath); end
        elseif strcmp(lastkey, "matvar")
            matvar = argval;
            if arg_debug; fprintf("MAT Variable Name: %s\n", matvar); end
        elseif strcmp(lastkey, "outdir")
            rna_spot_run.out_dir = argval;
            if arg_debug; fprintf("Output Directory Set: %s\n", rna_spot_run.out_dir); end
        elseif strcmp(lastkey, "cellseg")
            rna_spot_run.cellseg_path = argval;
            if arg_debug; fprintf("CellSeg Data Path Set: %s\n", rna_spot_run.cellseg_path); end
        elseif strcmp(lastkey, "chsamp")
            rna_spot_run.rna_ch = round(Force2Num(argval));
            if arg_debug; fprintf("Sample Channel Set: %d\n", rna_spot_run.rna_ch); end
        elseif strcmp(lastkey, "chtrans")
            rna_spot_run.light_ch = round(Force2Num(argval));
            if arg_debug; fprintf("Light/TRANS Channel Set: %d\n", rna_spot_run.light_ch); end
        elseif strcmp(lastkey, "chtotal")
            rna_spot_run.total_ch = round(Force2Num(argval));
            if arg_debug; fprintf("Channel Count Set: %d\n", rna_spot_run.total_ch); end
        elseif strcmp(lastkey, "ctrltif")
            rna_spot_run.ctrl_path = argval;
            if arg_debug; fprintf("Control TIF Path Set: %s\n", rna_spot_run.ctrl_path); end
        elseif strcmp(lastkey, "chctrsamp")
            rna_spot_run.ctrl_ch = round(Force2Num(argval));
            if arg_debug; fprintf("Control Channel Set: %d\n", rna_spot_run.ctrl_ch); end
        elseif strcmp(lastkey, "chctrtotal")
            rna_spot_run.ctrl_chcount = round(Force2Num(argval));
            if arg_debug; fprintf("Control Channel Count Set: %d\n", rna_spot_run.ctrl_chcount); end
        elseif strcmp(lastkey, "thmin")
            rna_spot_run.t_min = round(Force2Num(argval));
            if arg_debug; fprintf("Min Scan Threshold Set: %d\n", rna_spot_run.t_min); end
        elseif strcmp(lastkey, "thmax")
            rna_spot_run.t_max = round(Force2Num(argval));
            if arg_debug; fprintf("Max Scan Threshold Set: %d\n", rna_spot_run.t_max); end
        elseif strcmp(lastkey, "ztrim")
            rna_spot_run.ztrim = round(Force2Num(argval));
            if arg_debug; fprintf("Z Trim Set: %d\n", rna_spot_run.ztrim); end
        elseif strcmp(lastkey, "zmin")
            rna_spot_run.z_min = round(Force2Num(argval));
            if arg_debug; fprintf("Z Min Set: %d\n", rna_spot_run.z_min); end
        elseif strcmp(lastkey, "zmax")
            rna_spot_run.z_max = round(Force2Num(argval));
            if arg_debug; fprintf("Z Max Set: %d\n", rna_spot_run.z_max); end
        elseif strcmp(lastkey, "mfmin")
            rna_spot_run.ttune_madf_min = Force2Num(argval);
            if arg_debug; fprintf("mad Factor Minimum Set: %.3f\n", rna_spot_run.ttune_madf_min); end
        elseif strcmp(lastkey, "mfmax")
            rna_spot_run.ttune_madf_max = Force2Num(argval);
            if arg_debug; fprintf("mad Factor Maximum Set: %.3f\n", rna_spot_run.ttune_madf_max); end
        elseif strcmp(lastkey, "wszmin")
            rna_spot_run.ttune_winsz_min = round(Force2Num(argval));
            if arg_debug; fprintf("Window Size Minimum Set: %d\n", rna_spot_run.ttune_winsz_min); end
        elseif strcmp(lastkey, "wszmax")
            rna_spot_run.ttune_winsz_max = round(Force2Num(argval));
            if arg_debug; fprintf("Window Size Maximum Set: %d\n", rna_spot_run.ttune_winsz_max); end
        elseif strcmp(lastkey, "wszincr")
            rna_spot_run.ttune_winsz_incr = round(Force2Num(argval));
            if arg_debug; fprintf("Window Size Increment Set: %d\n", rna_spot_run.ttune_winsz_incr); end
        elseif strcmp(lastkey, "splitr")
            if rna_spot_run.ttune_spline_itr > 0
                %Ignore if it was previously set to 0 by something like 'nospline'
                rna_spot_run.ttune_spline_itr = round(Force2Num(argval));
                if arg_debug; fprintf("Spline Iterations: %d\n", rna_spot_run.ttune_spline_itr); end
            end
        elseif strcmp(lastkey, "thmw")
            rna_spot_run.ttune_thweight_med = round(Force2Num(argval));
            if arg_debug; fprintf("Med/MAD Thresholding Weight Set: %f\n", rna_spot_run.ttune_thweight_med); end
        elseif strcmp(lastkey, "thfw")
            rna_spot_run.ttune_thweight_fit = round(Force2Num(argval));
            if arg_debug; fprintf("Fit X Thresholding Weight Set: %f\n", rna_spot_run.ttune_thweight_fit); end
        elseif strcmp(lastkey, "thiw")
            rna_spot_run.ttune_thweight_fisect = round(Force2Num(argval));
            if arg_debug; fprintf("Fit Intersect Thresholding Weight Set: %f\n", rna_spot_run.ttune_thweight_fisect); end
        elseif strcmp(lastkey, "fitlog")
            rna_spot_run.ttune_fit_to_log = Force2Bool(argval);
            if arg_debug; fprintf("Fit Piecewise to Log Plot: %d\n", rna_spot_run.ttune_fit_to_log); end
        elseif strcmp(lastkey, "stdfac")
            rna_spot_run.ttune_std_factor = Force2Num(argval);
            if arg_debug; fprintf("StDev Add Factor Set: %f\n", rna_spot_run.ttune_std_factor); end
        elseif strcmp(lastkey, "sensitivity")
            specval = Force2Num(argval);
            if specval == 0
                rna_spot_run = setDefaultParams(rna_spot_run);
            elseif specval == 1
                rna_spot_run = setSensitive1(rna_spot_run);
            elseif specval >= 2
                rna_spot_run = setSensitive2(rna_spot_run);
            else
                %Treated as 0
                specval = 0;
                rna_spot_run = setDefaultParams(rna_spot_run);
            end
            if arg_debug; fprintf("Sensitivity Preset Level Set: %d\n", specval); end
        elseif strcmp(lastkey, "specificity")
            specval = Force2Num(argval);
            if specval == 0
                rna_spot_run = setDefaultParams(rna_spot_run);
            elseif specval == 1
                rna_spot_run = setSpecific1(rna_spot_run);
            elseif specval >= 2
                rna_spot_run = setSpecific2(rna_spot_run);
            else
                %Treated as 0
                specval = 0;
                rna_spot_run = setDefaultParams(rna_spot_run);
            end
            if arg_debug; fprintf("Specificity Preset Level Set: %d\n", specval); end
        elseif strcmp(lastkey, "voxelsize") | strcmp(lastkey, "pixelsize")
            rna_spot_run.idims_voxel = parseDimsTo(argval, rna_spot_run.idims_voxel);
            if rna_spot_run.idims_voxel.z > 0
                if arg_debug; fprintf("Voxel Size: %d x %d x %d nm\n", ...
                        rna_spot_run.idims_voxel.x, rna_spot_run.idims_voxel.y, rna_spot_run.idims_voxel.z); end
            else
                if arg_debug; fprintf("Pixel Size: %d x %d nm\n", rna_spot_run.idims_voxel.x, rna_spot_run.idims_voxel.y); end
            end
        elseif strcmp(lastkey, "expspotsize")
            rna_spot_run.idims_expspot = parseDimsTo(argval, rna_spot_run.idims_expspot);
            if rna_spot_run.idims_expspot.z > 0
                if arg_debug; fprintf("Estimated Spot Size: %d x %d x %d nm\n", ...
                        rna_spot_run.idims_expspot.x, rna_spot_run.idims_expspot.y, rna_spot_run.idims_expspot.z); end
            else
                if arg_debug; fprintf("Estimated Spot Size: %d x %d nm\n", rna_spot_run.idims_expspot.x, rna_spot_run.idims_expspot.y); end
            end
        elseif strcmp(lastkey, "probetype")
            rna_spot_run.type_probe = argval;
            if arg_debug; fprintf("Probe Set: %s\n", rna_spot_run.type_probe); end
        elseif strcmp(lastkey, "target")
            rna_spot_run.type_target = argval;
            if arg_debug; fprintf("Target Set: %s\n", rna_spot_run.type_target); end
        elseif strcmp(lastkey, "targettype")
            rna_spot_run.type_targetmol = argval;
            if arg_debug; fprintf("Target Type Set: %s\n", rna_spot_run.type_targetmol); end
        elseif strcmp(lastkey, "species")
            rna_spot_run.type_species = argval;
            if arg_debug; fprintf("Species Set: %s\n", rna_spot_run.type_species); end
        elseif strcmp(lastkey, "celltype")
            rna_spot_run.type_cell = argval;
            if arg_debug; fprintf("Cell Type Set: %s\n", rna_spot_run.type_cell); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

% ========================== Arg Check ==========================

if rna_spot_run.ttune_fit_strat == 3
    rna_spot_run.ttune_reweight_fit = false;
end

% ========================== Run ==========================

% rna_spot_run = Adapter_RNASpots(img_name, tif_path, rna_ch, light_ch, total_ch,...
%     out_dir, t_min, t_max, ztrim, cellseg_path, ctrl_path, ctrl_ch, ctrl_chcount, ttune_winsize,...
%     ttune_madfactor, overwrite_output, false);

if ~isempty(rna_spot_run.sample_matpath)
    rna_spot_run.total_ch = 1;
    rna_spot_run.rna_ch = 1;
    if ~isempty(matvar)
        img_ch_set = MatImages.loadImageChannels(rna_spot_run.sample_matpath, matvar);
    else
        img_ch_set = MatImages.loadImageChannels(rna_spot_run.sample_matpath);
    end
    rna_spot_run = Adapter_RNASpots(rna_spot_run, false, img_ch_set, debug_lvl);
else
    rna_spot_run = Adapter_RNASpots(rna_spot_run, false, [], debug_lvl);
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

function spotsrun = setDefaultParams(spotsrun)
    spotsrun.ttune_winsz_min = 3;
    spotsrun.ttune_fit_strat = 0;
    spotsrun.ttune_reweight_fit = false;
    spotsrun.ttune_fit_to_log = true;
    spotsrun.ttune_thweight_med = 0.25;
    spotsrun.ttune_thweight_fit = 0.0;
    spotsrun.ttune_thweight_fisect = 0.75;
    spotsrun.ttune_std_factor = 1.0;
end

function spotsrun = setSensitive1(spotsrun)
    spotsrun.ttune_winsz_min = 3;
    spotsrun.ttune_fit_strat = 0;
    spotsrun.ttune_reweight_fit = false;
    spotsrun.ttune_fit_to_log = true;
    spotsrun.ttune_thweight_med = 0.0;
    spotsrun.ttune_thweight_fit = 0.2;
    spotsrun.ttune_thweight_fisect = 0.8;
    spotsrun.ttune_std_factor = 0.0;
end

function spotsrun = setSensitive2(spotsrun)
    spotsrun.ttune_winsz_min = 3;
    spotsrun.ttune_fit_strat = 0;
    spotsrun.ttune_reweight_fit = false;
    spotsrun.ttune_fit_to_log = true;
    spotsrun.ttune_thweight_med = 0.0;
    spotsrun.ttune_thweight_fit = 1.0;
    spotsrun.ttune_thweight_fisect = 0.0;
    spotsrun.ttune_std_factor = 0.0;
end

function spotsrun = setSpecific1(spotsrun)
    spotsrun.ttune_winsz_min = 6;
    spotsrun.ttune_fit_strat = 0;
    spotsrun.ttune_reweight_fit = false;
    spotsrun.ttune_fit_to_log = true;
    spotsrun.ttune_thweight_med = 0.5;
    spotsrun.ttune_thweight_fit = 0.0;
    spotsrun.ttune_thweight_fisect = 0.5;
    spotsrun.ttune_std_factor = 1.0;
end

function spotsrun = setSpecific2(spotsrun)
    spotsrun.ttune_winsz_min = 6;
    spotsrun.ttune_fit_strat = 0;
    spotsrun.ttune_reweight_fit = false;
    spotsrun.ttune_fit_to_log = true;
    spotsrun.ttune_thweight_med = 0.5;
    spotsrun.ttune_thweight_fit = 0.0;
    spotsrun.ttune_thweight_fisect = 0.5;
    spotsrun.ttune_std_factor = 2.0;
end