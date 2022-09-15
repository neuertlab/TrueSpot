%%
%%

function rna_spot_run = Main_RNASpots(varargin)

% ========================== Process args ==========================
rna_spot_run = RNASpotsRun.initDefault();

%Highest level set to true is what is used, regardless of other flags.
arg_debug = true; %CONSTANT used for debugging arg parser.
debug_lvl = 0;
senspe_set = false;

lastkey = [];
for i = 1:nargin
    argval = varargin{i};
    if startsWith(argval, "-")
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
                rna_spot_run.ttune_fit_strat = 0;
                rna_spot_run.ttune_reweight_fit = false;
                rna_spot_run.ttune_fit_to_log = true;
                rna_spot_run.ttune_thweight_med = 0.0;
                rna_spot_run.ttune_thweight_fit = 1.0;
                rna_spot_run.ttune_thweight_fisect = 0.0;
                lastkey = [];
            end
        elseif strcmp(lastkey, "specific")
            if ~senspe_set
                if arg_debug; fprintf("Tuning Preset: Specificity\n"); end
                senspe_set = true;
                rna_spot_run.ttune_fit_strat = 3;
                rna_spot_run.ttune_reweight_fit = false;
                rna_spot_run.ttune_fit_to_log = true;
                rna_spot_run.ttune_thweight_med = 0.0;
                rna_spot_run.ttune_thweight_fit = 1.0;
                rna_spot_run.ttune_thweight_fisect = 0.0;
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

rna_spot_run = Adapter_RNASpots(rna_spot_run, false, debug_lvl);

end