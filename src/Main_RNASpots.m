%%
%%

function rna_spot_run = Main_RNASpots(varargin)

% ========================== Process args ==========================
rna_spot_run = RNASpotsRun.initDefault();

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
            fprintf("Overwrite Output: On\n");
            lastkey = [];
        elseif strcmp(lastkey, "usespc")
            rna_spot_run.ttune_use_rawcurve = true;
            fprintf("Use Raw Spot Count Curve: On\n");
            lastkey = [];
        elseif strcmp(lastkey, "usedfc")
            rna_spot_run.ttune_use_diffcurve = true;
            fprintf("Use Diff Curve: On\n");
            lastkey = [];
        elseif strcmp(lastkey, "nospline")
            rna_spot_run.ttune_spline_itr = 0;
            fprintf("Use Spline: Off\n");
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
            fprintf("Image Name Set: %s\n", rna_spot_run.img_name);
        elseif strcmp(lastkey, "tif")
            rna_spot_run.tif_path = argval;
            fprintf("TIF Path Set: %s\n", rna_spot_run.tif_path);
        elseif strcmp(lastkey, "outdir")
            rna_spot_run.out_dir = argval;
            fprintf("Output Directory Set: %s\n", rna_spot_run.out_dir);
        elseif strcmp(lastkey, "cellseg")
            rna_spot_run.cellseg_path = argval;
            fprintf("CellSeg Data Path Set: %s\n", rna_spot_run.cellseg_path);
        elseif strcmp(lastkey, "chsamp")
            rna_spot_run.rna_ch = round(Force2Num(argval));
            fprintf("Sample Channel Set: %d\n", rna_spot_run.rna_ch);
        elseif strcmp(lastkey, "chtrans")
            rna_spot_run.light_ch = round(Force2Num(argval));
            fprintf("Light/TRANS Channel Set: %d\n", rna_spot_run.light_ch);
        elseif strcmp(lastkey, "chtotal")
            rna_spot_run.total_ch = round(Force2Num(argval));
            fprintf("Channel Count Set: %d\n", rna_spot_run.total_ch);
        elseif strcmp(lastkey, "ctrltif")
            rna_spot_run.ctrl_path = argval;
            fprintf("Control TIF Path Set: %s\n", rna_spot_run.ctrl_path);
        elseif strcmp(lastkey, "chctrsamp")
            rna_spot_run.ctrl_ch = round(Force2Num(argval));
            fprintf("Control Channel Set: %d\n", rna_spot_run.ctrl_ch);
        elseif strcmp(lastkey, "chctrtotal")
            rna_spot_run.ctrl_chcount = round(Force2Num(argval));
            fprintf("Control Channel Count Set: %d\n", rna_spot_run.ctrl_chcount);
        elseif strcmp(lastkey, "thmin")
            rna_spot_run.t_min = round(Force2Num(argval));
            fprintf("Min Scan Threshold Set: %d\n", rna_spot_run.t_min);
        elseif strcmp(lastkey, "thmax")
            rna_spot_run.t_max = round(Force2Num(argval));
            fprintf("Max Scan Threshold Set: %d\n", rna_spot_run.t_max);
        elseif strcmp(lastkey, "ztrim")
            rna_spot_run.ztrim = round(Force2Num(argval));
            fprintf("Z Trim Set: %d\n", rna_spot_run.ztrim);
        elseif strcmp(lastkey, "zmin")
            rna_spot_run.z_min = round(Force2Num(argval));
            fprintf("Z Min Set: %d\n", rna_spot_run.z_min);
        elseif strcmp(lastkey, "zmax")
            rna_spot_run.z_max = round(Force2Num(argval));
            fprintf("Z Max Set: %d\n", rna_spot_run.z_max);
        elseif strcmp(lastkey, "mfmin")
            rna_spot_run.ttune_madf_min = Force2Num(argval);
            fprintf("mad Factor Minimum Set: %.3f\n", rna_spot_run.ttune_madf_min);
        elseif strcmp(lastkey, "mfmax")
            rna_spot_run.ttune_madf_max = Force2Num(argval);
            fprintf("mad Factor Maximum Set: %.3f\n", rna_spot_run.ttune_madf_max);
        elseif strcmp(lastkey, "wszmin")
            rna_spot_run.ttune_winsz_min = round(Force2Num(argval));
            fprintf("Window Size Minimum Set: %d\n", rna_spot_run.ttune_winsz_min);
        elseif strcmp(lastkey, "wszmax")
            rna_spot_run.ttune_winsz_max = round(Force2Num(argval));
            fprintf("Window Size Maximum Set: %d\n", rna_spot_run.ttune_winsz_max);
        elseif strcmp(lastkey, "wszincr")
            rna_spot_run.ttune_winsz_incr = round(Force2Num(argval));
            fprintf("Window Size Increment Set: %d\n", rna_spot_run.ttune_winsz_incr);
        elseif strcmp(lastkey, "splitr")
            if rna_spot_run.ttune_spline_itr > 0
                %Ignore if it was previously set to 0 by something like 'nospline'
                rna_spot_run.ttune_spline_itr = round(Force2Num(argval));
                fprintf("Spline Iterations: %d\n", rna_spot_run.ttune_spline_itr);
            end
        elseif strcmp(lastkey, "probetype")
            rna_spot_run.type_probe = argval;
            fprintf("Probe Set: %s\n", rna_spot_run.type_probe);
        elseif strcmp(lastkey, "target")
            rna_spot_run.type_target = argval;
            fprintf("Target Set: %s\n", rna_spot_run.type_target);
        elseif strcmp(lastkey, "targettype")
            rna_spot_run.type_targetmol = argval;
            fprintf("Target Type Set: %s\n", rna_spot_run.type_targetmol);
        elseif strcmp(lastkey, "species")
            rna_spot_run.type_species = argval;
            fprintf("Species Set: %s\n", rna_spot_run.type_species);
        elseif strcmp(lastkey, "celltype")
            rna_spot_run.type_cell = argval;
            fprintf("Cell Type Set: %s\n", rna_spot_run.type_cell);
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

% ========================== Run ==========================

% rna_spot_run = Adapter_RNASpots(img_name, tif_path, rna_ch, light_ch, total_ch,...
%     out_dir, t_min, t_max, ztrim, cellseg_path, ctrl_path, ctrl_ch, ctrl_chcount, ttune_winsize,...
%     ttune_madfactor, overwrite_output, false);

rna_spot_run = Adapter_RNASpots(rna_spot_run, false);

end