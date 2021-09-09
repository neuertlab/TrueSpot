%%
%%

function rna_spot_run = Main_RNASpots(img_name, tif_path, rna_ch, light_ch, total_ch,...
    out_dir, t_min, t_max, ztrim, cellseg_path, ctrl_path, ctrl_ch, ctrl_chcount, ttune_winsize, ttune_wscorethresh,...
    overwrite_output)

addpath('./core');

%Build a Run object for easy save/load
spotsrun = RNASpotsRun;
rna_spot_run = spotsrun;

%Save string args to run obj
spotsrun.img_name = img_name;
spotsrun.tif_path = tif_path;
spotsrun.out_dir = out_dir;
spotsrun.cellseg_path = cellseg_path;
spotsrun.ctrl_path = ctrl_path;

%Check arguments
%Force to numbers
spotsrun.rna_ch = Force2Num(rna_ch);
spotsrun.light_ch = Force2Num(light_ch);
spotsrun.total_ch = Force2Num(total_ch);
spotsrun.ctrl_ch = Force2Num(ctrl_ch);
spotsrun.ctrl_chcount = Force2Num(ctrl_chcount);
spotsrun.t_min = Force2Num(t_min);
spotsrun.t_max = Force2Num(t_max);
spotsrun.ztrim = Force2Num(ztrim);
spotsrun.ttune_winsize = Force2Num(ttune_winsize);
spotsrun.ttune_wscorethresh = Force2Num(ttune_wscorethresh);
spotsrun.overwrite_output = Force2Bool(overwrite_output);

%Check required arguments
if (isempty(tif_path)) || (~ischar(tif_path))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Input image path is required."), true);
    return;
end
if (~isempty(out_dir)) && (~ischar(out_dir))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Output directory argument is invalid."), true);
    return;
end
if (~isempty(cellseg_path)) && (~ischar(cellseg_path))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Cellseg data path argument is invalid."), true);
    return;
end
if total_ch < 1
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Input image must have at least one channel."), true);
    return;
end
if light_ch > total_ch
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Light (TRANS) channel index is invalid."), true);
    return;
end
if rna_ch > total_ch
    RNA_Fisher_State.outputMessageLineStatic(sprintf("RNA/Signal channel index is invalid."), true);
    return;
end

%Set defaults
if ttune_winsize < 1
    spotsrun.ttune_winsize = 10;
end
if ttune_wscorethresh < 0
    spotsrun.ttune_wscorethresh = 0.9;
end
if t_min < 1
    spotsrun.t_min = 1;
end
if t_max < 1
    spotsrun.t_max = 300;
end
if isempty(spotsrun.out_dir)
    %Defaults to input directory
    [spotsrun.out_dir, ~, ~] = fileparts(spotsrun.tif_path);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Output path not provided. Set to %s", spotsrun.out_dir), true);
end
if isempty(spotsrun.img_name)
    %Defaults to input file name
    [~, spotsrun.img_name, ~] = fileparts(spotsrun.tif_path);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Image name not provided. Set to %s", spotsrun.img_name), true);
end
if ~isempty(spotsrun.ctrl_path)
    if ~ischar(spotsrun.ctrl_path)
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Control path is not valid. Setting to empty."), true);
        spotsrun.ctrl_path = "";
    end
end
if spotsrun.ctrl_ch > spotsrun.ctrl_chcount
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control RNA/signal channel index is invalid. Control will not be used."), true);
    spotsrun.ctrl_path = "";
end

%Debug print
RNA_Fisher_State.outputMessageLineStatic(sprintf("Running RNASpots with the following parameters..."), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("img_name = %s", spotsrun.img_name), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("tif_path = %s", spotsrun.tif_path), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("rna_ch = %d", spotsrun.rna_ch), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("light_ch = %d", spotsrun.light_ch), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("total_ch = %d", spotsrun.total_ch), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("out_dir = %s", spotsrun.out_dir), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("t_min = %d", spotsrun.t_min), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("t_max = %d", spotsrun.t_max), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("ztrim = %d", spotsrun.ztrim), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("cellseg_path = %s", spotsrun.cellseg_path), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("ctrl_path = %s", spotsrun.ctrl_path), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("ctrl_ch = %d", spotsrun.ctrl_ch), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("ctrl_chcount = %d", spotsrun.ctrl_chcount), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("ttune_winsize = %d", spotsrun.ttune_winsize), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("ttune_wscorethresh = %f", spotsrun.ttune_wscorethresh), false);
RNA_Fisher_State.outputMessageLineStatic(sprintf("overwrite_output = %d", spotsrun.overwrite_output), false);

%Do background extraction if arguments provided.
%!! Don't redo if target exists and overwrite output is false!
spotsrun.bkg_path = [spotsrun.out_dir filesep 'bkgmask' filesep spotsrun.img_name '_bkg'];
if isempty(spotsrun.ctrl_path)
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control not provided. Will attempt to use background."), true);
    if ~spotsrun.overwrite_output && isfile(spotsrun.bkg_path)
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Background mask already exists at %s! Skipping extraction...", spotsrun.bkg_path), true);
    else
        if (~isempty(spotsrun.cellseg_path)) && (spotsrun.light_ch > 0)
            %fprintf("Extracting background...\n");
            %Make sure cellseg data file exists
            if (isfile(spotsrun.cellseg_path))
                Main_BackgroundMask(spotsrun.tif_path, spotsrun.cellseg_path, spotsrun.bkg_path, spotsrun.total_ch, spotsrun.light_ch, true);
            else
                RNA_Fisher_State.outputMessageLineStatic(sprintf("Cellseg file %s does not exist! Skipping background extraction...", spotsrun.cellseg_path), true);
            end
        else
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Control path was empty or light channel index is invalid. Skipping background extraction..."), true);
        end
    end
else
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control image provided - background extraction skipped."), true);
end

%Run spot detect
%!! Don't redo if target exists and overwrite output is false!
RNA_Fisher_State.outputMessageLineStatic(sprintf("Running spot detect... (This may take a few hours on large files)"), true);
spotsrun.out_stem = Main_RNASpotDetect(spotsrun.img_name, spotsrun.tif_path, spotsrun.out_dir,...
    spotsrun.rna_ch, spotsrun.total_ch, spotsrun.t_min, spotsrun.t_max, true, spotsrun.overwrite_output);

%Run spot detect on control (if a tif path was provided)
%   Detect if input control path is a tif. If not, assume the input is a spot
%       detect results path stem.
%!! Don't redo if target exists and overwrite output is false!
if ~isempty(spotsrun.ctrl_path)
    if endsWith(spotsrun.ctrl_path, ".tif")
        if ~isfile(spotsrun.ctrl_path)
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Control path %s does not exist. Aborting...", spotsrun.ctrl_path), true);
            return;
        end
        
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Running spot detect on control image... (This may take a few hours on large files)"), true);
        spotsrun.ctrl_stem = Main_RNASpotDetect([spotsrun.img_name '_Control'], spotsrun.ctrl_path, spotsrun.out_dir,...
            spotsrun.ctrl_ch, spotsrun.ctrl_chcount, spotsrun.t_min, spotsrun.t_max, true, spotsrun.overwrite_output);
    else
        spotsrun.ctrl_stem = spotsrun.ctrl_path;
    end
end

%Filter spots for bkg (if mask exists)
spotsrun.bkg_filter_stem = [spotsrun.out_stem '_bkgmasked'];
if isfile(spotsrun.bkg_path)
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Filtering detected coordinates through background mask for control..."), true);
    
    spotsrun.updateBackgroundFilteredCoords();
    
    %Set as control if there isn't one.
    if isempty(spotsrun.ctrl_stem)
        spotsrun.ctrl_stem = spotsrun.bkg_filter_stem;
    end
end

%Get image dimensions (for ztrim)
spotsrun.updateImageDimensions();
idims = spotsrun.idims_sample;
z_min_trimmed = 1;
z_max_trimmed = idims.z;

%Load spot count tables and apply ztrim
if ztrim > 0
    %Create ztrim mask.
    ztrim_mask = true(idims.y, idims.x, idims.z);
    z_min_trimmed = ztrim;
    ztrim_mask(:,:,1:z_min_trimmed) = false;
    z_max_trimmed = idims.z - ztrim +  1;
    ztrim_mask(:,:,z_max_trimmed:idims.z) = false;
    
    %Load the tables
    [~, coords_sample] = spotsrun.loadCoordinateTable();
    
    %Mask the tables
    [spots_sample, ~] = RNA_Threshold_Common.mask_spots(ztrim_mask, coords_sample);
    
    %--- And repeat with the control
    if ~isfile(bkg_path)
        %Need to retrieve dims again and make a new mask.
        idims = spotsrun.idims_ctrl;
        ztrim_mask = true(idims.y, idims.x, idims.z);
        zbot = ztrim;
        ztrim_mask(:,:,1:zbot) = false;
        ztop = idims.z - ztrim +  1;
        ztrim_mask(:,:,ztop:idims.z) = false;
    end
    [~, coords_control] = spotsrun.loadControlCoordinateTable();
    [spots_control, ~] = RNA_Threshold_Common.mask_spots(ztrim_mask, coords_control);
else
    %Just load the coord and spot tables
    [~, spots_sample] = spotsrun.loadSpotsTable();
    [~, spots_control] = spotsrun.loadControlSpotsTable();
end

%Detect threshold
[thresh, win_stdevs] = RNA_Threshold_Common.estimateThreshold(spots_sample, spots_control, spotsrun.ttune_winsize, 0.0, spotsrun.ttune_wscorethresh);
RNA_Fisher_State.outputMessageLineStatic(sprintf("Auto threshold selected: %d", thresh), true);
spotsrun.intensity_threshold = thresh;

%Print plots & image representation
    %Spot plots (log and linear scale) - w/ chosen threshold marked
    %Max projection of sample and control w circled spots (also max proj)
    
plots_dir = [spotsrun.out_stem filesep 'plots'];
RNA_Fisher_State.outputMessageLineStatic(sprintf("Now generating plots..."), true);
probeNames = [{spotsrun.img_name};
              {'Control'}];
RNA_Threshold_Plotter.plotPreprocessedData(spotsrun.out_stem, [{spotsrun.ctrl_stem}], probeNames,...
                thresh, plots_dir, true, z_min_trimmed, z_max_trimmed, false);
    
    %Window score plot
fighandle = RNA_Threshold_Common.drawWindowscorePlot(spots_sample(:,1), win_stdevs, spotsrun.ttune_wscorethresh, thresh);
saveas(fighandle, [plots_dir filesep 'autothresh_windowscore.png']);
close(fighandle);
    
%Save THIS module's run info (including paths, parameters, chosen threshold etc)
    %Save ztrimmed coords/spotcounts here too. Remove coord tables below
    %   threshold 10 to save space.
spotsrun.saveMe();

end