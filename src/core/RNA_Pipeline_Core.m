%
%%

function RNA_Pipeline_Core(spotsrun, verbosity)
%TestCommmentBK1
%Debug print
if verbosity > 0
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
end

%Load sample TIF
[spotsrun, sample_tif] = spotsrun.loadSampleTif(verbosity);
sample_rna_ch = sample_tif{spotsrun.rna_ch,1};

%Do background extraction if arguments provided.
%!! Don't redo if target exists and overwrite output is false!
bkg_mask_dir = [spotsrun.out_dir filesep 'bkgmask'];
mkdir(bkg_mask_dir);
spotsrun.bkg_path = [bkg_mask_dir filesep spotsrun.img_name '_bkg'];
if isempty(spotsrun.ctrl_path)
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control not provided. Will attempt to use background."), true);
    if ~spotsrun.overwrite_output && isfile([spotsrun.bkg_path '.mat'])
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Background mask already exists at %s! Skipping extraction...", spotsrun.bkg_path), true);
    else
        if (~isempty(spotsrun.cellseg_path)) && (spotsrun.light_ch > 0)
            %fprintf("Extracting background...\n");
            %Make sure cellseg data file exists
            if (isfile(spotsrun.cellseg_path))
                sample_light_ch = sample_tif{spotsrun.light_ch,1};
                Bkg_Mask_Core(sample_light_ch, spotsrun.cellseg_path, spotsrun.bkg_path, true);
                %Main_BackgroundMask(spotsrun.tif_path, spotsrun.cellseg_path, spotsrun.bkg_path, spotsrun.total_ch, spotsrun.light_ch, true);
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
%spotsrun.out_stem = Main_RNASpotDetect(spotsrun.img_name, spotsrun.tif_path, spotsrun.out_dir,...
%    spotsrun.rna_ch, spotsrun.total_ch, spotsrun.t_min, spotsrun.t_max, true, spotsrun.overwrite_output);
strat = 'all_3d';
%outdir = [spotsrun.out_dir filesep strat];
mkdir(spotsrun.out_dir);
spotsrun.out_stem = [spotsrun.out_dir filesep spotsrun.img_name '_' strat];
runme = true;
if ~spotsrun.overwrite_output
    if isfile([spotsrun.out_stem '_coordTable.mat'])
        fprintf("Spot detection output at %s already exists! Skipping spot detection...\n", spotsrun.out_stem);
        runme = false;
    end
end

if runme
    spotdec = RNA_Threshold_SpotDetector;
    spotdec.run_spot_detection(sample_rna_ch, spotsrun.out_stem, strat, spotsrun.t_min, spotsrun.t_max, true, (verbosity > 0));
end

%Run spot detect on control (if a tif path was provided)
%   Detect if input control path is a tif. If not, assume the input is a spot
%       detect results path stem.
%!! Don't redo if target exists and overwrite output is false!
spotsrun.ctrl_stem = [];
if ~isempty(spotsrun.ctrl_path)
    if endsWith(spotsrun.ctrl_path, ".tif")
        if ~isfile(spotsrun.ctrl_path)
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Control path %s does not exist. Aborting...", spotsrun.ctrl_path), true);
            return;
        end
        
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Running spot detect on control image... (This may take a few hours on large files)"), true);
        %spotsrun.ctrl_stem = Main_RNASpotDetect([spotsrun.img_name '_Control'], spotsrun.ctrl_path, spotsrun.out_dir,...
        %    spotsrun.ctrl_ch, spotsrun.ctrl_chcount, spotsrun.t_min, spotsrun.t_max, true, spotsrun.overwrite_output);
        [spotsrun, ctrl_image_channel] = spotsrun.loadControlChannel(verbosity);
        %outdir = [spotsrun.out_dir filesep strat];
        outstem = [spotsrun.out_dir filesep spotsrun.img_name '_Control_' strat];
        runme = true;
        if ~spotsrun.overwrite_output
            if isfile([outstem '_coordTable.mat'])
                fprintf("Spot detection output at %s already exists! Skipping spot detection...\n", outstem);
                runme = false;
            end
        end

        if runme
            spotdec = RNA_Threshold_SpotDetector;
            spotdec.run_spot_detection(ctrl_image_channel, outstem, strat, spotsrun.t_min, spotsrun.t_max, true, (verbosity > 0));
        end
        spotsrun.ctrl_stem = outstem;
    else
        spotsrun.ctrl_stem = spotsrun.ctrl_path;
    end
end

[~, th_list] = spotsrun.loadThresholdTable();

%Filter spots for bkg (if mask exists)
spotsrun.bkg_filter_stem = [spotsrun.out_stem '_bkgmasked'];
if isfile([spotsrun.bkg_path '.mat'])
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Filtering detected coordinates through background mask for control..."), true);
    
    spotsrun.updateBackgroundFilteredCoords();
    
    %Set as control if there isn't one.
    if isempty(spotsrun.ctrl_stem)
        spotsrun.ctrl_stem = spotsrun.bkg_filter_stem;
    end
end
spotsrun.saveMe();

%Get image dimensions (for ztrim)
idims = spotsrun.idims_sample;
z_min_trimmed = 1;
z_max_trimmed = idims.z;

%Load spot count tables and apply ztrim
%TODO - ztrim not working - I think it's just trimming out everything XD
if spotsrun.ztrim > 0
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Applying z trim..."), true);
    %Create ztrim mask.
    ztrim_mask = true(idims.y, idims.x, idims.z);
    z_min_trimmed = spotsrun.ztrim;
    ztrim_mask(:,:,1:z_min_trimmed) = false;
    z_max_trimmed = idims.z - spotsrun.ztrim +  1;
    ztrim_mask(:,:,z_max_trimmed:idims.z) = false;
    
    %Load the tables
    [~, coords_sample] = spotsrun.loadCoordinateTable();
    
    %Mask the tables
    [spots_sample, trimmed_coords] = RNA_Threshold_Common.mask_spots(ztrim_mask, coords_sample, th_list);
    save([spotsrun.out_stem '_ztrim' num2str(spotsrun.ztrim)], 'trimmed_coords');
    
    %--- And repeat with the control
    if ~isfile([spotsrun.bkg_path '.mat'])
        %Need to retrieve dims again and make a new mask.
        idims = spotsrun.idims_ctrl;
        ztrim_mask = true(idims.y, idims.x, idims.z);
        zbot = spotsrun.ztrim;
        ztrim_mask(:,:,1:zbot) = false;
        ztop = idims.z - spotsrun.ztrim +  1;
        ztrim_mask(:,:,ztop:idims.z) = false;
    end
    [~, coords_control] = spotsrun.loadControlCoordinateTable();
    [spots_control, trimmed_coords] = RNA_Threshold_Common.mask_spots(ztrim_mask, coords_control, th_list);
    save([spotsrun.ctrl_stem '_ztrim' num2str(spotsrun.ztrim)], 'trimmed_coords');
else
    %Just load the coord and spot tables
    [~, spots_sample] = spotsrun.loadSpotsTable();
    [~, spots_control] = spotsrun.loadControlSpotsTable();
end

%Detect threshold
RNA_Fisher_State.outputMessageLineStatic(sprintf("Finding a good threshold..."), true);
[thresh, win_stdevs] = RNA_Threshold_Common.estimateThreshold(spots_sample, spots_control, spotsrun.ttune_winsize, 0.0, spotsrun.ttune_wscorethresh);
RNA_Fisher_State.outputMessageLineStatic(sprintf("Auto threshold selected: %d", thresh), true);
spotsrun.intensity_threshold = thresh;

%Print plots & image representation
    %Spot plots (log and linear scale) - w/ chosen threshold marked
    %Max projection of sample and control w circled spots (also max proj)
    
plots_dir = [spotsrun.out_dir filesep 'plots'];
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