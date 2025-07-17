%
%%

%Debug levels: 
%   0 - None
%   1 - Regular verbosity
%   2 - Regular verbosity + output plots
%   3 - High verbosity + output plots

%%
function spotsrun = RNA_Pipeline_Core(spotsrun, debug_lvl, preloaded_imgs, thread_request)

if nargin > 2
    bPreloaded = ~isempty(preloaded_imgs);
else
    bPreloaded = false;
end

if nargin < 4
    thread_request = 1;
end

if isempty(spotsrun)
    return;
end

spotsrun.options.threads = thread_request;
spotsrun.options.debug_level = debug_lvl;

%Debug print
if debug_lvl > 0
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Running RNASpots with the following parameters..."), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("img_name = %s", spotsrun.img_name), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("tif_path = %s", spotsrun.paths.img_path), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("rna_ch = %d", spotsrun.channels.rna_ch), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("light_ch = %d", spotsrun.channels.light_ch), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("total_ch = %d", spotsrun.channels.total_ch), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("out_dir = %s", spotsrun.paths.out_dir), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("out_stem = %s", spotsrun.paths.out_namestem), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("t_min = %d", spotsrun.options.t_min), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("t_max = %d", spotsrun.options.t_max), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("ztrim = %d", spotsrun.dims.ztrim), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("cellseg_path = %s", spotsrun.paths.cellseg_path), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("ctrl_path = %s", spotsrun.paths.ctrl_img_path), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("ctrl_ch = %d", spotsrun.channels.ctrl_ch), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("ctrl_chcount = %d", spotsrun.channels.ctrl_chcount), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("overwrite_output = %d", spotsrun.options.overwrite_output), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("thread_request = %d", thread_request), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Use preloaded images? = %d", bPreloaded), false);
end

if debug_lvl > 2
	tif_v = 2;
elseif debug_lvl == 0
	tif_v = 0;
else
	tif_v = 1;
end

if debug_lvl > 1
    spotsrun.options.save_maxproj = true;
end

%Declare some vars so it doesn't complain later...
sample_tif = [];
sample_light_ch = [];

%Do background extraction if arguments provided.
%!! Don't redo if target exists and overwrite output is false!
bkg_mask_dir = [spotsrun.paths.out_dir filesep 'bkgmask'];
spotsrun.paths.bkg_mask_path = [bkg_mask_dir filesep spotsrun.img_name '_bkg'];
if (~bPreloaded & isempty(spotsrun.paths.ctrl_img_path)) | (bPreloaded & isempty(preloaded_imgs.dat_rna_control))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control not provided. Will attempt to use background."), true);
    if ~spotsrun.options.overwrite_output && isfile([spotsrun.paths.bkg_mask_path '.mat'])
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Background mask already exists at %s! Skipping extraction...", spotsrun.paths.bkg_mask_path), true);
    else
        if (~isempty(spotsrun.paths.cellseg_path)) & (spotsrun.channels.light_ch > 0)
            %fprintf("Extracting background...\n");
            %Make sure cellseg data file exists
            if (isfile(spotsrun.paths.cellseg_path))
                if bPreloaded
                    sample_light_ch = preloaded_imgs.dat_trans_sample;
                else
                    [spotsrun, sample_tif] = spotsrun.loadSampleTif(tif_v);
                    sample_light_ch = sample_tif{spotsrun.channels.light_ch,1};
                end
                if ~isfolder(bkg_mask_dir)
                    mkdir(bkg_mask_dir);
                end

                Bkg_Mask_Core(sample_light_ch, spotsrun.paths.cellseg_path, spotsrun.paths.bkg_mask_path, debug_lvl > 1);
                %Main_BackgroundMask(spotsrun.tif_path, spotsrun.cellseg_path, spotsrun.bkg_path, spotsrun.total_ch, spotsrun.light_ch, true);
            else
                RNA_Fisher_State.outputMessageLineStatic(sprintf("Cellseg file %s does not exist! Skipping background extraction...", spotsrun.paths.cellseg_path), true);
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

if spotsrun.options.use_max_proj 
    strat = 'max_proj';
else
    strat = 'all_3d'; 
end

if ~isfolder(spotsrun.paths.out_dir)
    mkdir(spotsrun.paths.out_dir);
end
if isempty(spotsrun.paths.out_namestem)
    spotsrun.paths.out_namestem = [spotsrun.img_name '_' strat];
end
sample_outstem = spotsrun.getFullOutStem();

%Check for existing spotsrun file to copy over some params
spotsrun_path = [sample_outstem '_rnaspotsrun.mat'];
old_spotsrun = [];
if isfile(spotsrun_path)
    old_spotsrun = RNASpotsRun.loadFrom(spotsrun_path);
end

runme = true;
if ~spotsrun.options.overwrite_output
    if isfile([sample_outstem '_callTable.mat'])
        fprintf("Spot detection output at %s already exists! Skipping spot detection...\n", sample_outstem);
        runme = false;
    end
end

if runme
    %Load image channel
    if bPreloaded
        spotsrun.dims.idims_sample = struct('x', 0, 'y', 0, 'z', 0);
        spotsrun.dims.idims_sample.x = size(preloaded_imgs.dat_rna_sample,2);
        spotsrun.dims.idims_sample.y = size(preloaded_imgs.dat_rna_sample,1);
        if ndims(preloaded_imgs.dat_rna_sample) > 2
            spotsrun.dims.idims_sample.z = size(preloaded_imgs.dat_rna_sample,3);
        else 
            spotsrun.dims.idims_sample.z = 1;
        end
        sample_rna_ch = double(preloaded_imgs.dat_rna_sample);
    else
        %Load sample TIF
        if isempty(sample_tif)
            [spotsrun, sample_tif] = spotsrun.loadSampleTif(tif_v);
        end
        sample_rna_ch = sample_tif{spotsrun.channels.rna_ch,1};
        spotsrun.dims.idims_sample = struct('x', 0, 'y', 0, 'z', 0);
        spotsrun.dims.idims_sample.x = size(sample_rna_ch,2);
        spotsrun.dims.idims_sample.y = size(sample_rna_ch,1);
        if ndims(sample_rna_ch) > 2
            spotsrun.dims.idims_sample.z = size(sample_rna_ch,3);
        else 
            spotsrun.dims.idims_sample.z = 1;
        end
    end
    
    spotdec = RNA_Threshold_SpotDetector;
    [img_f] = spotdec.run_spot_detection_pre(sample_rna_ch, sample_outstem, spotsrun.options.deadpix_detect, spotsrun.options.dtune_gaussrad, spotsrun.options.save_maxproj);
    
    %Clear original image to free (a ton of) memory
    if ~bPreloaded & ~isempty(sample_tif)
        clear sample_tif;
    end
    if ~bPreloaded & ~isempty(sample_light_ch)
        clear sample_light_ch;
    end
    if ~bPreloaded
        chdat_path = [sample_outstem '_samplech.mat'];
        save(chdat_path, 'sample_rna_ch', '-v7.3');
        clear sample_rna_ch; 
    end
    
    %Suggest scan min or max if not set
    if spotsrun.options.t_max < 1 | spotsrun.options.t_min < 1
        [th_min, th_max] = RNA_Threshold_Common.suggestScanThreshold(img_f);
        if spotsrun.options.t_min < 1
            spotsrun.options.t_min = th_min; 
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Threshold scan min auto-set to: %d", spotsrun.options.t_min), true);
        end
        if spotsrun.options.t_max < 1
            spotsrun.options.t_max = th_max; 
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Threshold scan max auto-set to: %d", spotsrun.options.t_max), true);
        end
    end
    
    if debug_lvl > 0
        RNA_Fisher_State.outputMessageLineStatic('Pre-processing complete. Now continuing to detection...', true);
    end

    zMin = spotsrun.dims.z_min;
    zMax = spotsrun.dims.z_max;
    if (spotsrun.dims.z_min < 1) & (spotsrun.dims.ztrim > 0)
        zMin = spotsrun.dims.ztrim + 1;
    end
    if (spotsrun.dims.z_max < 1) & (spotsrun.dims.ztrim > 0)
        zMax = spotsrun.dims.idims_sample.z - spotsrun.dims.ztrim;
    end
    [spotsrun.dims.z_min, spotsrun.dims.z_max, call_table] = spotdec.run_spot_detection_main(img_f, sample_outstem, strat, spotsrun.options.t_min, spotsrun.options.t_max, zMin, zMax, (debug_lvl > 0), thread_request);
    %spotsrun.dims.ztrim_auto = auto_zt;
    %spotsrun.t_min = new_th_min;
    if spotsrun.options.use_max_proj
        spotsrun.dims.ztrim = 0;
        spotsrun.dims.ztrim_auto = 0;
    end

    if ~isempty(call_table)
        if ~bPreloaded
            load(chdat_path, 'sample_rna_ch');
            delete(chdat_path);
        end

        if spotsrun.options.use_max_proj 
            max_proj = max(double(sample_rna_ch),[],3);
            call_table{:, 'intensity'} = single(max_proj(call_table{:, 'coord_1d'}));
            clear max_proj
        else
            call_table{:, 'intensity'} = sample_rna_ch(call_table{:, 'coord_1d'});
        end

        RNASpotsRun.saveCallTable(call_table, sample_outstem);
        clear call_table;

        if ~bPreloaded
            clear sample_rna_ch; 
        end
    end

    clear img_f;
else
    %Copy over params from old run, if available.
    if ~isempty(old_spotsrun)
        fprintf('Previous run was found. Copying back spot detection parameters!\n');
        spotsrun.channels = old_spotsrun.channels;
        spotsrun.options.dtune_gaussrad = old_spotsrun.options.dtune_gaussrad;
        spotsrun.options.deadpix_detect = old_spotsrun.options.deadpix_detect;
        spotsrun.options.use_max_proj = old_spotsrun.options.use_max_proj;
        spotsrun.options.t_min = old_spotsrun.options.t_min;
        spotsrun.options.t_max = old_spotsrun.options.t_max;
        spotsrun.dims.idims_sample = old_spotsrun.dims.idims_sample;
        spotsrun.dims.idims_ctrl = old_spotsrun.dims.idims_ctrl;
    end
end

%Run spot detect on control (if a tif path was provided)
%   Detect if input control path is a tif. If not, assume the input is a spot
%       detect results path stem.
%!! Don't redo if target exists and overwrite output is false!
ctrl_stem = [];
if ~isempty(spotsrun.paths.ctrl_img_path)
    if endsWith(spotsrun.paths.ctrl_img_path, ".tif")
        
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Running spot detect on control image... (This may take a few hours on large files)"), true);
        %spotsrun.ctrl_stem = Main_RNASpotDetect([spotsrun.img_name '_Control'], spotsrun.ctrl_path, spotsrun.out_dir,...
        %    spotsrun.ctrl_ch, spotsrun.ctrl_chcount, spotsrun.t_min, spotsrun.t_max, true, spotsrun.overwrite_output);
 
        %outdir = [spotsrun.out_dir filesep strat];
        if isempty(spotsrun.paths.ctrl_out_namestem)
            spotsrun.paths.ctrl_out_namestem = [spotsrun.img_name '_Control_' strat];
        end
        ctrl_stem = spotsrun.getFullCtrlOutStem();
        runme = true;
        if ~spotsrun.options.overwrite_output
            %TODO check this run params against saved - if not match, throw
            %an error and return!!
            if isfile([ctrl_stem '_callTable.mat'])
                fprintf("Spot detection output at %s already exists! Skipping spot detection...\n", ctrl_stem);
                runme = false;
            end
        end

        if runme
            if bPreloaded
                spotsrun.dims.idims_control = struct('x', 0, 'y', 0, 'z', 0);
                spotsrun.dims.idims_control.x = size(preloaded_imgs.dat_rna_control,2);
                spotsrun.dims.idims_control.y = size(preloaded_imgs.dat_rna_control,1);
                if ndims(preloaded_imgs.dat_rna_control) > 2
                    spotsrun.dims.idims_control.z = size(preloaded_imgs.dat_rna_control,3);
                else 
                    spotsrun.dims.idims_control.z = 1;
                end
                ctrl_image_channel = preloaded_imgs.dat_rna_control;
            else
                if ~isfile(spotsrun.paths.ctrl_img_path)
                    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control path %s does not exist. Aborting...", spotsrun.paths.ctrl_img_path), true);
                    return;
                end

                [spotsrun, ctrl_image_channel] = spotsrun.loadControlChannel(tif_v);
            end

            zMin = spotsrun.dims.z_min;
            zMax = spotsrun.dims.z_max;
            if (spotsrun.dims.z_min < 1) & (spotsrun.dims.ztrim > 0)
                zMin = spotsrun.dims.ztrim + 1;
            end
            if (spotsrun.dims.z_max < 1) & (spotsrun.dims.ztrim > 0)
                zMax = spotsrun.dims.idims_sample.z - spotsrun.dims.ztrim;
            end

            spotdec = RNA_Threshold_SpotDetector;
            [ctrl_f] = spotdec.run_spot_detection_pre(ctrl_image_channel, ctrl_stem, spotsrun.options.deadpix_detect, spotsrun.options.dtune_gaussrad, spotsrun.options.save_maxproj);
            [~, ~, call_table] = spotdec.run_spot_detection_main(ctrl_f, ctrl_stem, strat, spotsrun.options.t_min, spotsrun.options.t_max, zMin, zMax, (debug_lvl > 0), thread_request);
            if ~isempty(call_table)
                call_table{:, 'intensity'} = ctrl_image_channel(call_table{:, 'coord_1d'});
                RNASpotsRun.saveCallTable(call_table, ctrl_stem);
            end
            if ~bPreloaded; clear ctrl_image_channel; end
            clear call_table;
            clear ctrl_f;
        end
    else
        %spotsrun.ctrl_stem = spotsrun.ctrl_path;
    end
else
    spotsrun.paths.ctrl_out_namestem = [];
end

%Filter spots for bkg (if mask exists)
spotsrun.paths.bkg_filter_stem = [sample_outstem '_bkgmasked'];
if isfile([spotsrun.paths.bkg_mask_path '.mat'])
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Filtering detected coordinates through background mask for control..."), true);
    
    spotsrun.updateBackgroundFilteredCoords();
    
    %Set as control if there isn't one.
    if isempty(spotsrun.paths.ctrl_out_namestem)
        [~,ff,~] = fileparts(spotsrun.paths.bkg_filter_stem);
        [~,dd,~] = fileparts(bkg_mask_dir);
        spotsrun.paths.ctrl_out_namestem = [dd filesep ff];
    end
end
spotsrun.saveMe();

%Get image dimensions (for ztrim)
idims = spotsrun.dims.idims_sample;
if isempty(idims) | idims.x < 1
    [spotsrun,~] = spotsrun.loadSampleTif(tif_v);
    idims = spotsrun.dims.idims_sample;
end
spotsrun = spotsrun.updateZTrimParams();

if spotsrun.dims.z_min_apply > 1 | spotsrun.dims.z_max_apply < idims.z
    RNA_Fisher_State.outputMessageLineStatic(sprintf("ZTrim will be applied."), true);
end

%Load spot count tables and apply ztrim
%Ztrim is now handled in RNASpotsRun

%Detect threshold
% RNA_Fisher_State.outputMessageLineStatic(sprintf("Finding a good threshold..."), true);
% [thresh, win_scores, score_thresh, scanst] = RNA_Threshold_Common.estimateThreshold(spots_sample, spots_control, spotsrun.ttune_winsize, 0.5, spotsrun.ttune_madfactor);
% RNA_Fisher_State.outputMessageLineStatic(sprintf("Auto threshold selected: %d", thresh), true);
% spotsrun.intensity_threshold = thresh;

%Updated interface:
RNA_Fisher_State.outputMessageLineStatic(sprintf("Finding a good threshold..."), true);
spotsrun.threshold_results = RNAThreshold.runSavedParameters(spotsrun, tif_v);
if spotsrun.threshold_results.lowNoiseFlag
    spotsrun.threshold_results.threshold = 1;
end
spotsrun.intensity_threshold = spotsrun.threshold_results.threshold;

%Run through all threshold presets
RNA_Fisher_State.outputMessageLineStatic(sprintf("Testing all threshold presets..."), true);
thpresets = RNAThreshold.loadPresets();
thpreCount = size(thpresets, 2);
spotsrun.th_alt.preset_count = thpreCount;
spotsrun.th_alt.altThResults = cell(1, thpreCount);
spotsrun.th_alt.thPresetSugg = NaN(thpreCount, 7);

[spotsrun, sptTableSmpl] = spotsrun.loadSpotsTable();
[spotsrun, sptTableCtrl] = spotsrun.loadControlSpotsTable();
for thp = 1:thpreCount
    thpreRes = RNAThreshold.runWithPreset(sptTableSmpl, sptTableCtrl, thp);
    if ~isempty(thpreRes)
        spotsrun.th_alt.altThResults{thp} = thpreRes;
        spotsrun.th_alt.thPresetSugg(thp, 1) = thpreRes.threshold;
        if isfield(thpreRes, 'mean_w')
            spotsrun.th_alt.thPresetSugg(thp, 2) = thpreRes.mean_w - thpreRes.std_w;
            spotsrun.th_alt.thPresetSugg(thp, 3) = thpreRes.mean_w;
            spotsrun.th_alt.thPresetSugg(thp, 4) = thpreRes.mean_w + thpreRes.std_w;
        end
        if isfield(thpreRes, 'median_w')
            spotsrun.th_alt.thPresetSugg(thp, 5) = thpreRes.median_w - thpreRes.mad_w;
            spotsrun.th_alt.thPresetSugg(thp, 6) = thpreRes.median_w;
            spotsrun.th_alt.thPresetSugg(thp, 7) = thpreRes.median_w + thpreRes.mad_w;
        end
    end
end
clear thp thpresets thpreCount sptTableSmpl sptTableCtrl thpreRes

%Print plots & image representation
    %Spot plots (log and linear scale) - w/ chosen threshold marked
    %Max projection of sample and control w circled spots (also max proj)
   
if debug_lvl > 1
     plots_dir = [spotsrun.paths.out_dir filesep 'plots'];
     RNA_Fisher_State.outputMessageLineStatic(sprintf("Now generating plots..."), true);
     if ~isfolder(plots_dir); mkdir(plots_dir); end
% 
%     %These are all very broken at the moment. I wouldn't use them.
%     if ~isempty(ctrl_stem)
%         probeNames = [{spotsrun.img_name};
%                     {'Control'}];
%         RNA_Threshold_Plotter.plotPreprocessedData(sample_outstem, [{ctrl_stem}], probeNames,...
%                 spotsrun.intensity_threshold, plots_dir, true, spotsrun.dims.z_min_apply, spotsrun.dims.z_max_apply, false);
%     else
%         probeNames = [{spotsrun.img_name}];
%         RNA_Threshold_Plotter.plotPreprocessedData(sample_outstem, [], probeNames,...
%                 spotsrun.intensity_threshold, plots_dir, true, spotsrun.dims.z_min_apply, spotsrun.dims.z_max_apply, false);
%     end

    %Window score plots and thresholding results.
    figh = RNAThreshold.resultPlotCombine(spotsrun, 615);
    if ~isempty(figh)
        saveas(figh, [plots_dir filesep 'thres_all.png']);
        close(figh);
    end
    
    fig_handles = RNAThreshold.resultPlotsIndiv(spotsrun, 100);
    fig_count = size(fig_handles,2);
    for i = 1:fig_count
        saveas(fig_handles(1,i), [plots_dir filesep sprintf('thres_plot_%02d.png', i)]);
        close(fig_handles(1,i));
    end
end
    
%Save THIS module's run info (including paths, parameters, chosen threshold etc)
    %Save ztrimmed coords/spotcounts here too. Remove coord tables below
    %   threshold 10 to save space.
spotsrun.saveMe();

%Dump tables, if applicable.
if ~isempty(spotsrun.paths.params_out_path)
    if debug_lvl > 0
        RNA_Fisher_State.outputMessageLineStatic("Dumping run parameters to text...", true);
    end
    spotsrun.toTextFile(spotsrun.paths.params_out_path);
end

if ~isempty(spotsrun.paths.csv_out_path)
    if debug_lvl > 0
        RNA_Fisher_State.outputMessageLineStatic("Dumping results to csv...", true);
    end
    [spotsrun, ~, call_table] = spotsrun.loadZTrimmedTables_Sample();

    th_min = spotsrun.options.t_min;

    if spotsrun.options.csv_output_level == 2
        th_min = min(RNAThreshold.getAllThresholdSuggestions(spotsrun.threshold_results), []);
    elseif spotsrun.options.csv_output_level == 3
        th_min = spotsrun.intensity_threshold;
    end

    RNADetection.callTable2csv(spotsrun.paths.csv_out_path, call_table, th_min, spotsrun.options.csv_zero_based_coords);
end

if debug_lvl > 0
    RNA_Fisher_State.outputMessageLineStatic("Core pipeline done!", true);
end

end