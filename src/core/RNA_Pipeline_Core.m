%
%%

%Debug levels: 
%   0 - None
%   1 - Regular verbosity
%   2 - Regular verbosity + output plots
%   3 - High verbosity + output plots

%%
function spotsrun = RNA_Pipeline_Core(spotsrun, debug_lvl, preloaded_imgs, limitSaveSize, thread_request, use_max_proj)

if nargin > 2
    bPreloaded = ~isempty(preloaded_imgs);
else
    bPreloaded = false;
end

if nargin < 4
    limitSaveSize = true;
end

if nargin < 5
    thread_request = 1;
end

if nargin < 6
    use_max_proj = false;
end

%Debug print
if debug_lvl > 0
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
    %RNA_Fisher_State.outputMessageLineStatic(sprintf("ttune_winsize = %d", spotsrun.ttune_winsize), false);
    %RNA_Fisher_State.outputMessageLineStatic(sprintf("ttune_wscorethresh = %f", spotsrun.ttune_wscorethresh), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("overwrite_output = %d", spotsrun.overwrite_output), false);
    RNA_Fisher_State.outputMessageLineStatic(sprintf("limitSaveSize = %d", limitSaveSize), false);
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

%Declare some vars so it doesn't complain later...
sample_tif = [];
sample_light_ch = [];

%Do background extraction if arguments provided.
%!! Don't redo if target exists and overwrite output is false!
bkg_mask_dir = [spotsrun.out_dir filesep 'bkgmask'];
spotsrun.bkg_path = [bkg_mask_dir filesep spotsrun.img_name '_bkg'];
if (~bPreloaded & isempty(spotsrun.ctrl_path)) | (bPreloaded & isempty(preloaded_imgs.dat_rna_control))
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control not provided. Will attempt to use background."), true);
    if ~spotsrun.overwrite_output && isfile([spotsrun.bkg_path '.mat'])
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Background mask already exists at %s! Skipping extraction...", spotsrun.bkg_path), true);
    else
        if (~isempty(spotsrun.cellseg_path)) & (spotsrun.light_ch > 0)
            %fprintf("Extracting background...\n");
            %Make sure cellseg data file exists
            if (isfile(spotsrun.cellseg_path))
                if bPreloaded
                    sample_light_ch = preloaded_imgs.dat_trans_sample;
                else
                    [spotsrun, sample_tif] = spotsrun.loadSampleTif(tif_v);
                    sample_light_ch = sample_tif{spotsrun.light_ch,1};
                end
                mkdir(bkg_mask_dir);
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

if use_max_proj 
    strat = 'max_proj';
else
    strat = 'all_3d'; 
end
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
    %Load image channel
    if bPreloaded
        spotsrun.idims_sample = struct('x', 0, 'y', 0, 'z', 0);
        spotsrun.idims_sample.x = size(preloaded_imgs.dat_rna_sample,2);
        spotsrun.idims_sample.y = size(preloaded_imgs.dat_rna_sample,1);
        if ndims(preloaded_imgs.dat_rna_sample) > 2
            spotsrun.idims_sample.z = size(preloaded_imgs.dat_rna_sample,3);
        else 
            spotsrun.idims_sample.z = 1;
        end
        sample_rna_ch = double(preloaded_imgs.dat_rna_sample);
    else
        %Load sample TIF
        if isempty(sample_tif)
            [spotsrun, sample_tif] = spotsrun.loadSampleTif(tif_v);
        end
        sample_rna_ch = sample_tif{spotsrun.rna_ch,1};
        spotsrun.idims_sample = struct('x', 0, 'y', 0, 'z', 0);
        spotsrun.idims_sample.x = size(sample_rna_ch,2);
        spotsrun.idims_sample.y = size(sample_rna_ch,1);
        if ndims(sample_rna_ch) > 2
            spotsrun.idims_sample.z = size(sample_rna_ch,3);
        else 
            spotsrun.idims_sample.z = 1;
        end
    end
    
    spotdec = RNA_Threshold_SpotDetector;
    [img_f] = spotdec.run_spot_detection_pre(sample_rna_ch, spotsrun.out_stem, spotsrun.deadpix_detect, spotsrun.dtune_gaussrad);
    
    %Clear original image to free (a ton of) memory
    if ~isempty(sample_tif)
        clear sample_tif;
    end
    if ~isempty(sample_light_ch)
        clear sample_light_ch;
    end
    clear sample_rna_ch;
    
    %Suggest scan min or max if not set
    if spotsrun.t_max < 1 | spotsrun.t_min < 1
        [th_min, th_max] = RNA_Threshold_Common.suggestScanThreshold(img_f);
        if spotsrun.t_min < 1
            spotsrun.t_min = th_min; 
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Threshold scan min auto-set to: %d", spotsrun.t_min), true);
        end
        if spotsrun.t_max < 1
            spotsrun.t_max = th_max; 
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Threshold scan max auto-set to: %d", spotsrun.t_max), true);
        end
    end
    
    [auto_zt, new_th_min] = spotdec.run_spot_detection_main(img_f, spotsrun.out_stem, strat, spotsrun.t_min, spotsrun.t_max, spotsrun.ztrim, limitSaveSize, (debug_lvl > 0), thread_request);
    spotsrun.ztrim_auto = auto_zt;
    spotsrun.t_min = new_th_min;
    clear img_f;
end

%Run spot detect on control (if a tif path was provided)
%   Detect if input control path is a tif. If not, assume the input is a spot
%       detect results path stem.
%!! Don't redo if target exists and overwrite output is false!
spotsrun.ctrl_stem = [];
if ~isempty(spotsrun.ctrl_path)
    if endsWith(spotsrun.ctrl_path, ".tif")
        
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Running spot detect on control image... (This may take a few hours on large files)"), true);
        %spotsrun.ctrl_stem = Main_RNASpotDetect([spotsrun.img_name '_Control'], spotsrun.ctrl_path, spotsrun.out_dir,...
        %    spotsrun.ctrl_ch, spotsrun.ctrl_chcount, spotsrun.t_min, spotsrun.t_max, true, spotsrun.overwrite_output);
 
        %outdir = [spotsrun.out_dir filesep strat];
        outstem = [spotsrun.out_dir filesep spotsrun.img_name '_Control_' strat];
        runme = true;
        if ~spotsrun.overwrite_output
            %TODO check this run params against saved - if not match, throw
            %an error and return!!
            if isfile([outstem '_coordTable.mat'])
                fprintf("Spot detection output at %s already exists! Skipping spot detection...\n", outstem);
                runme = false;
            end
        end

        if runme
            if bPreloaded
                spotsrun.idims_control = struct('x', 0, 'y', 0, 'z', 0);
                spotsrun.idims_control.x = size(preloaded_imgs.dat_rna_control,2);
                spotsrun.idims_control.y = size(preloaded_imgs.dat_rna_control,1);
                if ndims(preloaded_imgs.dat_rna_control) > 2
                    spotsrun.idims_control.z = size(preloaded_imgs.dat_rna_control,3);
                else 
                    spotsrun.idims_control.z = 1;
                end
                ctrl_image_channel = preloaded_imgs.dat_rna_control;
            else
                if ~isfile(spotsrun.ctrl_path)
                    RNA_Fisher_State.outputMessageLineStatic(sprintf("Control path %s does not exist. Aborting...", spotsrun.ctrl_path), true);
                    return;
                end
                
                [spotsrun, ctrl_image_channel] = spotsrun.loadControlChannel(tif_v);
            end
            
            spotdec = RNA_Threshold_SpotDetector;
            [ctrl_f] = spotdec.run_spot_detection_pre(ctrl_image_channel, outstem, spotsrun.deadpix_detect, spotsrun.dtune_gaussrad);
            clear ctrl_image_channel;
            spotdec.run_spot_detection_main(ctrl_f, outstem, strat, spotsrun.t_min, spotsrun.t_max, spotsrun.ztrim, limitSaveSize, (debug_lvl > 0), thread_request);
            clear ctrl_f;
        end
        spotsrun.ctrl_stem = outstem;
    else
        spotsrun.ctrl_stem = spotsrun.ctrl_path;
    end
end

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
if isempty(idims)
    [spotsrun,~] = spotsrun.loadSampleTif(tif_v);
    idims = spotsrun.idims_sample;
end
spotsrun = spotsrun.updateZTrimParams();

if spotsrun.z_min_apply > 1 | spotsrun.z_max_apply < idims.z
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
spotsrun.intensity_threshold = spotsrun.threshold_results.threshold;

%Print plots & image representation
    %Spot plots (log and linear scale) - w/ chosen threshold marked
    %Max projection of sample and control w circled spots (also max proj)
   
if debug_lvl > 1
    plots_dir = [spotsrun.out_dir filesep 'plots'];
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Now generating plots..."), true);

    if ~isempty(spotsrun.ctrl_stem)
        probeNames = [{spotsrun.img_name};
                    {'Control'}];
        RNA_Threshold_Plotter.plotPreprocessedData(spotsrun.out_stem, [{spotsrun.ctrl_stem}], probeNames,...
                spotsrun.intensity_threshold, plots_dir, true, spotsrun.z_min_apply, spotsrun.z_max_apply, false);
    else
        probeNames = [{spotsrun.img_name}];
        RNA_Threshold_Plotter.plotPreprocessedData(spotsrun.out_stem, [], probeNames,...
                spotsrun.intensity_threshold, plots_dir, true, spotsrun.z_min_apply, spotsrun.z_max_apply, false);
    end

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

end