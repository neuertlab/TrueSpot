%
%%
function quant_results = RNAQuantPipe(param_struct, guimode)

use_nuc_mask = param_struct.use_nuc_mask; %0 = 2d, 1 = lo, 2 = mid, 3 = hi
gaussrad = param_struct.gaussian_radius;

quant_results = [];
RNA_Fisher_State.setGUIMode(guimode);
workdir = param_struct.outdir;
if isempty(workdir)
    workdir = pwd;
end

%-- Load RNA Spotsrun, if provided.
%If not provided, check for required parameters.
rnaspots_run = [];
tifpath = [];
if ~isempty(param_struct.runpath)
    rnaspots_run = RNASpotsRun.loadFrom(param_struct.runpath);
    
    if isempty(param_struct.tifpath)
        %Try to use the run's (though it may be no good).
        tifpath = rnaspots_run.paths.img_path;
    else
        %Override
        tifpath = param_struct.tifpath;
    end
    trgch = rnaspots_run.channels.rna_ch;
    totch = rnaspots_run.channels.total_ch;
    if gaussrad < 1
        gaussrad = rnaspots_run.options.dtune_gaussrad;
    end

    if isempty(param_struct.cellsegpath)
        param_struct.cellsegpath = rnaspots_run.paths.cellseg_path;
    end
    
    %Rethreshold, if needed or requested.
    if param_struct.rethresh | isempty(rnaspots_run.threshold_results)
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Rerunning spot thresholding..."), true);
        if param_struct.rethreshPreset > 0
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Using thresholding preset %d...", param_struct.rethreshPreset), true);
            rnaspots_run = RNAThreshold.applyPreset(rnaspots_run, param_struct.rethreshPreset);
        end
        rnaspots_run.threshold_results = RNAThreshold.runSavedParameters(rnaspots_run, 0);
        rnaspots_run.intensity_threshold = rnaspots_run.threshold_results.threshold;
        rnaspots_run = rnaspots_run.saveMeTo(param_struct.runpath);
        use_thresh = rnaspots_run.threshold_results.threshold;
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Rethresholding complete. New threshold: %d", use_thresh), true);
    else
        use_thresh = rnaspots_run.threshold_results.threshold;
    end
    
    if param_struct.man_thresh > 0
        %Manual override
        use_thresh = param_struct.man_thresh;
    end
    
    %Update z_adj and load coords
    if ~isempty(rnaspots_run.meta.idims_voxel)
        if rnaspots_run.meta.idims_voxel.z > 0
            param_struct.z_adj = rnaspots_run.meta.idims_voxel.z / rnaspots_run.meta.idims_voxel.x;
        end
    end
    
    %TODO handle new call_table structure
%     [~, coord_table] = rnaspots_run.loadCoordinateTable();
%     if ~isempty(coord_table)
%         th_idx = find(rnaspots_run.threshold_results.x(:,1) == use_thresh, 1);
%         if ~isempty(th_idx)
%             param_struct.coords = coord_table{th_idx,1};
%         else
%             param_struct.noclouds = true;
%             param_struct.no_refilter = false;
%         end
%     else 
%         param_struct.noclouds = true;
%         param_struct.no_refilter = false;
%     end
%     clear coord_table

    [~, call_table] = rnaspots_run.loadCallTable();
    if ~isempty(call_table)
        param_struct.coords = RNACoords.getThresholdCalls(call_table, use_thresh, true);
        if isempty(param_struct.coords)
            param_struct.noclouds = true;
            param_struct.no_refilter = false;
        end
    else
        RNA_Fisher_State.outputMessageLineStatic(sprintf("WARNING: Run call table could not be loaded!"), true);
        param_struct.noclouds = true;
        param_struct.no_refilter = false;
    end
    clear call_table

else
    %Need to make sure all the channel info is there...
    if param_struct.rna_channel < 1
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Target channel index is required for TIF loading!"), true);
        return;
    end
    if param_struct.channel_count < 1
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Channel count is required for TIF loading!"), true);
        return;
    end
    if param_struct.man_thresh < 1
        RNA_Fisher_State.outputMessageLineStatic(sprintf("If run file is not provided, manual threshold value must be specified."), true);
        return;
    end
    trgch = param_struct.rna_channel;
    totch = param_struct.channel_count;
    use_thresh = param_struct.man_thresh;
    %param_struct.noclouds = true; %Not supported w/o coord table (for now)
end

%Try to load the image (if not preloaded)
image_raw = [];
if isempty(param_struct.preloaded_image)
    if isempty(tifpath)
        RNA_Fisher_State.outputMessageLineStatic(sprintf("TIF image data are required! Exiting..."), true);
        return;
    end
    if ~isfile(tifpath)
        RNA_Fisher_State.outputMessageLineStatic(sprintf("TIF file %s could not be found! Exiting...", tifpath), true);
        return;
    end
    
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Loading TIF channel..."), true);
    [channels, ~] = LoadTif(tifpath, totch, [trgch], 1);
    image_raw = channels{trgch,1};
    clear channels;
else
    image_raw = param_struct.preloaded_image;
    param_struct.preloaded_image = [];
end

if isempty(image_raw)
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Raw image could not be loaded! Exiting..."), true);
	return;
end

%Try to load masks (if not preloaded)
cellmask = [];
if ~param_struct.nocells
    if isempty(param_struct.preloaded_cellmask)
        if ~isempty(param_struct.cellsegpath)
            cellmask_path = param_struct.cellsegpath;
        elseif ~isempty(param_struct.cellsegdir)
            cellmask_path = [param_struct.cellsegdir filesep 'Lab_' param_struct.cellsegname '.mat'];
        else
            RNA_Fisher_State.outputMessageLineStatic(sprintf("WARNING: Cell mask not provided. Quant will not be split by cell."), true);
            param_struct.nocells = true;
        end

        if ~param_struct.nocells
            if ~isfile(cellmask_path)
                RNA_Fisher_State.outputMessageLineStatic(sprintf("WARNING: Cell mask file %s could not be found! Cell seg will not be used.", cellmask_path), true);
                param_struct.nocells = true;
            else
                cellmask = CellSeg.openCellMask(cellmask_path);
            end
        end
    else
        cellmask = param_struct.preloaded_cellmask;
        param_struct.preloaded_cellmask = [];
    end
end

if isempty(cellmask) & ~param_struct.nocells
    RNA_Fisher_State.outputMessageLineStatic(sprintf("Cell mask could not be loaded! Cell seg will not be used."), true);
	param_struct.nocells = true;
end

nucmask = [];
if ~param_struct.nocells
    if isempty(param_struct.preloaded_nucmask)
        if ~isempty(param_struct.nucsegpath)
            nucmask_path = param_struct.nucsegpath;
        elseif ~isempty(param_struct.cellsegdir)
            nucmask_path = [param_struct.cellsegdir filesep 'nuclei_' param_struct.cellsegname '.mat'];
        elseif ~isempty(param_struct.cellsegpath)
            nucmask_path = param_struct.cellsegpath;
        else
            %Maybe update this to just use cell mask instead?
            RNA_Fisher_State.outputMessageLineStatic(sprintf("WARNING: Nuc mask not provided. Quant will not be split by cell."), true);
            param_struct.nocells = true;
        end

        if ~param_struct.nocells
            if ~isfile(nucmask_path)
                RNA_Fisher_State.outputMessageLineStatic(sprintf("WARNING: Nucleus mask file %s could not be found! Cell seg will not be used.", nucmask_path), true);
                param_struct.nocells = true;
            else
                nucmask = CellSeg.openNucMask(nucmask_path, use_nuc_mask, true); 
                if ndims(image_raw) < 3
                    nucmask = max(nucmask, [], 3);
                end
            end
        end
    else
        nucmask = param_struct.preloaded_nucmask;
        param_struct.preloaded_nucmask = [];
    end
end

if isempty(nucmask) & ~param_struct.nocells
    RNA_Fisher_State.outputMessageLineStatic(sprintf("WARNING: Nucleus mask could not be loaded! Cell seg will not be used."), true);
	param_struct.nocells = true;
end

if(param_struct.nocells)
    %Generate a dummy cell and nuc mask that's just the size of the whole
    %image.
    idim_x = size(image_raw,2);
    idim_y = size(image_raw,1);
    if ndims(image_raw) > 2
        idim_z = size(image_raw,3);
    else
        idim_z = 1;
    end
    
    cellmask = uint16(ones(idim_y,idim_x));
    nucmask = true(idim_y,idim_x,idim_z);
end

%Load coord table, if needed
if param_struct.no_refilter
    if isempty(param_struct.coords)
        %Look for table
        if ~isempty(param_struct.coord_tbl_path)
            if isfile(param_struct.coord_tbl_path)
                if endsWith(param_struct.coord_tbl_path, '.mat')
                    %Try to load it as an output of this tool
                    finfo = who('-file', param_struct.coord_tbl_path);
                    if ~isempty(find(ismember(finfo, 'coord_table'),1))
                        load(param_struct.coord_tbl_path, 'coord_table');
                        if ~isempty(coord_table)
                            if iscell(coord_table)
                                %We have to find the right threshold...
                                %Check for a spot table...
                                spot_table_path = replace(param_struct.coord_tbl_path, 'coordTable', 'spotTable');
                                if isfile(spot_table_path)
                                    load(spot_table_path, 'spot_table');
                                    thidx = RNAUtils.findThresholdIndex(use_thresh, transpose(spot_table(:,1)));
                                    param_struct.coords = coord_table{thidx,1};
                                    clear spot_table;
                                else
                                    %Just use as index...
                                    param_struct.coords = coord_table{use_thresh,1};
                                end
                            else
                                param_struct.coords = coord_table;
                            end
                            clear coord_table
                        end
                    elseif ~isempty(find(ismember(finfo, 'call_table'),1))
                        load(param_struct.coord_tbl_path, 'call_table');
                        param_struct.coords = RNACoords.getThresholdCalls(call_table, use_thresh, true);
                        clear call_table 
                    else
                        RNA_Fisher_State.outputMessageLineStatic(sprintf("Refilter skip requested, but coord table file not recognized! Doing refilter..."), true);
                        param_struct.no_refilter = false;
                    end
                else
                    %Try to load it as a table of x,y,z
                    param_struct.coords = readmatrix(param_struct.coord_tbl_path);
                end
            else
                RNA_Fisher_State.outputMessageLineStatic(sprintf("Refilter skip requested, but coord table path is invalid! Doing refilter..."), true);
                param_struct.no_refilter = false;
            end
        else
            RNA_Fisher_State.outputMessageLineStatic(sprintf("Refilter skip requested, but no coord table path provided! Doing refilter..."), true);
            param_struct.no_refilter = false;
        end
    end
end

%Setup proper parameter struct
quant_results = RNAQuant.genRNAQuantInfoStruct();
quant_results.img_raw = image_raw;
quant_results.threshold = use_thresh;
quant_results.t_coord_table = param_struct.coords;
quant_results.workdir = workdir;
quant_results.cell_mask = cellmask;
quant_results.nuc_mask = nucmask;
quant_results.do_clouds = ~param_struct.noclouds;
quant_results.do_refilter = ~param_struct.no_refilter;
quant_results.z_adj = param_struct.z_adj;
quant_results.workers = param_struct.workers;
quant_results.dbgcell = param_struct.dbgcell;
quant_results.no_bkg_subtract = param_struct.no_bkg_subtract;

if ~isempty(quant_results.t_coord_table)
   quant_results.t_coord_table = int32(quant_results.t_coord_table);
end

if gaussrad < 1; gaussrad = 7; end
quant_results.gaussian_radius = gaussrad;
quant_results.small_obj_size = param_struct.small_obj_size;
quant_results.spotzoom_r_xy = param_struct.spotzoom_r_xy;
quant_results.spotzoom_r_z = param_struct.spotzoom_r_z;

if ~quant_results.do_refilter
    if isempty(quant_results.t_coord_table)
        RNA_Fisher_State.outputMessageLineStatic(sprintf("Coordinate table is required if refiltering is to be skipped! Exiting..."), true);
        quant_results = [];
        return;
    end
end

%Run
quant_results = RNAQuant.FitRNA(quant_results);

end