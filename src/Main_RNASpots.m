%%
%%

function [bkg_path, out_stem, intensity_threshold] = Main_RNASpots(img_name, tif_path, rna_ch, light_ch, total_ch,...
    out_dir, t_min, t_max, ztrim, cellseg_path, ctrl_path, ctrl_ch, ctrl_chcount, ttune_winsize, ttune_wscorethresh,...
    overwrite_output)

addpath('./core');
intensity_threshold = 0;

%Check arguments
%Force to numbers
rna_ch = Force2Num(rna_ch);
light_ch = Force2Num(light_ch);
total_ch = Force2Num(total_ch);
ctrl_ch = Force2Num(ctrl_ch);
ctrl_chcount = Force2Num(ctrl_chcount);
t_min = Force2Num(t_min);
t_max = Force2Num(t_max);
ztrim = Force2Num(ztrim);
ttune_winsize = Force2Num(ttune_winsize);
ttune_wscorethresh = Force2Num(ttune_wscorethresh);

%Check required arguments
if (isempty(tif_path)) || (~ischar(tif_path))
    fprintf("Input image path is required.\n");
    return;
end
if (~isempty(out_dir)) && (~ischar(out_dir))
    fprintf("Output directory argument is invalid.\n");
    return;
end
if (~isempty(cellseg_path)) && (~ischar(cellseg_path))
    fprintf("Cellseg data path argument is invalid.\n");
    return;
end
if total_ch < 1
    fprintf("Input image must have at least one channel.\n");
    return;
end
if light_ch > total_ch
    fprintf("Light (TRANS) channel index is invalid.\n");
    return;
end
if rna_ch > total_ch
    fprintf("RNA/Signal channel index is invalid.\n");
    return;
end
if ~islogical(overwrite_output)
    if ischar(overwrite_output)
        overwrite_output = startsWith(overwrite_output, "true", 'IgnoreCase', true);
    else
        if isnumeric(overwrite_output)
            if overwrite_output == 0
                overwrite_output = false;
            else
                overwrite_output = true;
            end
        else
            overwrite_output = false;
        end
    end
end

%Set defaults
if ttune_winsize < 1
    ttune_winsize = 10;
end
if ttune_wscorethresh < 0
    ttune_wscorethresh = 0.9;
end
if t_min < 1
    t_min = 1;
end
if t_max < 1
    t_max = 300;
end
if isempty(out_dir)
    %Defaults to input directory
    [out_dir, ~, ~] = fileparts(tif_path);
    fprintf("Output path not provided. Set to %s\n", out_dir);
end
if isempty(img_name)
    %Defaults to input file name
    [~, img_name, ~] = fileparts(tif_path);
    fprintf("Image name not provided. Set to %s\n", img_name);
end
if ~isempty(ctrl_path)
    if ~ischar(ctrl_path)
        fprintf("Control path is not valid. Setting to empty.\n");
        ctrl_path = "";
    end
end
if ctrl_ch > ctrl_chcount
    fprintf("Control RNA/signal channel index is invalid. Control will not be used.\n");
    ctrl_path = "";
end

%Debug print
fprintf("Running RNASpots with the following parameters...\n");
fprintf("img_name = %s\n", img_name);
fprintf("tif_path = %s\n", tif_path);
fprintf("rna_ch = %d\n", rna_ch);
fprintf("light_ch = %d\n", light_ch);
fprintf("total_ch = %d\n", total_ch);
fprintf("out_dir = %s\n", out_dir);
fprintf("t_min = %d\n", t_min);
fprintf("t_max = %d\n", t_max);
fprintf("ztrim = %d\n", ztrim);
fprintf("cellseg_path = %s\n", cellseg_path);
fprintf("ctrl_path = %s\n", ctrl_path);
fprintf("ctrl_ch = %d\n", ctrl_ch);
fprintf("ctrl_chcount = %d\n", ctrl_chcount);
fprintf("ttune_winsize = %d\n", ttune_winsize);
fprintf("ttune_wscorethresh = %f\n", ttune_wscorethresh);
fprintf("overwrite_output = %d\n", overwrite_output);

%Do background extraction if arguments provided.
%!! Don't redo if target exists and overwrite output is false!
bkg_path = [out_dir filesep 'bkgmask' filesep img_name '_bkg'];
if isempty(ctrl_path)
    if ~overwrite_output && isfile(bkg_path)
        fprintf("Background mask already exists at %s! Skipping extraction...\n", bkg_path);
    else
        if (~isempty(cellseg_path)) && (light_ch > 0)
            fprintf("Extracting background...\n");
            %Make sure cellseg data file exists
            if (isfile(cellseg_path))
                Main_BackgroundMask(tif_path, cellseg_path, bkg_path, total_ch, light_ch, true);
            else
                fprintf("Cellseg file %s does not exist! Skipping background extraction...\n", cellseg_path);
            end
        end
    end
else
    fprintf("Control image provided - background extraction skipped.\n");
end

%Run spot detect
%!! Don't redo if target exists and overwrite output is false!
fprintf("Running spot detect... (This may take a few hours on large files)\n");
out_stem = Main_RNASpotDetect(img_name, tif_path, out_dir, rna_ch, total_ch, t_min, t_max, true, overwrite_output);

%Run spot detect on control (if a tif path was provided)
%   Detect if input control path is a tif. If not, assume the input is a spot
%       detect results path stem.
%!! Don't redo if target exists and overwrite output is false!
if ~isempty(ctrl_path)
    if endsWith(ctrl_path, ".tif")
        if ~isfile(ctrl_path)
            fprintf("Control path %s does not exist. Aborting...\n", ctrl_path);
            return;
        end
        
        fprintf("Running spot detect on control image... (This may take a few hours on large files)\n");
        ctrl_stem = Main_RNASpotDetect(img_name, tif_path, out_dir, rna_ch, total_ch, t_min, t_max, true, overwrite_output);
    else
        ctrl_stem = ctrl_path;
    end
end

%Filter spots for bkg (if mask exists)
bkg_filter_stem = [out_stem '_bkgmasked'];
if isfile(bkg_path)
    fprintf("Filtering detected coordinates through background mask for control...\n");
    
    %Coord & spot count tables
    load(bkg_path, 'background_mask');
    tbl_path_RNA = [out_stem '_coordTable'];
    load(tbl_path_RNA, 'coord_table');
    [masked_spot_table, masked_coord_table] = RNA_Threshold_Common.mask_spots(background_mask, coord_table);
    coord_table = masked_coord_table;
    spot_table = masked_spot_table;
    save([bkg_filter_stem '_coordTable'], 'coord_table');
    save([bkg_filter_stem '_spotTable'], 'spot_table');
    
    %Masked image views
    load([out_stem '_imgviewstructs'], 'my_images');
    my_images(1).image = immultiply(my_images(1).image, background_mask);
    my_images(2).image = immultiply(my_images(2).image, background_mask);
    save([bkg_filter_stem '_imgviewstructs'], 'my_images');
    
    if isempty(ctrl_stem)
        ctrl_stem = bkg_filter_stem;
    end
end

%Load spot count tables and apply ztrim
runinfo_path = [out_stem '_runinfo'];
if ztrim > 0
    %Find image dimensions
    if isfile(runinfo_path)
        %Load the run info.
        load(runinfo_path, 'idims');
    else
        %Reload TIF and count manually.
        [X,Y,Z] = GetTifDims(tif_path, total_ch);
        idims = struct('x', X, 'y', Y, 'z', Z);
    end
    %Create ztrim mask.
    ztrim_mask = true(idims.y, idims.x, idims.z);
    ztrim_mask(:,:,1:ztrim) = false;
    ztop = idims.z - ztrim +  1;
    ztrim_mask(:,:,ztop:idims.z) = false;
    
    %Load the tables
    load([out_stem '_coordTable'], 'coord_table');
    coords_sample = coord_table;
    
    %Mask the tables
    [spots_sample, coords_sample] = RNA_Threshold_Common.mask_spots(ztrim_mask, coords_sample);
    
    %--- And repeat with the control
    if ~isfile(bkg_path)
        %Need to retrieve dims again and make a new mask.
        ctrl_runinfo_path = [ctrl_stem '_runinfo'];
        if isfile(ctrl_runinfo_path)
            %Load the run info.
            load(ctrl_runinfo_path, 'idims');
        else
            %Reload TIF and count manually.
            [X,Y,Z] = GetTifDims(tif_path, total_ch);
            idims = struct('x', X, 'y', Y, 'z', Z);
        end
        ztrim_mask = true(idims.y, idims.x, idims.z);
        ztrim_mask(:,:,1:ztrim) = false;
        ztop = idims.z - ztrim +  1;
        ztrim_mask(:,:,ztop:idims.z) = false;
    end
    load([ctrl_stem '_coordTable'], 'coord_table');
    coords_control = coord_table;
    [spots_control, coords_control] = RNA_Threshold_Common.mask_spots(ztrim_mask, coords_control);
else
    %Just load the coord and spot tables
    load([out_stem '_coordTable'], 'coord_table');
    load([out_stem 'spot_table'], 'spot_table');
    coords_sample = coord_table;
    spots_sample = spot_table;
    load([ctrl_stem '_coordTable'], 'coord_table');
    load([ctrl_stem 'spot_table'], 'spot_table');
    coords_control = coord_table;
    spots_control = spot_table;
end

%Detect threshold
[thresh, win_stdevs] = RNA_Threshold_Common.estimateThreshold(spots_sample, spots_control, ttune_winsize, 0.0, ttune_wscorethresh);

%Print plots & image representation
    %Spot plots (log and linear scale) - w/ chosen threshold marked
    %Window score plot
    %Max projection of sample and control w circled spots (also max proj)

%Save THIS module's run info (including paths, parameters, chosen threshold etc)
    %Save ztrimmed coords/spotcounts here too. Remove coord tables below
    %   threshold 10 to save space.

end