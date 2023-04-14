%GUI module for interactive manual curation of RNA spot detection results.
%Blythe Hospelhorn
%Version 1.10.0
%Updated April 13, 2023

%Update Log:
%   1.0.0 | 21.03.12
%       Init doc
%   1.1.0 | 21.12.14
%       Large images getting too big to save to one file
%       Moved saves for pos table and neg table to their own files
%   1.2.0 | 22.03.25
%       Version 4 save (save toggles)
%       Force ptbl and ntbl to uint16 (TOO BIG!!!)
%       Fixed load + change savestem bug
%   1.3.0 | 22.07.01
%       Version 7 save - Ref coords to separate file
%   1.4.0 | 22.08.26
%       "ztrim" deprecated. Instead specify min and max Z.
%       Can now generate a selector without coordinate table (for agnostic
%           truthset creation)
%   1.5.0 | 22.09.01
%       Added some functions for quickly handling fscores
%   1.6.0 | 22.09.09
%       Added ability to change color of crosshair
%       Debugging some z masking and max proj rendering issues
%   1.6.1 | 22.09.28
%       Added timestamp to save file
%   1.6.2 | 22.10.26
%       Fixed bug with initializing zmin/zmax of new annoobj
%   1.7.0 | 22.12.15
%       Added ability to swap between views of raw or filtered images in
%       ref mod (max proj only - need to rearrange some stuff to be able to
%       look at full 3D raw image...)
%   1.8.0 | 23.02.24
%       Misc bug fixes, added a couple of functions
%   1.9.0 | 23.03.07
%       Overhauled snap function to make it faster.
%       Also allowed it to take parameters specifying search 3d and z radii
%       3D radius defaults to 4, z radius defaults to 2 (reduced from
%       before - also before didn't limit on z alone)
%   1.9.1 | 23.03.08
%       (Begrudgingly) updated to save ptbl to v7.3 so can store enormous
%       tables.
%   1.10.0 | 23.04.13
%       Added genScoreResponseTable


%%
classdef RNA_Threshold_SpotSelector
    
    %TODO - There is a bug with the z trim command, so can't change z trim
    %from GUI. 
    %TODO: Clean out 0,0,0 points that randomly appear. I think rendering
    %stops after encountering these?
    %TODO: Some spots apparently not rendering in ref select? Is this
    %related to ztrim or 0,0,0 problem?
    
    %GUI Controls:
            %   space (or any key) - Ready listener
            %   s - Save
            %   x - Exit
            %   c - Clear
            %   e - Clear all but current threshold
            %   a - Auto-reset
            %   r - Check against ref truth set
            %   p - Print counts to console
            %   P - Print ALL counts (all thresholds) to console
            %   4 - Set sample mask (1/4 image size)
            %   6 - Set sample mask (1/6 image size)
            %   8 - Set sample mask (1/8 image size)
            %   1 - Remove sampler mask
            %   S - (Ref mode) snap to computer calls
            %   U - (Ref mode) toggle remove unsnapped spots from ref set
            %   Left Click - Set (or unset) false positive
            %   Right Click - Set (or unset) false negative
            %   Enter - Input done (return to defo fig controls)
            %   > - Increase threshold by 1
            %   < - Decrease threshold by 1
            %   ] - Increase threshold by 10
            %   [ - Decrease threshold by 10
            %   Z - Toggle Z min/max control
            %   ^ - Increase Z trim
            %   v - Decrease Z trim
            
            %   + - Up one slice
            %   - - Down one slice
            %   = - Up 10 slices
            %   _ - Down 10 slices

            %   f - Switch between filtered and raw images. Right now, only
            %       works for max proj.
            %   m - Toggle between max-proj/single z slice (if 3D)
            %   z - Toggle on/off all-z selection
            %   # - Toggle on/off z spot counts
            %   % - Toggle z-color scale (all z, or just nearest slices?)
            %   C - Toggle contrast scale (max proj versus indiv slice -
            %       only available for slice view)
    
    %%
    properties
        
        save_stem; %Path stem for saving figures as images
        img_structs; %Hold images and intensity scaling factors
        imgdat_path; %Path to mat file containing the filtered channel data (for slice view)
        
        %tif_path; %Path to source tif (only used if loading orig. is requested by user)
        %tif_chcount; %Total channels in source tif (only used if loading orig. is requested by user)
        %tif_ch; %Channel of interest in source tif (only used if loading orig. is requested by user)
        loaded_ch; %Loaded full source channel data
        slice_drawn; %I think this is a pseudo-boolean for keeping track of whether the specified slice has been drawn yet
        current_slice; %Current z-slice in view (either for orig image or 3D view) 
        max_slice; %Max z-slice or total number of z slices
        crosshair_color;
        
        ztrim; %[DEPR] Number of z planes from top and bottom to mask out in spot counting and selection.
        z_min; %Index of lowest z plane to include
        z_max; %Index of highest z plane to include.
        selmcoords; %Optional mask def for spots and selection area (because some images just have too damn many spots)
        
        mode_3d; %true if in 3D  mode (determined from coord table)
        threshold_idx; %Set threshold for spot detection lookup
       	threshold_table; %Table of actual threshold values used (index may or may not be actual value)
        
        positives; %cell(t_count) of int[n][4] - x,y,z,t/f/e (Table of positives, whether true, false, or excluded)
        false_negs; %cell(t_count) of int[n][3] - x,y,z (Table of false negatives)
        ref_coords; %int[n][3] - x,y,z (Table of manually selected ref spots)
        f_scores; %double[t][4] (fscore[t][1], truepos[t][2], falsepos[t][3], falseneg[t][4] at each thresh value: recalculated whenever snap to ref coords is performed)
        f_scores_dirty; %Flag for when fscores need to be updated.
        ref_last_modified; %Datetime when reference set was last modified
        
        fh_filter; %Figure handle for filtered image max z projection
        fh_raw; %Figure handle for max z proj of (contrast enhanced) raw image
        fh_origslice; %
        fh_3d; % a 3D plot of points
        
        alloc_size;
        n_table;
        n_alloc;
        n_used;
        
        toggle_singleSlice; %(3D mode) whether to render only selected z slice or max projection
        toggle_allz;
        toggle_3dcount;
        toggle_clr_local;
        toggle_del_unsnapped;
        toggle_cscale_max; % (3D mode) whether contrast scale for indiv z slices uses the full image scale or a scale for that slice
        toggle_change_zminmax; %Not saved. Toggles whether ^v update the min (false) or max (true) Z
        toggle_raw_view = false; %If true, view raw image in ref mode. If false, filtered image (default).
        
        loop_breaker;
        
        colortbl_red;
        colortbl_magenta;
        colortbl_yellow;
        colortbl_cyan;
        colortbl_white;
        
    end
    
    %%
    methods
        
        %%
        % Initialize a new spot selector instance using the provided
        % threshold table and the spot detection data found at the provided
        % save path.
        %
        %ARGS
        %   obj - RNA_Threshold_SpotSelector to perform method on (applies
        %       to all instance methods)
        %   save_stem (String) - Path prefix to previously generated spot
        %       detection save data
        %   init_thresh_idx (int) - Index of initial threshold value to set
        %       the spot selector to
        %   th_table (int[t_count][1]) - Table of threshold values
        %       (OPTIONAL)
        %
        function obj = initializeNew(obj, save_stem, init_thresh_idx, th_table)
            
            %First, look for run specs...
            spotsrun = [];
            if isfile([save_stem '_rnaspotsrun.mat'])
                spotsrun = RNASpotsRun.loadFrom(save_stem);
                spotsrun.out_stem = save_stem;
                t_count = spotsrun.t_max - spotsrun.t_min + 1;
                [~,coord_table] = spotsrun.loadCoordinateTable();
            else
                t_count = size(th_table,1); 
               %Load coordinate table
                coord_table_suffix = '_coordTable';
                tbl_path = [save_stem coord_table_suffix];
                load(tbl_path, 'coord_table');
            end

            obj.positives = cell(t_count, 1);
            ct_sub = coord_table{1,1};
            dimcount = size(ct_sub,2);
            for t = 1:t_count
                src_tbl = coord_table{t};
                tbl = ones([size(src_tbl,1) 4]);
                if dimcount == 2
                    tbl(:,1:2) = src_tbl(:,1:2);
                    obj.mode_3d = false;
                else
                    tbl(:,1:3) = src_tbl(:,1:3);
                    obj.mode_3d = true;
                end
                tbl = uint16(tbl);
                obj.positives{t} = tbl;
            end
            obj.false_negs = cell(t_count, 1);
            
            %Load image structure
            imgstructs_suffix = '_imgviewstructs';
            istruct_path = [save_stem imgstructs_suffix '.mat'];

            Z = 1;
            if isfile(istruct_path)
                load(istruct_path, 'my_images');
                obj.img_structs = my_images;
                s_img = my_images(1).image;
                Y = size(s_img, 1);
                X = size(s_img, 2);
            else
                %Estimate size from coords and generate dummy
                fprintf('WARNING: Image struct file could not be found. GUI will be unavailable for selector.\n');
                th_coords = coord_table{1,1};
                X = max(th_coords(:,1),[],'all');
                Y = max(th_coords(:,2),[],'all');
                Z = max(th_coords(:,3),[],'all');
                dummy_img = NaN(Y,X);
                my_images(2) = struct('image', dummy_img, 'Lmin', 0, 'Lmax', 0);
                my_images(1) = struct('image', dummy_img, 'Lmin', 0, 'Lmax', 0);
                obj.img_structs = my_images;
                save(istruct_path, 'my_images');
            end
            
            %Load filtered image channel
            obj.imgdat_path = [save_stem '_prefilteredIMG.mat'];
            if isfile(obj.imgdat_path)
                load(obj.imgdat_path, 'img_filter');
                obj.loaded_ch = double(img_filter);
            else
                obj.loaded_ch = NaN(Y,X,Z);
            end
            
            %Set defaults
            obj.threshold_table = zeros([t_count, 2]);
            %obj.threshold_table(:,1) = th_table(:,:);
            obj.save_stem = save_stem;
            obj.threshold_idx = init_thresh_idx;
            obj.toggle_del_unsnapped = false;
            obj.toggle_cscale_max = false;
            %obj.ztrim = 0;
            obj.selmcoords = [];
            
            obj.alloc_size = X*Y;
            obj.n_alloc = 0;
            obj.n_used = 0;
            
            %obj.tif_path = [];
            obj.slice_drawn = 0;
            obj.current_slice = 1;
            obj.max_slice = 1;
            
            %TODO update z trim
            if ~isempty(spotsrun)
                obj.ztrim = spotsrun.ztrim;
                spotsrun = spotsrun.updateZTrimParams;
                obj.threshold_table(:,1) = transpose(spotsrun.t_min:1:spotsrun.t_max);
                if obj.mode_3d
                    %Find the max z
                    obj.max_slice = spotsrun.idims_sample.z;
                    obj.current_slice = uint16(obj.max_slice./2);
                end
                obj.z_min = spotsrun.z_min_apply;
                obj.z_max = spotsrun.z_max_apply;
            else
                obj.ztrim = 0;
                obj.threshold_table(:,1) = th_table(:,:);
                if obj.mode_3d
                    %Find the max z
                    obj = obj.setMaxZ();
                    obj.current_slice = uint16(obj.max_slice./2);
                end
                obj.z_min = 1;
                obj.z_max = obj.max_slice;
            end
            
            obj.f_scores = NaN(t_count,4);
            
            obj.toggle_singleSlice = false;
            obj.toggle_allz = ~obj.toggle_singleSlice;
            obj.toggle_3dcount = false;
            obj.toggle_clr_local = false;
            obj.toggle_del_unsnapped = false;
            obj.f_scores_dirty = false;
            obj.toggle_raw_view = false;
            
            obj.toggle_change_zminmax = false;
            obj.crosshair_color = [0.000, 0.000, 0.000];
            obj.ref_last_modified = datetime;
        end
        
        %%
        % Launch the spot selector GUI as specified by this spot selector
        % instance. This method launches the main version of the GUI that
        % allows for curation of the automatically detected spot set.
        %
        function obj = launchGUI(obj)
            
            %Block if this is a ref-only selector.
            if isempty(obj.positives)
                fprintf("Selector is for agnostic reference set generation only. GUI will not be launched.\n");
                return;
            end
            
            %Set some flags
            obj.toggle_allz = ~obj.toggle_singleSlice;
            %obj.toggle_3dcount = false;
            %obj.toggle_clr_local = false;
            %obj.toggle_singleSlice = false;
            %obj.toggle_del_unsnapped = false;
            
            %Generate color tables
            obj = obj.generateColorTables();
            
            %Set some stuff...
            obj.loop_breaker = 0;
            obj.threshold_table(obj.threshold_idx, 2) = 1;
            
            %Draw initial images
            obj = obj.drawImages();
            
            %Interaction Loop
            while obj.loop_breaker == 0
                w = waitforbuttonpress;
                if w
                    obj = obj.onReadyKey();
                end
            end
            
        end
        
        %%
        % Launch the GUI for manual reference set specification from this
        % spot selector instance. This GUI version does not display auto
        % detection data, instead allowing the user to mark spots on the
        % image.
        %
        function obj = launchRefSelectGUI(obj)
            
            %Set some flags
            obj.toggle_allz = ~obj.toggle_singleSlice;
            %obj.toggle_3dcount = false;
            %obj.toggle_clr_local = false;
            %obj.toggle_singleSlice = true;
            %obj.toggle_del_unsnapped = false;
            
            %Generate color tables
            obj = obj.generateColorTables();
            
            %Set some stuff...
            obj.loop_breaker = 0;
            
            %Draw initial images
            obj = obj.drawRefImage();
            
            %Interaction Loop
            while obj.loop_breaker == 0
                w = waitforbuttonpress;
                if w
                    obj = obj.onReadyKey_RefMode();
                end
            end
            
        end
        
        %%
        % Save the state of this spot selector instance including all ref,
        % pos, and neg tables. This allows for later reloading of curated
        % spot set, even including between ref GUI and regular GUI.
        %
        function obj = saveMe(obj)
            
            save_path = [obj.save_stem 'spotAnnoObj'];
            %save(save_path, 'spot_annotator');
            fprintf("Saving to %s...\n", save_path)
            
            istructs = obj.img_structs;
            th_idx = obj.threshold_idx;
            th_tbl = obj.threshold_table;
            pos_tbl = obj.positives;
            neg_tbl = obj.false_negs;
            
            %tiff_path = obj.tif_path;
            %tiff_channels = obj.tif_chcount;
            %tiff_ch_selected = obj.tif_ch;
            ref_coord_tbl = obj.ref_coords;
            bool3d = obj.mode_3d;
            lastz = obj.current_slice;
            maxz = obj.max_slice;
            filimg_path = obj.imgdat_path;
            z_trim = obj.ztrim;
            mask_selection = obj.selmcoords;
            save_ver = 13;
            
            %Version 4+
            toggle_ss = obj.toggle_singleSlice;
            toggle_az = obj.toggle_allz;
            toggle_3dc = obj.toggle_3dcount;
            toggle_cl = obj.toggle_clr_local;
            toggle_du = obj.toggle_del_unsnapped;
            
            %Version 5+
            toggle_cs = obj.toggle_cscale_max;
            
            %Version 6+
            ftable = obj.f_scores;
            
            %Version 8+
            refonly = isempty(obj.positives);
            zmin = obj.z_min;
            zmax = obj.z_max;
            
            %Version 9+
            flag_fscores_dirty = obj.f_scores_dirty;
            
            %Version 10+
            crossclr = obj.crosshair_color;

            %Version 11+
            timestamp = datetime;
            rtimestamp = obj.ref_last_modified;
            
            %Version 12+
            toggle_rvr = obj.toggle_raw_view;
            
            %save(save_path, 'istructs', 'th_idx', 'th_tbl', 'pos_tbl', 'neg_tbl', 'tiff_path', 'tiff_channels', 'tiff_ch_selected', 'ref_coord_tbl', 'bool3d', 'lastz', 'maxz');
            %save(save_path, 'istructs', 'th_idx', 'th_tbl', 'pos_tbl', 'neg_tbl', 'ref_coord_tbl', 'bool3d', 'lastz', 'maxz', 'filimg_path', 'z_trim', 'mask_selection', 'save_ver');
            save(save_path, 'istructs', 'th_idx', 'th_tbl', 'bool3d', 'lastz', 'maxz', 'filimg_path', 'z_trim', 'mask_selection', 'save_ver', 'toggle_ss', 'toggle_az' ,'toggle_3dc','toggle_cl','toggle_du','toggle_cs','ftable','refonly','zmin','zmax','flag_fscores_dirty','crossclr','timestamp','toggle_rvr');
            save([save_path '_ptbl'], 'pos_tbl', '-v7.3'); %Ver 13+, upped to v7.3. This will be very bad for disk space but oh well.
            save([save_path '_ntbl'], 'neg_tbl');
            save([save_path '_refset'], 'ref_coord_tbl', 'rtimestamp'); %Ver 7+
            fprintf("Save complete!\n");
        end
        
        %%
        % (Internal method)
        % Store a temporary false negative call to this Selector. All
        % temp false negatives are saved to false neg table when selection
        % mode is ended in GUI. (z coord will be object's current slice
        % coordinate)
        %
        % ARGS
        %   x (double) - appr. x coordinate of spot to add (rounded)
        %   y (double) - appr. y coordinate of spot to add (rounded)
        %
        function obj = addNeg(obj, x, y)
            %fprintf("n_alloc: %d, n_used: %d\n", obj.n_alloc, obj.n_used)
            
            if obj.n_alloc == obj.n_used
                %Reallocate
                %fprintf("Reallocation required.\n")
                if obj.n_alloc == 0
                    %fprintf("New realloc\n")
                    obj.n_table = zeros([256 4]);
                    obj.n_table(:,4) = 1;
                    obj.n_alloc = 256;
                else
                    now_sz = size(obj.n_table,1);
                    %Increase by 256
                    new_size = now_sz + 256;
                    %fprintf("Increase size to %d...\n", new_size)
                    new_table = zeros([new_size 4]);
                    obj.n_table(:,4) = 1;
                    new_table(1:now_sz,:) = obj.n_table(:,:);
                    obj.n_table = new_table;
                    obj.n_alloc = now_sz + 256;
                end
            end
            
            add_idx = obj.n_used + 1;
            %fprintf("Adding at index %d\n", add_idx);
            obj.n_table(add_idx, 1) = round(x);
            obj.n_table(add_idx, 2) = round(y);
            obj.n_table(add_idx, 3) = obj.current_slice;
            obj.n_table(add_idx, 4) = 0;
            obj.n_used = add_idx+1;
            
        end
        
        %%
        % (Internal method)
        % Copy the entries in the temporary false negative table to the
        % object's main storage false negative table and clear the temp
        % table.
        %
        function obj = flushNegTable(obj)
            
            if isempty(obj.n_table)
                return;
            end
            
            [~, fneg_idx] = RNA_Threshold_SpotSelector.isolateTFIndices(obj.n_table);
            fneg_tbl = obj.n_table(fneg_idx,:);
            good_idx = RNA_Threshold_SpotSelector.findNonzeroCoords(fneg_tbl, obj.mode_3d);
            fneg_tbl = fneg_tbl(good_idx,:);
            obj.false_negs{obj.threshold_idx} = fneg_tbl;
            
            obj.n_table = [];
            obj.n_alloc = 0;
            obj.n_used = 0;
            
        end
        
        %%
        % (Internal method)
        % Flush the temporary coordinate table to the ref coord table (copy
        % coordinates, clear temp table)
        %
        function obj = flushRefTable(obj)
            
            if isempty(obj.n_table)
                return;
            end
            
            [~, fneg_idx] = RNA_Threshold_SpotSelector.isolateTFIndices(obj.n_table);
            obj.ref_coords = obj.n_table(fneg_idx,1:3);
            good_idx = RNA_Threshold_SpotSelector.findNonzeroCoords(obj.ref_coords, obj.mode_3d);
            obj.ref_coords = obj.ref_coords(good_idx,:);
            
            obj.n_table = [];
            obj.n_alloc = 0;
            obj.n_used = 0;
            
        end
        
        %%
        % (Internal method)
        % Break the coordinate table for the given threshold into
        % sub-tables for true positives, false positives, false negatives,
        % and excluded points to make it easier to draw and count them.
        %
        % ARGS
        %   th_idx (int) - Index of threshold to pull coord table from
        %   filter_bool (bool) - Whether to apply selection mask
        %
        function [obj, tpos, fpos, fneg, epos] = splitCoordTables(obj, th_idx, filter_bool)
            
            %Debug lines (it's not drawing all circles on some slices for
            %some reason???)
            %fprintf("[DEBUG] splitCoordTables -- INPUT:\n");
            %fprintf("\tThreshold Index = %d\n", th_idx);
            %fprintf("\tz = %d\n", obj.current_slice);
            
            if nargin < 3
                filter_bool = false;
            end
            
            %Get mask
            smask = [];
            if filter_bool
                if obj.ztrim > 0 | obj.z_min > 1 | obj.z_max < obj.max_slice | ~isempty(obj.selmcoords)
                    [obj, smask] = obj.genCountMask();
                end
            end
            
            if isempty(obj.false_negs{th_idx})
                fneg = [];
            else
                fneg = obj.false_negs{th_idx}(:,1:3);
                %fprintf("\tFNeg Count: %d\n", size(fneg,1));
                if ~isempty(smask)
                    fneg = RNA_Threshold_SpotSelector.maskFilterSpotTable(fneg, smask);
                end
            end
            
            pos_tbl = obj.positives{th_idx};
            [tpos_idx, fpos_idx, epos_idx] = RNA_Threshold_SpotSelector.isolateTFEIndices(pos_tbl);
            %fprintf("\tAll Pos Count: %d\n", size(pos_tbl,1));
            
            if isempty(tpos_idx)
                tpos = [];
            else
                tpos = pos_tbl(tpos_idx,1:3);
                %fprintf("\tTPos Count: %d\n", size(tpos,1));
                if ~isempty(smask)
                    tpos = RNA_Threshold_SpotSelector.maskFilterSpotTable(tpos, smask);
                    %fprintf("DEBUG smask not empty\n");
                end
            end
            
            if isempty(fpos_idx)
                fpos = [];
            else
                fpos = pos_tbl(fpos_idx,1:3);
                %fprintf("\tFPos Count: %d\n", size(fpos,1));
                if ~isempty(smask)
                    fpos = RNA_Threshold_SpotSelector.maskFilterSpotTable(fpos, smask);
                end
            end
            
            if isempty(epos_idx)
                epos = [];
            else
                epos = pos_tbl(epos_idx,1:3);
                %fprintf("\tEPos Count: %d\n", size(epos,1));
                if ~isempty(smask)
                    epos = RNA_Threshold_SpotSelector.maskFilterSpotTable(epos, smask);
                end
            end
            
            %fprintf("[DEBUG] splitCoordTables -- OUTPUT:\n");
            %fprintf("\tFNeg Count: %d\n", size(fneg,1));
            %fprintf("\tTPos Count: %d\n", size(tpos,1));
            %fprintf("\tFPos Count: %d\n", size(fpos,1));
            %fprintf("\tEPos Count: %d\n", size(epos,1));
            
        end
        
        %%
        % (Internal method)
        % Draw a single marker at x,y on the current raw & filtered image
        % figures.
        %
        % ARGS
        %   x (int) - x coordinate of point to draw
        %   y (int) - y coordiante of point to draw
        %   marker_str (string) - string specifying string marker type (eg.
        %       'or' for a red circle)
        %
        function obj = drawSingleCircle(obj, x, y, marker_str)
            
            obj.fh_raw;
            hold on;
            plot(x, y, marker_str,'markersize',10); 
            
            obj.fh_filter;
            hold on;
            plot(x, y, marker_str,'markersize',10); 
            
        end
        
        %%
        % (Internal method)
        % Set the current axes to axes scaled properly for the 3D plot.
        %
        function obj = get3DPlotAxes(obj)
            img = obj.img_structs(2);
            X = size(img.image,2);
            Y = size(img.image,1);
            Z = obj.max_slice;
            
            %ax = gca;
            ax = axes();
            ax.XLim = [1, X];
            ax.YLim = [1, Y];
            ax.ZLim = [1, Z];
            ax.XAxisLocation = 'top';
            ax.YDir = 'reverse';
            ax.XLabel.String = 'X Axis';
            ax.YLabel.String = 'Y Axis';
            ax.ZLabel.String = 'Z Axis';
            
        end
        
        %%
        % (Internal method)
        % Plot the provided table of xyz coordinates to the instance's
        % currently loaded 3D plot using the provided color.
        %
        % ARGS
        %   tbl (int[n][3+]) - List of coordinates of points to plot 
        %   base_color (double[1][3]) - [r,g,b] (0.0 - 1.0) of base color
        %       to use for points
        %
        function obj = draw3DPlot(obj, tbl, base_color)  
            figure(obj.fh_3d);
            hold on;
            %plot3(tbl(:,1), tbl(:,2), tbl(:,3), 'o', 'MarkerEdgeColor', e_color, 'MarkerFaceColor', f_color);  
            
            tbl = double(tbl);
            img = obj.img_structs(2);
            X = double(size(img.image,2));
            Y = double(size(img.image,1));
            Z = double(obj.max_slice);
            max_rad = sqrt(X^2 + Y^2 + Z^2);
            
            count = size(tbl, 1);
            for i = 1:count
                x = tbl(i,1);
                y = tbl(i,2);
                z = tbl(i,3);
                rad = sqrt(x^2 + y^2 + z^2);
                rfrac = rad/max_rad;
                
                r = base_color(1,1) * rfrac;
                g = base_color(1,2) * rfrac;
                b = base_color(1,3) * rfrac;
                
                plot3(x, y, z, 'o','MarkerEdgeColor', base_color, 'MarkerFaceColor', [r,g,b]);  
            end
            
            grid on;
            
        end
        
        %%
        % (Internal method)
        % Plot the provided table of xyz coordinates to the instance's
        % currently loaded 2D raw and filtered projections. Z coordinates
        % are either mapped by color (further from current plane = darker
        % color) in max proj mode.
        %
        % ARGS
        %   tbl (int[n][2+]) - List of coordinates of points to plot 
        %   color_tbl (double[n][3]) - [r,g,b] (0.0 - 1.0) table of marker
        %       colors to use for gradient - higher indices for points
        %       further away in z
        %
        function obj = drawMulti(obj, tbl, color_tbl)
            
            %count = size(all_x,1);
            %colors
            
            %separate by z diff...
            Z = obj.max_slice+1;
            cz = double(obj.current_slice);
            z_groups = cell(Z,1);
            diffs = abs(tbl(:,3) - cz);
            if obj.toggle_clr_local
                incr = Z./5;
                difflvl = (diffs*incr);
                diffs = min(difflvl, obj.max_slice);
            end
            
            for i = 1:Z
                val = i-1;
                [rows, ~] = find(diffs == val);
                z_groups{i} = tbl(rows, :);
            end
            
            figure(obj.fh_raw);
            hold on;
            for i = Z:-1:1
                %plot(all_x(i), all_y(i),'Marker','o','MarkerEdgeColor',colors(i,:),'markersize',10); 
                ztbl = z_groups{i};
                if ~isempty(ztbl)
                    plot(ztbl(:,1), ztbl(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',color_tbl(i,:),'markersize',10);
                end 
            end
            
            figure(obj.fh_filter);
            hold on;
            if (obj.toggle_singleSlice)
                %Just the spots on this z plane
                [rows, ~] = find(tbl(:,3) == obj.current_slice);
                subtbl = tbl(rows,:);
                if ~isempty(subtbl)
                    plot(subtbl(:,1), subtbl(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',color_tbl(1,:),'markersize',10);
                end
            else
                for i = Z:-1:1
                    %plot(all_x(i), all_y(i),'Marker','o','MarkerEdgeColor',colors(i,:),'markersize',10); 
                    ztbl = z_groups{i};
                
                    if ~isempty(ztbl)
                        plot(ztbl(:,1), ztbl(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',color_tbl(i,:),'markersize',10);
                    end 
                end
            end
            
        end
        
        %%
        % (Internal method)
        % Draw/redraw the figures for the main spot selector GUI based upon
        % the current state.
        %
        function obj = drawImages(obj)
            
            th_idx = obj.threshold_idx;
            th = obj.threshold_table(th_idx, 1);
            
            [obj, t_pos, f_pos, f_neg, e_pos] = obj.splitCoordTables(th_idx, true);
            
            tp_count = RNA_Threshold_SpotSelector.countFromTable(obj, t_pos);
            fp_count = RNA_Threshold_SpotSelector.countFromTable(obj, f_pos);
            fn_count = RNA_Threshold_SpotSelector.countFromTable(obj, f_neg);
            
            Y = size(obj.img_structs(1).image,1);
            X = size(obj.img_structs(1).image,2);
            immul_mask = double(ones(Y,X));
            dval = 0.25;
            if ~isempty(obj.selmcoords)
                %Darken pixels outside the mask.
                immul_mask(:,1:obj.selmcoords(1,1)) = dval;
                immul_mask(:,obj.selmcoords(2,1):X) = dval;
                immul_mask(1:obj.selmcoords(3,1),:) = dval;
                immul_mask(obj.selmcoords(4,1):Y,:) = dval;
            end
            
            %Raw - Update only draw those where z is also current_slice
            obj.fh_raw = figure(11);
            clf;
            img2 = obj.img_structs(2);
            imshow(immultiply(img2.image, immul_mask),[img2.Lmin img2.Lmax]);
            hold on;
            title(['Threshold = ' num2str(th) ' | t+: ' num2str(tp_count) ', f+: ' num2str(fp_count) ', f-: ' num2str(fn_count)]);
            impixelinfo;
            
            %Filtered
            obj.fh_filter = figure(10);
            clf;
            if (obj.toggle_singleSlice)
                slice = obj.loaded_ch(:,:,obj.current_slice);
                Lmin = min(slice(:));
                Lmax = median(slice(:)) + round(10 * std(slice(:)));
                
                %See if slice is within mask, if mask is present.
                if (obj.current_slice < obj.z_min) | (obj.current_slice >obj.z_max)
                    slice = slice * dval;
                else
                    slice = immultiply(slice, immul_mask);
                end
                
                imshow(slice,[Lmin Lmax]);
            else
                img1 = obj.img_structs(1);
                imshow(immultiply(img1.image, immul_mask),[img1.Lmin img1.Lmax]);
            end
            hold on;
            if ~isempty(f_pos)
                %plot(f_pos(:,1), f_pos(:,2),'oy','markersize',10); 
                %plot(f_pos(:,1), f_pos(:,2),'marker',o,'MarkerEdgeColor',fp_clrmtx,'markersize',10); 
                %obj.drawMulti(f_pos(:,1), f_pos(:,2), fp_clrmtx);
                obj = obj.drawMulti(f_pos, obj.colortbl_yellow);
            end
            if ~isempty(e_pos)
                %plot(e_pos(:,1), e_pos(:,2),'om','markersize',10); 
                %plot(e_pos(:,1), e_pos(:,2),'marker',o,'MarkerEdgeColor',ep_clrmtx,'markersize',10); 
                %obj.drawMulti(e_pos(:,1), e_pos(:,2), ep_clrmtx);
                obj = obj.drawMulti(e_pos, obj.colortbl_magenta);
            end
            if ~isempty(f_neg)
                %plot(f_neg(:,1), f_neg(:,2),'oc','markersize',10); 
                %plot(f_neg(:,1), f_neg(:,2),'marker',o,'MarkerEdgeColor',fn_clrmtx,'markersize',10); 
                %obj.drawMulti(f_neg(:,1), f_neg(:,2), fn_clrmtx);
                obj = obj.drawMulti(f_neg, obj.colortbl_cyan);
            end
            if ~isempty(t_pos)
                %plot(t_pos(:,1), t_pos(:,2),'or','markersize',10); 
                %plot(t_pos(:,1), t_pos(:,2),'marker',o,'MarkerEdgeColor',tp_clrmtx,'markersize',10); 
                %obj.drawMulti(t_pos(:,1), t_pos(:,2), tp_clrmtx);
                obj = obj.drawMulti(t_pos, obj.colortbl_red);
            end
            
            %impixelinfo;
            if obj.mode_3d
                title(['z = ' num2str(obj.current_slice) ' |  Threshold = ' num2str(th) ' | t+: ' num2str(tp_count) ', f+: ' num2str(fp_count) ', f-: ' num2str(fn_count)]);
                
                %Do 3D plot
                %obj = obj.get3DPlotAxes();
                stdax = gca;
                obj.fh_3d = figure(3);
                clf;
                obj.fh_3d = figure(3);
                obj = obj.get3DPlotAxes();
                %hold on;
                if ~isempty(f_pos)
                    %obj.draw3DPlot(f_pos, [1.0, 1.0, 0.0], [0.6, 0.6, 0.0]);
                    obj.draw3DPlot(f_pos, [1.0, 1.0, 0.0]);
                end
                if ~isempty(e_pos)
                    %obj.draw3DPlot(e_pos, [1.0, 0.0, 1.0], [0.6, 0.0, 0.6]);
                    obj.draw3DPlot(e_pos, [1.0, 0.0, 1.0]);
                end
                if ~isempty(f_neg)
                    %obj.draw3DPlot(f_neg, [0.0, 1.0, 1.0], [0.0, 0.6, 0.6]);
                    obj.draw3DPlot(f_neg, [0.0, 1.0, 1.0]);
                end
                if ~isempty(t_pos)
                    %obj.draw3DPlot(t_pos, [1.0, 0.0, 0.0], [0.6, 0.0, 0.0]);
                    obj.draw3DPlot(t_pos, [1.0, 0.0, 0.0]);
                end
                
                axes(stdax);
                figure(obj.fh_filter);
                
            else
                title(['Threshold = ' num2str(th) ' | t+: ' num2str(tp_count) ', f+: ' num2str(fp_count) ', f-: ' num2str(fn_count)]);
            end
            impixelinfo;
            
        end
        
        %%
        %TODO
        function obj = drawOnlySlice(obj, f_pos, e_pos, f_neg, t_pos)
            if ~isempty(obj.loaded_ch)
                obj.slice_drawn = 1;
                obj.fh_origslice = figure(12);
                clf;
                slice = obj.loaded_ch(:,:,obj.current_slice);
                
                if (obj.toggle_cscale_max)
                    img1 = obj.img_structs(1);
                    imshow(slice,[img1.Lmin img1.Lmax]);
                else
                    Lmin = min(slice(:));
                    Lmax = median(slice(:)) + round(10 * std(slice(:)));
                    imshow(slice,[Lmin Lmax]);
                end
                
                hold on;
                if ~isempty(f_pos)
                    plot(f_pos(:,1), f_pos(:,2),'oy','markersize',10); 
                end
                if ~isempty(e_pos)
                    plot(e_pos(:,1), e_pos(:,2),'om','markersize',10); 
                end
                if ~isempty(f_neg)
                    plot(f_neg(:,1), f_neg(:,2),'oc','markersize',10); 
                end
                if ~isempty(t_pos)
                    plot(t_pos(:,1), t_pos(:,2),'or','markersize',10); 
                end  
                title(['Source Image Z Slice: ' num2str(obj.current_slice)]);
                impixelinfo;
            end
        end
        
        %%
        function obj = drawRefImage(obj)
            
            green1 = [0.0, 1.0, 0.0];
            green2 = [0.0, 0.50, 0.20];
            yellow1 = [1.0, 1.0, 0.0];
            yellow2 = [0.60, 0.50, 0.0];
            
            filtered_coords = obj.ref_coords;
            if (obj.ztrim > 0) | (obj.z_min > 1) | (obj.z_max < obj.max_slice) | (~isempty(obj.selmcoords))
                [obj,smask] = obj.genCountMask();
                filtered_coords = RNA_Threshold_SpotSelector.maskFilterSpotTable(filtered_coords,smask);
            end
            
            Y = size(obj.img_structs(1).image,1);
            X = size(obj.img_structs(1).image,2);
            immul_mask = double(ones(Y,X));
            dval = 0.25;
            if ~isempty(obj.selmcoords)
                %Darken pixels outside the mask.
                immul_mask(:,1:obj.selmcoords(1,1)) = dval;
                immul_mask(:,obj.selmcoords(2,1):X) = dval;
                immul_mask(1:obj.selmcoords(3,1),:) = dval;
                immul_mask(obj.selmcoords(4,1):Y,:) = dval;
            end
            
            if isempty(filtered_coords)
                p_count = 0;
            else
                p_count = size(filtered_coords,1);
            end
            
            obj.fh_filter = figure(888);
            clf;
            if (obj.toggle_singleSlice)
                slice = obj.loaded_ch(:,:,obj.current_slice);
                
                if (obj.current_slice < obj.z_min) | (obj.current_slice > obj.z_max)
                    slice = slice * dval;
                else
                    slice = immultiply(slice, immul_mask);
                end
                
                if (obj.toggle_cscale_max)
                    img1 = obj.img_structs(1);
                    imshow(slice,[img1.Lmin img1.Lmax]);
                else
                    Lmin = min(slice(:));
                    Lmax = median(slice(:)) + round(10 * std(slice(:)));
                    imshow(slice,[Lmin Lmax]);
                end
                
                ref_coords_z = RNA_Threshold_SpotSelector.filterByZ(filtered_coords, obj.current_slice);
                
                %Draw the spots on other z planes in green (above) or
                %yellow (below) - lighter if within 5 planes
                if ~isempty(filtered_coords)
                    low5 = obj.current_slice - 5;
                    high5 = obj.current_slice + 5;
                    if low5 < 1
                        low5 = 1;
                    end
                    if high5 > obj.max_slice
                        high5 = obj.max_slice;
                    end
                    hold on;
                    [rowz, ~] = find(filtered_coords(:,3) < low5);
                    if ~isempty(rowz)
                        subset = filtered_coords(rowz, 1:3);
                        plot(subset(:,1), subset(:,2),'marker','o','color',yellow2,'markersize',10,'LineStyle','none'); 
                    end
                    [rowz, ~] = find(filtered_coords(:,3) >= low5 & filtered_coords(:,3) < obj.current_slice);
                    if ~isempty(rowz)
                        subset = filtered_coords(rowz, 1:3);
                        plot(subset(:,1), subset(:,2),'marker','o','color',yellow1,'markersize',10,'LineStyle','none'); 
                    end
                    [rowz, ~] = find(filtered_coords(:,3) > high5);
                    if ~isempty(rowz)
                        subset = filtered_coords(rowz, 1:3);
                        plot(subset(:,1), subset(:,2),'marker','o','color',green2,'markersize',10,'LineStyle','none'); 
                    end
                    [rowz, ~] = find(filtered_coords(:,3) <= high5 & filtered_coords(:,3) > obj.current_slice);
                    if ~isempty(rowz)
                        subset = filtered_coords(rowz, 1:3);
                        plot(subset(:,1), subset(:,2),'markersize',10,'marker','o','color',green1,'LineStyle','none'); 
                    end
                end
            else
                %Okay, let's actually re-render a max proj to eliminate
                %hidden z slices. Since it's been friggin confusing.
                imgstr = obj.img_structs(1);
                if obj.toggle_raw_view
                    imgstr = obj.img_structs(2);
                    max_proj = imgstr.image;
                else
                    visdata = obj.loaded_ch(:,:,obj.z_min:obj.z_max);
                    max_proj = max(visdata,[],3);
                end
                imshow(immultiply(max_proj,immul_mask),[imgstr.Lmin imgstr.Lmax]);
                ref_coords_z = filtered_coords;
            end
            hold on;
            %ref_coords_z = RNA_Threshold_SpotSelector.filterByZ(obj.ref_coords, obj.current_slice);
            if ~isempty(ref_coords_z)
                plot(ref_coords_z(:,1), ref_coords_z(:,2),'or','markersize',10); 
            end
            if obj.mode_3d & obj.toggle_singleSlice
                title(['Spot count: ' num2str(p_count) ' | z = ' num2str(obj.current_slice)]);
            else
                title(['Spot count: ' num2str(p_count)]);
            end
            
        end
        
        %%
        function obj = generateColorTables(obj)
            obj.colortbl_red = RNA_Threshold_SpotSelector.generateColors(obj.max_slice+1, [1.0, 0.0, 0.0]);
            obj.colortbl_magenta = RNA_Threshold_SpotSelector.generateColors(obj.max_slice+1, [1.0, 0.0, 1.0]);
            obj.colortbl_yellow = RNA_Threshold_SpotSelector.generateColors(obj.max_slice+1, [1.0, 1.0, 0.0]);
            obj.colortbl_cyan = RNA_Threshold_SpotSelector.generateColors(obj.max_slice+1, [0.0, 1.0, 1.0]);
            obj.colortbl_white = RNA_Threshold_SpotSelector.generateColors(obj.max_slice+1, [1.0, 1.0, 1.0]);
        end
        
        %%
        function obj = onReadyKey(obj)

            [x,y,btn] = ginput_color(1, obj.crosshair_color);
            %btn
            
            if btn == 1 %Mouse click (set/unset f-pos)
                %Loop until any key pressed that isn't mouse
                obj = obj.whileMouseListening(x,y,btn);
                obj.f_scores_dirty = true;
            elseif  btn == 3 %Right click (set/unset f-neg)
                obj = obj.whileMouseListening(x,y,btn);
                obj.f_scores_dirty = true;
            elseif btn == 'x' %120 'x' - exit
                close(obj.fh_filter);
                close(obj.fh_raw);
                if obj.slice_drawn ~= 0
                    close(obj.fh_origslice);
                end
                if obj.mode_3d
                    close(obj.fh_3d);
                end
                obj.loop_breaker = 1;
            elseif btn == 's' %115 's' - save
                obj = obj.saveMe();
            elseif btn == 'c' %99 'c' - clear (this threshold)
                obj = obj.clearAtThreshold(obj.threshold_idx);
                obj = obj.drawImages();
                obj.threshold_table(obj.threshold_idx, 2) = 1;
                obj.f_scores_dirty = true;
            elseif btn == 'e' %'e' - clear (all but this threshold)
                fprintf("Clearing all but current threshold...\n");
                t_count = size(obj.threshold_table,1);
                for i = 1:t_count
                    if i ~= obj.threshold_idx
                        obj = obj.clearAtThreshold(i);
                    end
                end
                obj.f_scores_dirty = true;
                fprintf("Done!\n");
            elseif btn == 'a' %'a' - auto-reset (this threshold)
                obj = obj.clearAtThreshold(obj.threshold_idx);
                ref_th_idx = RNA_Threshold_SpotSelector.findAReferenceIndex(obj.threshold_table, obj.threshold_idx);
                if ref_th_idx
                    obj = obj.autoClassify(ref_th_idx, obj.threshold_idx);
                else
                    obj = obj.clearAtThreshold(obj.threshold_idx);
                end
                obj = obj.drawImages();
                obj.threshold_table(obj.threshold_idx, 2) = 1;
                obj.f_scores_dirty = true;
            elseif btn == 'r' %'r' - ref
                obj = obj.checkAgainstRef();
                obj = obj.drawImages();
                obj.threshold_table(obj.threshold_idx, 2) = 1;
                obj.f_scores_dirty = true;
            elseif btn == 'p' %'p' - print counts
                obj = obj.printCounts();
            elseif btn == 'P' %'P' - print counts+
                [obj, ctmask] = obj.genCountMask();
                [obj, tpos, fpos, fneg, maskedout] = obj.takeCounts(ctmask);
                tpos
                fpos
                fneg
                maskedout
                if ~isempty(obj.selmcoords)
                    x1 = obj.selmcoords(1,1);
                    x2 = obj.selmcoords(2,1);
                    y1 = obj.selmcoords(3,1);
                    y2 = obj.selmcoords(4,1);
                    fprintf("Mask: (%d,%d)(%d,%d)", x1, y1, x2, y2);
                end
            elseif btn == '^' %'^' - Increase Z trim
                if obj.toggle_change_zminmax
                    if obj.z_max >= obj.max_slice
                        obj.z_max = obj.max_slice;
                    else
                        obj.z_max = obj.z_max + 1;
                    end
                else
                    if obj.z_min >= obj.z_max
                        obj.z_min = obj.z_max;
                    else
                        obj.z_min = obj.z_min + 1;
                    end
                end
                obj.f_scores_dirty = true;
                fprintf("Z Range Updated: %d - %d\n", obj.z_min, obj.z_max);
                obj = obj.drawImages();
            elseif btn == 'v' %'v' - Decrease Z trim
                if obj.toggle_change_zminmax
                    if obj.z_max <= obj.z_min
                        obj.z_max = obj.z_min;
                    else
                        obj.z_max = obj.z_max - 1;
                    end
                else
                    if obj.z_min <= 1
                        obj.z_min = 1;
                    else
                        obj.z_min = obj.z_min - 1;
                    end
                end
                obj.f_scores_dirty = true;
                fprintf("Z Range Updated: %d - %d\n", obj.z_min, obj.z_max);
                obj = obj.drawImages();
            elseif btn == 'Z' %'Z' - Toggle z trim control min/max
                obj.toggle_change_zminmax = ~obj.toggle_change_zminmax;
            elseif btn == '1' %'1' - Remove sample mask
                obj.selmcoords = [];
                obj.f_scores_dirty = true;
                obj = obj.drawImages();
            elseif btn == '4' %'4' - Apply 1/4 sample mask at clicked point
                obj = obj.whileMouseListening_selectMask(2);
                obj.f_scores_dirty = true;
                obj = obj.drawImages();
            elseif btn == '6' %'6' - Apply 1/6 sample mask at clicked point
                obj = obj.whileMouseListening_selectMask(2.45);
                obj.f_scores_dirty = true;
                obj = obj.drawImages();
            elseif btn == '8' %'8' - Apply 1/8 sample mask at clicked point
                obj = obj.whileMouseListening_selectMask(2.83);
                obj.f_scores_dirty = true;
                obj = obj.drawImages();
            elseif btn == 'm' %'m' - toggle single plane/max proj view
                obj.toggle_singleSlice = ~obj.toggle_singleSlice;
                obj.toggle_allz = ~obj.toggle_singleSlice;
                obj = obj.drawImages();
            elseif btn == 'z' %'z' - toggle all-z selection
                obj.toggle_allz = ~obj.toggle_allz;
                obj = obj.drawImages();
            elseif btn == '#' %'#' - toggle spot counts 2D/3D? (I think, I'm not sure what I meant :P)
                obj.toggle_3dcount = ~obj.toggle_3dcount;
                obj = obj.drawImages();
            elseif btn == '%' %'%' - toggle z color scale (nearest 5 planes/all planes)
                obj.toggle_clr_local = ~obj.toggle_clr_local;
                obj = obj.drawImages();
            elseif btn == 'C' %'C' - toggle contrast scale
                obj.toggle_cscale_max = ~obj.toggle_cscale_max;
                obj = obj.drawImages();
            elseif btn == '+' %'+' - Up one slice
                if ~isempty(obj.loaded_ch)
                    Z = size(obj.loaded_ch, 3);
                    if (obj.current_slice < Z)
                        obj.current_slice = obj.current_slice + 1;
                        obj = obj.drawImages();
                    else
                        fprintf("Already at top Z slice!\n");
                    end
                else
                    %fprintf("No source image loaded!\n");
                    %Check to see if it's in 3D mode...
                    if obj.mode_3d
                        if obj.current_slice == obj.max_slice
                            fprintf("Already at top Z slice!\n");
                        else
                            obj.current_slice = obj.current_slice + 1;
                            obj = obj.drawImages();
                        end
                    else
                        fprintf("Source image required for 2D mode!\n");
                    end
                end
            elseif btn == '-' %'-' - Down one slice
                if ~isempty(obj.loaded_ch)
                    if (obj.current_slice > 1)
                        obj.current_slice = obj.current_slice - 1;
                        obj = obj.drawImages();
                        %th_idx = obj.threshold_idx;
                        %[obj, t_pos, f_pos, f_neg, e_pos] = obj.splitCoordTables(th_idx);
                        %obj = obj.drawOnlySlice(f_pos, e_pos, f_neg, t_pos);
                        %obj.fh_filter;
                    else
                        fprintf("Already at bottom Z slice!\n");
                    end
                else
                    %fprintf("No source image loaded!\n");
                    if obj.mode_3d
                        if obj.current_slice <= 1
                            fprintf("Already at bottom Z slice!\n");
                        else
                            obj.current_slice = obj.current_slice - 1;
                            obj = obj.drawImages();
                        end
                    else
                        fprintf("Source image required for 2D mode!\n");
                    end
                end
            elseif btn == '=' %'=' - Up 10 slices
                if ~isempty(obj.loaded_ch)
                    Z = size(obj.loaded_ch, 3);
                    if (obj.current_slice < Z)
                        obj.current_slice = min(obj.current_slice + 10, obj.max_slice);
                        obj = obj.drawImages();
                        %th_idx = obj.threshold_idx;
                        %[obj, t_pos, f_pos, f_neg, e_pos] = obj.splitCoordTables(th_idx);
                        %obj = obj.drawOnlySlice(f_pos, e_pos, f_neg, t_pos);
                        %obj.fh_filter;
                    else
                        fprintf("Already at top Z slice!\n");
                    end
                else
                    %fprintf("No source image loaded!\n");
                    %Check to see if it's in 3D mode...
                    if obj.mode_3d
                        if obj.current_slice == obj.max_slice
                            fprintf("Already at top Z slice!\n");
                        else
                            obj.current_slice = min(obj.current_slice + 10, obj.max_slice);
                            obj = obj.drawImages();
                        end
                    else
                        fprintf("Source image required for 2D mode!\n");
                    end
                end
            elseif btn == '_' %'_' - Down ten slices
                if ~isempty(obj.loaded_ch)
                    if (obj.current_slice > 1)
                        obj.current_slice = max(obj.current_slice - 10, 1);
                        obj = obj.drawImages();
                        %th_idx = obj.threshold_idx;
                        %[obj, t_pos, f_pos, f_neg, e_pos] = obj.splitCoordTables(th_idx);
                        %obj = obj.drawOnlySlice(f_pos, e_pos, f_neg, t_pos);
                        %obj.fh_filter;
                    else
                        fprintf("Already at bottom Z slice!\n");
                    end
                else
                    %fprintf("No source image loaded!\n");
                    if obj.mode_3d
                        if obj.current_slice <= 1
                            fprintf("Already at bottom Z slice!\n");
                        else
                            obj.current_slice = max(obj.current_slice - 10, 1);
                            obj = obj.drawImages();
                        end
                    else
                        fprintf("Source image required for 2D mode!\n");
                    end
                end
            elseif btn == '>' %'>' - th_idx + 1
                th_max = size(obj.threshold_table, 1);
                new_idx = obj.threshold_idx + 1;
                if new_idx > th_max
                    new_idx = th_max;
                end
                obj = obj.changeThreshold(new_idx);
            elseif btn == '<' %'<' - th_idx - 1
                new_idx = obj.threshold_idx - 1;
                if new_idx < 0
                    new_idx = 0;
                end
                obj = obj.changeThreshold(new_idx);
            elseif btn == ']' %']' - th_idx + 10
                th_max = size(obj.threshold_table, 1);
                new_idx = obj.threshold_idx + 10;
                if new_idx > th_max
                    new_idx = th_max;
                end
                obj = obj.changeThreshold(new_idx);
            elseif btn == '[' %'[' - th_idx - 10
                new_idx = obj.threshold_idx - 10;
                if new_idx < 0
                    new_idx = 0;
                end
                obj = obj.changeThreshold(new_idx);
            end

        end
        
        %%
        function obj = onReadyKey_RefMode(obj)
            [x,y,btn] = ginput_color(1, obj.crosshair_color);
            
            if btn == 1 %Mouse click (set/unset f-pos)
                %Loop until any key pressed that isn't mouse
                obj = obj.whileMouseListening_RefMode(x,y,btn);
                obj.ref_last_modified = datetime;
            elseif btn == 'm' %106 'm' - toggle single plane/max proj view
                obj.toggle_singleSlice = ~obj.toggle_singleSlice;
                obj.toggle_allz = ~obj.toggle_singleSlice;
                fprintf("Max proj mode: Single Slice Toggle: %d, All Z Toggle: %d\n", obj.toggle_singleSlice,obj.toggle_allz);
                obj = obj.drawRefImage();
            elseif btn == 'f' %'f' - toggle raw/filtered view
                obj.toggle_raw_view = ~obj.toggle_raw_view;
                if obj.toggle_raw_view
                    fprintf("Switching to raw image view...\n");
                else
                    fprintf("Switching to filtered image view...\n");
                end
                obj = obj.drawRefImage();
            elseif btn == 'S' %83 'S' - snap
                %Can only do if not refonly
                if isempty(obj.positives)
                    fprintf("Can't snap if no auto set to snap to!\n");
                else
                    obj = obj.refSnapToAutoSpots();
                    obj.ref_last_modified = datetime;
                    obj = obj.drawRefImage();
                end
            elseif btn == 'U' %'U' - toggle remove unsnapped
                obj.toggle_del_unsnapped = ~obj.toggle_del_unsnapped;
                obj = obj.drawRefImage();
            elseif btn == 'x' %'x' - exit
                close(obj.fh_filter);
                obj.loop_breaker = 1;
            elseif btn == 's' %'s' - save
                obj = obj.saveMe();
            elseif btn == 'c' %'c' - clear
                obj.ref_coords = [];
                obj.ref_last_modified = datetime;
                obj = obj.drawRefImage();
            elseif btn == '^' %'^' - Increase Z trim
                if obj.toggle_change_zminmax
                    if obj.z_max >= obj.max_slice
                        obj.z_max = obj.max_slice;
                    else
                        obj.z_max = obj.z_max + 1;
                    end
                else
                    if obj.z_min >= obj.z_max
                        obj.z_min = obj.z_max;
                    else
                        obj.z_min = obj.z_min + 1;
                    end
                end
                obj.f_scores_dirty = true;
                fprintf("Z Range Updated: %d - %d\n", obj.z_min, obj.z_max);
                obj = obj.drawRefImage();
            elseif btn == 'v' %'v' - Decrease Z trim
                if obj.toggle_change_zminmax
                    if obj.z_max <= obj.z_min
                        obj.z_max = obj.z_min;
                    else
                        obj.z_max = obj.z_max - 1;
                    end
                else
                    if obj.z_min <= 1
                        obj.z_min = 1;
                    else
                        obj.z_min = obj.z_min - 1;
                    end
                end
                obj.f_scores_dirty = true;
                fprintf("Z Range Updated: %d - %d\n", obj.z_min, obj.z_max);
                obj = obj.drawRefImage();
            elseif btn == 'Z' %'Z' - Toggle z trim control min/max
                obj.toggle_change_zminmax = ~obj.toggle_change_zminmax;
            elseif btn == '1' %'1' - Remove sample mask
                obj.selmcoords = [];
                obj.f_scores_dirty = true;
                obj = obj.drawRefImage();
            elseif btn == '4' %'4' - Apply 1/4 sample mask at clicked point
                obj = obj.whileMouseListening_selectMask(2);
                obj.f_scores_dirty = true;
                obj = obj.drawRefImage();
            elseif btn == '6' %'6' - Apply 1/6 sample mask at clicked point
                obj = obj.whileMouseListening_selectMask(2.45);
                obj.f_scores_dirty = true;
                obj = obj.drawRefImage();
            elseif btn == '8' %'8' - Apply 1/8 sample mask at clicked point
                obj = obj.whileMouseListening_selectMask(2.83);
                obj.f_scores_dirty = true;
                obj = obj.drawRefImage();
            elseif btn == 'C' %'C' - toggle contrast scale
                obj.toggle_cscale_max = ~obj.toggle_cscale_max;
                obj = obj.drawRefImage();
            elseif btn == '+' %'+' - Up one slice
                if obj.toggle_singleSlice
                    Z = size(obj.loaded_ch, 3);
                    if (obj.current_slice < Z)
                        obj.current_slice = obj.current_slice + 1;
                        obj = obj.drawRefImage();
                    else
                        fprintf("Already at top Z slice!\n");
                    end
                end
            elseif btn == '-' %'-' - Down one slice
                if obj.toggle_singleSlice
                    if (obj.current_slice > 1)
                        obj.current_slice = obj.current_slice - 1;
                        obj = obj.drawRefImage();
                    else
                        fprintf("Already at bottom Z slice!\n");
                    end
                end
            elseif btn == '=' %'=' - Up 10 slices
                if obj.toggle_singleSlice
                    Z = size(obj.loaded_ch, 3);
                    if (obj.current_slice < Z)
                        if obj.current_slice + 10 <= Z
                            obj.current_slice = obj.current_slice + 10;
                        else
                            obj.current_slice = Z;
                        end
                        obj = obj.drawRefImage();
                    else
                        fprintf("Already at top Z slice!\n");
                    end
                end
            elseif btn == '_' %'_' - Down ten slices
                if obj.toggle_singleSlice
                    if (obj.current_slice > 1)
                        if obj.current_slice - 10 > 1
                            obj.current_slice = obj.current_slice - 10;
                        else
                            obj.current_slice = 1;
                        end
                        obj = obj.drawRefImage();
                    else
                        fprintf("Already at bottom Z slice!\n");
                    end
                end
            end
        end
        
        %%
        function obj = changeThreshold(obj, new_th_idx)
            %See if new threshold has been manually modified
            if obj.threshold_table(new_th_idx, 2)
                %User modified. Load as is.
                obj.threshold_idx = new_th_idx;
                obj = obj.drawImages();
            else
                %Do auto classifying before loading
                %Find ref...
                ref_th_idx = RNA_Threshold_SpotSelector.findAReferenceIndex(obj.threshold_table, new_th_idx);
                if ref_th_idx
                    obj = obj.autoClassify(ref_th_idx, new_th_idx);
                end
                obj.threshold_idx = new_th_idx;
                obj = obj.drawImages();
                %Mark as user modified
                obj.threshold_table(new_th_idx, 2) = 1;
            end
            
        end
        
        %%
        function obj = whileMouseListening(obj, x1, y1, b1)
            
            fullX = size(obj.img_structs(1).image, 2);
            fullY = size(obj.img_structs(1).image, 1);
            rad = RNA_Threshold_SpotSelector.calculateClickRadius(gcf, fullX, fullY);
            
            if b1 == 1
                obj = obj.onLeftClick(x1, y1, rad);
            elseif b1 == 3
                obj = obj.onRightClick(x1, y1, rad);
            end
            
            %Loop
            loopy = 1;
            while loopy == 1
                [x,y,btn] = ginput_color(1, obj.crosshair_color);
                %fprintf("Click detected at %f,%f\n", x, y);
                if btn == 1
                    obj = obj.onLeftClick(x, y, rad);
                elseif btn == 3
                    obj = obj.onRightClick(x, y, rad);
                else
                    loopy = 0;
                end
            end
            
            %Save to obj
            obj = obj.flushNegTable();
            obj = obj.drawImages();
            
        end
        
        %%
        function obj = whileMouseListening_RefMode(obj, x1, y1, b1)
            
            fullX = size(obj.img_structs(1).image, 2);
            fullY = size(obj.img_structs(1).image, 1);
            rad = RNA_Threshold_SpotSelector.calculateClickRadius(gcf, fullX, fullY);
            
            if b1 == 1
                obj = obj.onLeftClick_RefMode(x1, y1, rad);
            end
            
            %Loop
            loopy = 1;
            while loopy == 1
                [x,y,btn] = ginput_color(1, obj.crosshair_color);
                %fprintf("Click detected at %f,%f\n", x, y);
                if btn == 1
                    obj = obj.onLeftClick_RefMode(x, y, rad);
                else
                    loopy = 0;
                end
            end
            
            %Save to obj
            obj = obj.flushRefTable();
            obj = obj.drawRefImage();
            
        end
        
        %%
        function obj = whileMouseListening_selectMask(obj, denom)

            [x1,y1,~] = ginput_color(1, obj.crosshair_color);
            fullX = size(obj.img_structs(1).image, 2);
            fullY = size(obj.img_structs(1).image, 1);
            
            %dimdenom = denom./2;
            dimdenom = denom;
            
            xsz = fullX./dimdenom;
            ysz = fullY./dimdenom;
            
            x_left = x1 - (xsz./2);
            y_up = y1 - (ysz./2);
            x_left = max(1, x_left);
            y_up = max(1, y_up);
            
            x_right = x_left + xsz;
            y_down = y_up + ysz;
            
            if(x_right > fullX)
                x_right = fullX;
                x_left = x_right - xsz;
            end
            
            if(y_down > fullY)
                y_down = fullY;
                y_up = y_down - ysz;
            end
            
            obj.selmcoords = zeros(4,1);
            obj.selmcoords(1,1) = uint16(x_left);
            obj.selmcoords(2,1) = uint16(x_right);
            obj.selmcoords(3,1) = uint16(y_up);
            obj.selmcoords(4,1) = uint16(y_down);
            
        end
        
        %%
        %TODO
        function obj = onLeftClick(obj, x, y, r)
            
            %Find matching points
            c_table = obj.positives{obj.threshold_idx};
            [~,idxs] = obj.findSelectedSpotsObj(x, y, c_table, r);
            
            %For all selected points...
            if ~isempty(idxs)
                m_count = size(idxs,1);
                for mi = 1:m_count
                    idx = idxs(mi);
                    %Reverse value in column 4 of table
                    old_val = c_table(idx,4);
                    if old_val == 0
                        %fpos -> epos
                        c_table(idx,4) = 2;
                        new_marker = 'om';
                    elseif old_val == 1
                        %tpos -> fpos
                        c_table(idx,4) = 0;
                        new_marker = 'oy';
                    elseif old_val == 2
                        %epos -> tpos
                        c_table(idx,4) = 1;
                        new_marker = 'or';
                    end
                    %Draw a new circle in the appropriate color
                    obj.drawSingleCircle(c_table(idx,1), c_table(idx,2), new_marker);
                end    
                obj.positives{obj.threshold_idx} = c_table;
            end
            
        end
        
        %%
        %TODO
        function obj = onRightClick(obj, x, y, r)
            
            %Check to see if there is a temp table in use. If not, fetch.
            if obj.n_alloc == 0
                if isempty(obj.false_negs{obj.threshold_idx})
                    obj.n_table = zeros([256 4]);
                    obj.n_table(:,4) = 1;
                    obj.n_alloc = 256;
                    obj.n_used = 0;
                else
                    fneg_tbl = obj.false_negs{obj.threshold_idx};
                    fn_tbl_sz = size(fneg_tbl,1);
                    %fprintf("fn_tbl_sz = %d\n", fn_tbl_sz);
                    tbl_sz = fn_tbl_sz + 256;
                    obj.n_table = zeros([tbl_sz 4]);
                    obj.n_table(1:fn_tbl_sz,1:3) = fneg_tbl(:,1:3);
                    obj.n_table(fn_tbl_sz+1:tbl_sz,4) = 1;
                    obj.n_alloc = tbl_sz;
                    obj.n_used = fn_tbl_sz;
                end
            end
            
            %Test for matching points
            match_found = 0;
            
            if obj.n_used > 0
                %idxs = RNA_Threshold_SpotSelector.findSelectedSpots(x, y, 0, obj.n_table, r);
                [~,idxs] = obj.findSelectedSpotsObj(x, y, c_table, r);
                if ~isempty(idxs)
                    match_found = 1;
                    %Reverse flag and re-render circle
                    m_count = size(idxs,1);
                    for mi = 1:m_count
                        idx = idxs(mi);
                        %Reverse value in column 3 of table
                        old_val = obj.n_table(idx,4);
                        if old_val == 0
                            %fneg -> tneg (Flag for del)
                            obj.n_table(idx,4) = 1;
                            new_marker = 'ow';
                        else
                            %tneg -> fneg (Unflag)
                            obj.n_table(idx,4) = 0;
                            new_marker = 'oc';
                        end
                        %Draw a new circle in the appropriate color
                        obj.drawSingleCircle(obj.n_table(idx,1), obj.n_table(idx,2), new_marker);
                    end 
                end
            end
            
            %If none, create a new false neg
            c_table = obj.positives{obj.threshold_idx};
            %idxs = RNA_Threshold_SpotSelector.findSelectedSpots(x, y, c_table, r);
            [~,idxs] = obj.findSelectedSpotsObj(x, y, c_table, r);
            if isempty(idxs) && match_found == 0
                obj = obj.addNeg(x,y);
                obj.drawSingleCircle(x, y, 'oc');
            end
            
        end
        
        %%
        function [obj, idxs] = findSelectedSpotsObj(obj, x, y, c_table, r)
            if obj.mode_3d
                if obj.toggle_allz
                    idxs = RNA_Threshold_SpotSelector.findSelectedSpots(x, y, 0, c_table, r);
                else
                    idxs = RNA_Threshold_SpotSelector.findSelectedSpots(x, y, obj.current_slice, c_table, r);
                end
            else
                idxs = RNA_Threshold_SpotSelector.findSelectedSpots(x, y, 0, c_table, r);
            end
        end
        
        %%
        %TODO TWO ISSUES: r is way too small?, toggle_allz is not set in
        %max proj mode, but it should be.
        function obj = onLeftClick_RefMode(obj, x, y, r)
            
            %Check to see if there is a temp table in use. If not, fetch.
            if obj.n_alloc == 0
                if isempty(obj.ref_coords)
                    obj.n_table = zeros([256 4]);
                    obj.n_table(:,4) = 1;
                    obj.n_alloc = 256;
                    obj.n_used = 0;
                else
                    temp_tbl = obj.ref_coords;
                    temp_tbl_sz = size(temp_tbl,1);
                    %fprintf("fn_tbl_sz = %d\n", fn_tbl_sz);
                    tbl_sz = temp_tbl_sz + 256;
                    obj.n_table = zeros([tbl_sz 4]);
                    obj.n_table(1:temp_tbl_sz,1:3) = temp_tbl(:,1:3);
                    obj.n_table(temp_tbl_sz+1:tbl_sz,4) = 1;
                    obj.n_alloc = tbl_sz;
                    obj.n_used = temp_tbl_sz;
                end
            end
            
            %Test for matching points
            %match_found = 0;
            
            if obj.n_used > 0
                if obj.toggle_allz
                    idxs = RNA_Threshold_SpotSelector.findSelectedSpots(x, y, 0, obj.n_table, r);
                else
                    idxs = RNA_Threshold_SpotSelector.findSelectedSpots(x, y, obj.current_slice, obj.n_table, r);
                end
                if ~isempty(idxs)
                    %match_found = 1;
                    %Reverse flag and re-render circle
                    m_count = size(idxs,1);
                    for mi = 1:m_count
                        idx = idxs(mi);
                        %Reverse value in column 3 of table
                        old_val = obj.n_table(idx,4);
                        if old_val == 0
                            %fneg -> tneg (Flag for del)
                            obj.n_table(idx,4) = 1;
                            new_marker = 'ow';
                        else
                            %tneg -> fneg (Unflag)
                            obj.n_table(idx,4) = 0;
                            new_marker = 'oc';
                        end
                        %Draw a new circle in the appropriate color
                        obj.drawSingleCircle(obj.n_table(idx,1), obj.n_table(idx,2), new_marker);
                    end 
                    return;
                end
            end
            
            if (x < 1) | (y < 1)
                return;
            end
            
            %If none, create a new false neg
            obj = obj.addNeg(x,y);
            obj.drawSingleCircle(x, y, 'oc');
            
        end
        
        %%
        function obj = autoClassify(obj, ref_th_idx, target_th_idx)
 
            %Clear
            obj = obj.clearAtThreshold(target_th_idx);
            
            %Common
            [obj, rt_pos, rf_pos, rf_neg, re_pos] = obj.splitCoordTables(ref_th_idx);
            full_pos = obj.positives{target_th_idx};
            tt_pos = full_pos(:,1:3);
            
            %Any false neg in ref that are not pos in targ copy to false
            %neg in targ
            if ~isempty(rf_neg)
                [~, not_idxs] = RNA_Threshold_SpotSelector.findRowMatches(rf_neg, tt_pos);
                tf_neg = zeros([size(not_idxs,1) 3]);
                tf_neg(:,1:3) = rf_neg(not_idxs,1:3);
                obj.false_negs{target_th_idx} = tf_neg;
            end
            
            if ref_th_idx < target_th_idx
               %If true positive at lower thresh ref...
               %And negative at higher thresh target...
               %Make a false negative at target threshold
               
               %Find true positives in ref that are not positive in targ
               [~, not_idxs] = RNA_Threshold_SpotSelector.findRowMatches(rt_pos, tt_pos);
               if ~isempty(not_idxs)
                    tf_neg_add = rt_pos(not_idxs,:);
                    tf_neg = obj.false_negs{target_th_idx};
                    if isempty(tf_neg)
                        %Set new additions
                        obj.false_negs{target_th_idx} = tf_neg_add;
                    else
                        %Add to table
                        sz_1 = size(tf_neg,1);
                        sz_2 = size(tf_neg_add,1);
                        
                        new_tbl = zeros(sz_1 + sz_2,3);
                        new_tbl(1:sz_1,:) = tf_neg(:,:);
                        new_tbl(sz_1+1:sz_1+sz_2,:) = tf_neg_add(:,:);
                        obj.false_negs{target_th_idx} = new_tbl;
                    end
               end
            elseif ref_th_idx > target_th_idx
                %If tpos in targ, but neither fneg or tpos in ref,
                % set to fpos
                %Find positives that aren't tpos in reference
                [~, not_idxs] = RNA_Threshold_SpotSelector.findRowMatches(tt_pos, rt_pos);
                if ~isempty(not_idxs)
                    %Now take out the ones in fneg...
                    p_list = tt_pos(not_idxs,:);
                    [~, not_idxs] = RNA_Threshold_SpotSelector.findRowMatches(p_list, rf_neg);
                    if ~isempty(not_idxs)
                        p_list = p_list(not_idxs,:);
                    end
                    [idx_list, ~] = RNA_Threshold_SpotSelector.findRowMatches(tt_pos, p_list);
                    obj.positives{target_th_idx}(idx_list,4) = 0;
                end
            end
            
            %Any false pos in ref that are true pos in targ switch to false
            %pos
            %Same for excluded
            %row_matches = ismember(tt_pos, rf_pos, 'rows');
            [match_idxs, ~] = RNA_Threshold_SpotSelector.findRowMatches(tt_pos, rf_pos);
            obj.positives{target_th_idx}(match_idxs,4) = 0;
            [match_idxs, ~] = RNA_Threshold_SpotSelector.findRowMatches(tt_pos, re_pos);
            obj.positives{target_th_idx}(match_idxs,4) = 2;

        end
        
        %%
        function obj = checkAgainstRef(obj)
            %Make this (and findRowMatches) more efficient for HUGE
            %spotcall sets
 
            %Clear
            target_th_idx = obj.threshold_idx;
            obj = obj.clearAtThreshold(target_th_idx);
            if isempty(obj.ref_coords)
                return;
            end
            fprintf('DEBUG -- Checking spots against ref for threshold index %d...\n', target_th_idx);
            
            r_tbl = obj.ref_coords;
            th_pos = obj.positives{target_th_idx};
            if ~isempty(th_pos)
                %Presort th_pos(?)
                [~, si] = sort(th_pos(:,3));
                th_pos = th_pos(si, :);
                [~, si] = sort(th_pos(:,2));
                th_pos = th_pos(si, :);
                [~, si] = sort(th_pos(:,1));
                th_pos = th_pos(si, :);
                
                %Anything in pos, but not in ref, set to false pos
                th_pos(:,4) = 1;
                mres = ~ismember(th_pos(:,1:3), r_tbl, 'rows');
                not_idxs = find(mres);
                if ~isempty(not_idxs)
                    th_pos(not_idxs,4) = 0;
                end
                obj.positives{target_th_idx} = th_pos;

                %Anything in ref, but not in pos, set to false neg
                mres = ~ismember(r_tbl, th_pos(:,1:3), 'rows');
                not_idxs = find(mres);
                if ~isempty(not_idxs)
                    neg_tbl = zeros(size(not_idxs,1), 4);
                    neg_tbl(:,1:3) = r_tbl(not_idxs,:);
                    obj.false_negs{target_th_idx} = neg_tbl;
                end
            else
                rcount = size(r_tbl,1);
                neg_tbl = zeros(size(rcount,1), 4);
                neg_tbl(:,1:3) = r_tbl(rcount,:);
                obj.false_negs{target_th_idx} = neg_tbl;
            end
        end
        
        %%
        function obj = clearAtThreshold(obj, target_th_idx)
 
            obj.false_negs{target_th_idx} = [];
            if ~isempty(obj.positives{target_th_idx})
                obj.positives{target_th_idx}(:,4) = 1;
            end

            obj.threshold_table(target_th_idx,2) = 0;
            
        end
        
        %%
        function obj = setMaxZ(obj)
            %Scan coord tables to find maximum z coordinate
            maxz = 0;
            tcount = size(obj.threshold_table, 1);
            
            for t = 1:tcount
                %Pos table
                ptbl = obj.positives{t};
                scount = size(ptbl, 1);
                for s = 1:scount
                    if ptbl(s,3) > maxz
                        maxz = ptbl(s,3);
                    end
                end
            
                %Neg table
                if ~isempty(obj.false_negs{t})
                    ntbl = obj.false_negs{t};
                    scount = size(ntbl, 1);
                    for s = 1:scount
                        if ntbl(s,3) > maxz
                            maxz = ntbl(s,3);
                        end
                    end
                end
            end
            
            obj.max_slice = maxz;
        end
        
        %%
        function obj = printCounts(obj)
            [obj, tpos, fpos, fneg, epos] = obj.splitCoordTables();
            
            if isempty(tpos)
                tp_count = 0;
            else
                tp_count = size(tpos,1);
            end
            
            if isempty(fpos)
                fp_count = 0;
            else
                fp_count = size(fpos,1);
            end
            
            if isempty(fneg)
                fn_count = 0;
            else
                fn_count = size(fneg,1);
            end
            
            if isempty(epos)
                ep_count = 0;
            else
                ep_count = size(epos,1);
            end
            
            
            fprintf("Threshold: %d\n", obj.threshold_table(obj.threshold_idx, 1));
            fprintf("True Positives: %d\n", tp_count);
            fprintf("False Positives: %d\n", fp_count);
            fprintf("False Negatives: %d\n", fn_count);
            fprintf("Excluded: %d\n", ep_count);
        end
        
        %%
        function obj = updateFTable(obj)
            [obj, ctmask] = obj.genCountMask();
            [obj, tpos, fpos, fneg, ~] = obj.takeCounts(ctmask);
            obj.f_scores(:,2) = tpos(:,1);
            obj.f_scores(:,3) = fpos(:,1);
            obj.f_scores(:,4) = fneg(:,1);
            totalpos = tpos + fneg;
            totalpos(totalpos == 0) = NaN;
            totalhits = tpos + fpos;
            totalhits(totalhits == 0) = NaN;
            sensitivity = tpos./totalpos;
            precision = tpos./totalhits;
            obj.f_scores(:,1) = (2 .* sensitivity .* precision)./(sensitivity + precision);
            obj.f_scores_dirty = false;
        end
        
        %%
        function obj = loadSourceImageChannel(obj)
            
            if isempty(obj.tif_path)
                fprintf("No TIFF path set!\n");
                return;
            end
            
            [stack, img_read] = tiffread2(obj.tif_path);
            Z = img_read/obj.tif_chcount;
            Y = size(stack(1,1).data,1);
            X = size(stack(1,1).data,2);
            
            obj.loaded_ch = NaN(Y,X,Z);
            idx = obj.tif_ch;
            for z = 1:Z
                 obj.loaded_ch(:,:,z) = stack(1,idx).data;
                idx = idx + obj.tif_chcount;
            end
            
            fprintf("TIFF channel loaded!\n");
            obj.current_slice = round(Z./2);
            
            RNA_Threshold_Common.saveDeadPixels(obj.loaded_ch);
            obj.loaded_ch = RNA_Threshold_Common.cleanDeadPixels(obj.loaded_ch);
            
        end
        
        %%
        function [obj, ctmask] = genCountMask(obj)

            Y = size(obj.img_structs(1).image,1);
            X = size(obj.img_structs(1).image,2);
            Z = obj.max_slice;
                
            if ~isempty(obj.selmcoords)
                %apply XY mask
                x1 = obj.selmcoords(1,1);
                x2 = obj.selmcoords(2,1);
                y1 = obj.selmcoords(3,1);
                y2 = obj.selmcoords(4,1);
                ctmask = false(Y,X,Z);
                ctmask(y1:y2,x1:x2,:) = true;
            else
                ctmask = true(Y,X,Z);
            end
            
            if obj.z_min > 1
                ctmask(:,:,1:obj.z_min) = false;
            end
            if obj.z_max < Z
                ctmask(:,:,obj.z_max:Z) = false;
            end
        end
          
        %%
        % (Internal function)
        % In ref mode, snap all spots specified as a reference spot to the
        % nearest auto-detected spot in 3D, starting at the highest
        % threshold a nearby spot can be found.
        %
        function obj = refSnapToAutoSpots(obj, stopAt, maxrad_3d, maxrad_z, verbose)
            %If not in single slice mode, ignore z...
            %Also make sure to ignore masks.
            if nargin < 2
                stopAt = 20;
            end
            if nargin < 3
                maxrad_3d = 4;
            end
            if nargin < 4
                maxrad_z = 2;
            end
            if nargin < 5
                verbose = true;
            end
            
            if isempty(obj.ref_coords)
                return;
            end
            
            T = size(obj.threshold_table, 1);
            for t = 1:T
                obj = obj.clearAtThreshold(t);
            end
            if verbose; fprintf("Max Threshold: %d\n",obj.threshold_table(T,1)); end
            
            ref_spot_count = size(obj.ref_coords,1);
            snapped_tbl = zeros(ref_spot_count,4);
            snapped_tbl(:,1:3) = obj.ref_coords(:,1:3);
            
            for t = T:-1:stopAt
                if verbose; fprintf("Trying threshold %d...\n",obj.threshold_table(t,1)); end
                if isempty(obj.positives{t,1}); continue; end
                
                %Find yet unsnapped spots
                unsnapped_count = nnz(~snapped_tbl(:,4));
                if unsnapped_count < 1
                    %All spots have been snapped. We're done here.
                    break;
                end
                
                unsnapped_rows = find(~snapped_tbl(:,4));
                already_snapped = snapped_tbl(find(snapped_tbl(:,4)), :);
                
                thpos = obj.positives{t,1};
                spotcount = size(thpos,1);
                pos_tbl = NaN(spotcount,6);
                pos_tbl(:,1:3) = thpos(:,1:3);
                
                for ri = 1:unsnapped_count
                    r = unsnapped_rows(ri);
                    
                    x_dist = pos_tbl(:,1) - snapped_tbl(r,1);
                    y_dist = pos_tbl(:,2) - snapped_tbl(r,2);
                    z_dist = pos_tbl(:,3) - snapped_tbl(r,3);
                    
                    %This behavior changes based on whether we are looking
                    %at all z...
                    if obj.toggle_allz
                        pos_tbl(:,4) = 0;
                        pos_tbl(:,5) = sqrt(x_dist.^2 + y_dist.^2);
                        pos_tbl(:,6) = pos_tbl(:,5);
                    else
                        pos_tbl(:,4) = abs(z_dist);
                        pos_tbl(:,5) = sqrt(x_dist.^2 + y_dist.^2);
                        pos_tbl(:,6) = sqrt(x_dist.^2 + y_dist.^2 + z_dist.^2);
                    end
                    
                    match_bool = (pos_tbl(:,6) <= maxrad_3d);
                    match_bool = and(match_bool, (pos_tbl(:,4) <= maxrad_z));
                    
                    if nnz(match_bool) > 0
                        %Isolate matches
                        [match_rows, ~] = find(match_bool);
                        rmatches = pos_tbl(match_rows,:);
                        
                        %Remove any of those already in snapped set
                        if ~isempty(already_snapped)
                            match_bool = ~ismember(rmatches(:,1:3), already_snapped(:,1:3),'rows');
                            if nnz(match_bool) < 1; continue; end
                            [match_rows, ~] = find(match_bool);
                            rmatches = rmatches(match_rows,:);
                        end
                        
                        %Now, snap to nearest match.
                        [~,I] = min(rmatches(:,6));
                        old_pt = snapped_tbl(r,1:3);
                        snapped_tbl(r,1:3) = rmatches(I,1:3);
                        if verbose
                            fprintf('DEBUG -- (%d,%d,%d) snapped to (%d,%d,%d)\n', ...
                            old_pt(1,1), old_pt(1,2), old_pt(1,3),...
                            snapped_tbl(r,1), snapped_tbl(r,2), snapped_tbl(r,3));
                        end
                        snapped_tbl(r,4) = 1;
                    end
                end
            end
            
            %Count still unsnapped spots, print message,
            %   and remove if toggle is on
            unsnapped_count = nnz(~snapped_tbl(:,4));
            if unsnapped_count > 0
                if obj.toggle_del_unsnapped
                    unsnapped_rows = find(~snapped_tbl(:,4));
                    snapped_tbl(unsnapped_rows,:) = [];
                    if verbose
                        fprintf('Removed %d unsnapped spots.\n', unsnapped_count);
                    end
                else
                    if verbose
                        fprintf('%d spots could not be snapped.\n', unsnapped_count);
                    end
                end
            end
            
            %Copy back to refset and return
            obj.ref_coords = snapped_tbl(:,1:3);
            obj.ref_last_modified = datetime;
            obj.f_scores_dirty = true;
        end
        
        %%
        function [obj, tpos, fpos, fneg, maskedout] = takeCounts(obj, mask)
            
            if isempty(obj.positives)
                tpos = 0;
                fpos = 0;
                fneg = 0;
                maskedout = 0;
                return;
            end
            
            T = size(obj.threshold_table, 1);
            tpos = zeros(T,1);
            fpos = zeros(T,1);
            fneg = zeros(T,1);
            maskedout = zeros(T,1);
            
            oldt = obj.threshold_idx;
            %size(mask)
            
            for t = 1:T
                obj.threshold_idx = t;
                obj = obj.checkAgainstRef();
                obj.threshold_table(obj.threshold_idx, 2) = 1;
                
%                 if t >= 533
%                     fprintf("Debug\n");
%                 end
                
                [obj, tp, fp, fn, ~] = obj.splitCoordTables(t, false);
                if isempty(mask)
                    %Can just use find and take sizes
                    if ~isempty(tp)
                        tpos(t) = size(tp, 1);
                    end
                    if ~isempty(fp)
                        fpos(t) = size(fp, 1);
                    end
                    if ~isempty(fn)
                        fneg(t) = size(fn, 1);
                    end
                else
                    %Have to go through each point to see if masked
                    if ~isempty(tp)
                        sz = size(tp, 1);
                        for i = 1:sz
                            if mask(tp(i,2), tp(i,1), tp(i,3))
                                tpos(t) = tpos(t)+1;
                            else
                                maskedout(t) = maskedout(t)+1; 
                            end
                        end
                    end
                    if ~isempty(fp)
                        sz = size(fp, 1);
                        for i = 1:sz
                            if mask(fp(i,2), fp(i,1), fp(i,3))
                                fpos(t) = fpos(t)+1;
                            else
                                maskedout(t) = maskedout(t)+1; 
                            end
                        end
                    end
                    if ~isempty(fn)
                        %Also use this as a cleaning opportunity.
                        fn = obj.false_negs{t}(:,:);
                        tmp = zeros(size(fn,1), size(fn,2));
                        goodcount = 0;
                        sz = size(fn, 1);
                        for i = 1:sz
                            if(fn(i,1) < 1) | (fn(i,2) < 1) | (fn(i,3) < 1)
                                fprintf("WARNING: Bad fneg found (%d,%d,%d) @ th = %d. Will be removed.\n",fn(i,1),fn(i,2),fn(i,3),t);
                            else
                               if mask(fn(i,2), fn(i,1), fn(i,3))
                                    fneg(t) = fneg(t)+1;
                               else
                                    maskedout(t) = maskedout(t)+1; 
                               end 
                               goodcount = goodcount+1;
                               tmp(goodcount,:) = fn(i,:);
                            end
                        end
                        obj.false_negs{t} = tmp(1:goodcount,:);
                    end
                end
            end
            
            obj.threshold_idx = oldt;
            
        end
        
        %%
        function [obj, copy] = makeCopy(obj)
            copy = RNA_Threshold_SpotSelector;
            copy.save_stem = obj.save_stem;
            
            copy.img_structs = obj.img_structs;
            copy.threshold_idx = obj.threshold_idx;
            copy.threshold_table = obj.threshold_table;
            copy.positives = obj.positives;
            copy.false_negs = obj.false_negs;
            
            copy.ref_coords = obj.ref_coords;
            copy.mode_3d = obj.mode_3d;
            copy.current_slice = obj.current_slice;
            copy.max_slice = obj.max_slice;
            copy.imgdat_path = obj.imgdat_path;
            copy.ztrim = obj.ztrim;
            copy.z_min = obj.z_min;
            copy.z_max = obj.z_max;
            copy.selmcoords = obj.selmcoords;
            
            copy.toggle_singleSlice = obj.toggle_singleSlice;
            copy.toggle_allz = obj.toggle_allz;
            copy.toggle_3dcount = obj.toggle_3dcount;
            copy.toggle_clr_local = obj.toggle_clr_local;
            copy.toggle_del_unsnapped = obj.toggle_del_unsnapped;
            copy.toggle_cscale_max = obj.toggle_cscale_max;
            copy.f_scores = obj.f_scores;
            
            copy.loaded_ch = [];
            s_img = copy.img_structs(1).image;
            copy.alloc_size = size(s_img, 1) * size(s_img,2);
            copy.n_alloc = 0;
            copy.n_used = 0;
            copy.n_table = [];
            copy.slice_drawn = 0;
            copy.crosshair_color = obj.crosshair_color;
        end
        
        %%
        function obj = loadNewSpotset(obj, spot_counts, coord_table)
            
            t_count = size(spot_counts,1);
            dimcount = size(coord_table{1});
            
            obj.positives = cell(t_count, 1);
            for t = 1:t_count
                src_tbl = coord_table{t};
                if isempty(src_tbl)
                    obj.positives{t} = [];
                    continue;
                end
                
                tbl = ones([size(src_tbl,1) 4]);
                if dimcount == 2
                    tbl(:,1:2) = src_tbl(:,1:2);
                    obj.mode_3d = false;
                else
                    tbl(:,1:3) = src_tbl(:,1:3);
                    obj.mode_3d = true;
                end
                tbl = uint16(tbl);
                obj.positives{t} = tbl;
            end
            obj.false_negs = cell(t_count, 1);
            
            obj.threshold_table = zeros(t_count,1);
            obj.threshold_table(:,1) = spot_counts(:,1);
            obj.threshold_idx = uint16(t_count/2);
            
            obj.f_scores = NaN(t_count,4);
            obj = obj.updateFTable();
        end
        
        %%
        function idims = getImageDimensions(obj)
            idims = struct('x', 0, 'y', 0, 'z', 0);
            idims.z = obj.max_slice;
            idims.y = size(obj.img_structs(1).image,1);
            idims.x = size(obj.img_structs(1).image,2);
        end
        
        %%
        function scoreTable = genScoreResponseTable(obj)
            %Return table columns:
            %   x,y,z,intensity,call(1 or 0)

            %This does not resnap or reclassify spots
            %For any calls that are false negatives at lowest threshold,
            %just set the intensity value to something below the lowest
            %checked.
            if isempty(obj.ref_coords) 
                scoreTable = [];
                return; 
            end

            %Grab calls at lowest threshold.
            ptbl = obj.positives{1,1};
            call_count = size(ptbl,1);
            ref_count = size(obj.ref_coords,1);

            temp_table = zeros(call_count + ref_count, 5);
            temp_table(1:call_count,1:3) = ptbl(:,1:3);
            temp_table((call_count + 1):(call_count + ref_count), 1:3) = obj.ref_coords(:,1:3);
            
            all_spots = unique(temp_table, 'rows');

            %Mark true/false
            inref = ismember(all_spots(:,1:3), obj.ref_coords(:,1:3), 'rows');
            all_spots(:,5) = inref(:,1);

            %Now go through all thresholds...
            T = size(obj.threshold_table, 1);
            for t = 1:T
                thval = obj.threshold_table(t,1);
                ptbl = obj.positives{t,1};
                if isempty(ptbl); break; end
                inset = ismember(all_spots(:,1:3), ptbl(:,1:3), 'rows');
                rowidxs = find(inset(:,1));
                all_spots(rowidxs, 4) = thval;
            end

            scoreTable = all_spots;
        end
    end
    
    %%
    methods(Static)
        
        %%
        function selector = createEmptyRefSelector(tif_path, total_ch, sample_ch, save_stem)
            selector = RNA_Threshold_SpotSelector;
            
            %Initialize paths
            selector.save_stem = save_stem;
            selector.imgdat_path = [save_stem '_prefilteredIMG'];
            
            %Create image projections
            [channels, idims] = LoadTif(tif_path, total_ch, sample_ch, 1);
            sample_channel = channels{sample_ch,1};
            [img_filter] = RNA_Threshold_SpotDetector.run_spot_detection_pre(sample_channel, save_stem, true);
            save(selector.imgdat_path, 'img_filter');
            load([save_stem '_imgviewstructs'], 'my_images');
            selector.img_structs = my_images;
            
            %Set remaining values
            selector.current_slice = 1;
            selector.max_slice = idims.z;
            selector.z_min = 1;
            selector.z_max = idims.z;
            selector.mode_3d = (idims.z > 1);
            selector.threshold_idx = 0;
            selector.threshold_table = [];
            
            selector.toggle_singleSlice = true;
            selector.toggle_allz = false;
            selector.toggle_3dcount = true;
            selector.toggle_clr_local = false;
            selector.toggle_del_unsnapped = true;
            selector.toggle_cscale_max = true;
            selector.toggle_change_zminmax = false;
            selector.f_scores_dirty = false;
            
            selector.crosshair_color = [0.000, 0.000, 0.000];
        end
        
        %%
        function count = countFromTable(obj, table)
            if isempty(table)
                count = 0;
            else
                if(obj.mode_3d)
                    if(obj.toggle_3dcount)
                        count = size(table,1);
                    else
                        %Filter to unique xy
                        t2 = RNA_Threshold_Common.collapse3DCoordTable(table);
                        if isempty(t2)
                            count = 0;
                        else
                            count = size(t2, 1);
                        end
                    end
                else
                    count = size(table,1);
                end
            end
        end
        
        %%
        function clrmtx = gen_z_color_matrix(color_tbl, coord_tbl, rel_z, range)
            maxidx = size(color_tbl,1);
            incr = 1;
            
            if range > 0
                incr = maxidx./range;
            end
            
            color_idxs = min((abs(rel_z - coord_tbl(:,3)) + 1) * incr, maxidx);
            clrmtx = color_tbl(color_idxs,:);
        end
        
        %%
        function color_tbl = generateColors(levels, base_color)
            
            levels = double(levels+1);
            color_tbl = NaN(levels,3);
            konst = 1/log10(levels);
            
            r = base_color(1,1);
            g = base_color(1,2);
            b = base_color(1,3);
            
            color_tbl(1,1) = r;
            color_tbl(1,2) = g;
            color_tbl(1,3) = b;
            
            %r_div = r/levels;
            %g_div = g/levels;
            %b_div = b/levels;
            
            for l = 2:levels
                %r = r-r_div;
                r = (1 - (log10(l) * konst)) * base_color(1,1);
                color_tbl(l,1) = r;
                
                %g = g-g_div;
                g = (1 - (log10(l) * konst)) * base_color(1,2);
                color_tbl(l,2) = g;
                
                %b = b-b_div;
                b = (1 - (log10(l) * konst)) * base_color(1,3);
                color_tbl(l,3) = b;
            end
            
        end
                
        %%
        function obj = openSelector(save_stem, updatepaths)
            
            if nargin < 2
                updatepaths = false;
            else
                updatepaths = Force2Bool(updatepaths);
            end
            
            save_path = [save_stem 'spotAnnoObj'];
            save_ver = 1;
            mask_selection = [];
            pos_tbl = [];
            neg_tbl = [];
            ref_coord_tbl = [];
            bool3d = false;
            
            %load(save_path, 'spot_annotator');
            %obj = spot_annotator;
            %load(save_path, 'istructs', 'th_idx', 'th_tbl', 'pos_tbl', 'neg_tbl');
            %load(save_path, 'istructs', 'th_idx', 'th_tbl', 'pos_tbl', 'neg_tbl', 'tiff_path', 'tiff_channels', 'tiff_ch_selected', 'ref_coord_tbl', 'bool3d', 'lastz', 'maxz');
            %load(save_path, 'istructs', 'th_idx', 'th_tbl', 'pos_tbl', 'neg_tbl', 'ref_coord_tbl', 'bool3d', 'lastz', 'maxz', 'filimg_path', 'z_trim', 'mask_selection', 'save_ver'); 
            load(save_path);
            if save_ver > 1
                if save_ver >= 8
                    if ~refonly
                        load([save_path '_ptbl'], 'pos_tbl');
                        load([save_path '_ntbl'], 'neg_tbl');
                    end
                else
                    load([save_path '_ptbl'], 'pos_tbl');
                    load([save_path '_ntbl'], 'neg_tbl');
                end
            end
            
            if save_ver < 4
                %TODO
                if ~isempty(pos_tbl)
                    T = size(pos_tbl,1);
                    for t = 1:T
                        pos_tbl{t,1} = uint16(pos_tbl{t,1});
                    end
                end
            end
            
            if save_ver < 3
                %Do some cleanups
                %Force mask coords to ints
                if ~isempty(mask_selection)
                    temp_selmcoords = zeros(4,1);
                    temp_selmcoords(1,1) = uint16(mask_selection(1,1));
                    temp_selmcoords(2,1) = uint16(mask_selection(2,1));
                    temp_selmcoords(3,1) = uint16(mask_selection(3,1));
                    temp_selmcoords(4,1) = uint16(mask_selection(4,1));
                    mask_selection = temp_selmcoords;
                end
                
                if isfile([save_stem '_rnaspotsrun.mat'])
                    spotsrun = RNASpotsRun.loadFrom(save_stem);
                    maxz = spotsrun.idims_sample.z;
                end
                
                %Clean out 0,0,0 coordinates
                if ~isempty(pos_tbl)
                    T = size(pos_tbl,1);
                    for t = 1:T
                        subtbl = pos_tbl{t,1};
                        if ~isempty(subtbl)
                            good_rows = RNA_Threshold_SpotSelector.findNonzeroCoords(subtbl, bool3d);
                            if isempty(good_rows)
                                subtbl = [];
                            else
                                subtbl = subtbl(good_rows,:);
                            end
                            pos_tbl{t,1} = subtbl;
                        end
                    end
                end
                
                if ~isempty(neg_tbl)
                    T = size(neg_tbl,1);
                    for t = 1:T
                        subtbl = neg_tbl{t,1};
                        if ~isempty(subtbl)
                            good_rows = RNA_Threshold_SpotSelector.findNonzeroCoords(subtbl, bool3d);
                            if isempty(good_rows)
                                subtbl = [];
                            else
                                subtbl = subtbl(good_rows,:);
                            end
                            neg_tbl{t,1} = subtbl;
                        end
                    end
                end
                
                if ~isempty(ref_coord_tbl)
                    good_rows = RNA_Threshold_SpotSelector.findNonzeroCoords(ref_coord_tbl, bool3d);
                    if isempty(good_rows)
                        ref_coord_tbl = [];
                    else
                        ref_coord_tbl = ref_coord_tbl(good_rows,:);
                    end
                end
                
            end
            
            obj = RNA_Threshold_SpotSelector;
            obj.save_stem = save_stem;
            obj.img_structs = istructs;
            obj.threshold_idx = th_idx;
            obj.threshold_table = th_tbl;
            obj.positives = pos_tbl;
            obj.false_negs = neg_tbl;
            %obj.false_negs{obj.threshold_idx}
            %obj.tif_path = tiff_path;
            %obj.tif_chcount = tiff_channels;
            %obj.tif_ch = tiff_ch_selected;
            obj.mode_3d = bool3d;
            obj.current_slice = lastz;
            obj.max_slice = maxz;
            obj.imgdat_path = filimg_path;
            
            %TODO is there an alternative?
            if updatepaths
                obj.imgdat_path = [save_stem '_prefilteredIMG.mat'];
            end
            
            if isfile(obj.imgdat_path)
                load(obj.imgdat_path, 'img_filter');
                obj.loaded_ch = double(img_filter);
            else
                fprintf("Warning: Image data could not be found. GUI will be non-functional.\n");
                obj.loaded_ch = NaN(5,5,5);
            end
            
            refsetpath = [save_path '_refset.mat'];
            rfinfo = who('-file', refsetpath);
            if save_ver >= 7
                if isfile(refsetpath)
                    if save_ver >= 11
                        %Check if the timestamp is there.
                        %   Sometimes when using an old refset with a new
                        %   annoobj, the versions don't match.
                        if ~isempty(find(ismember(rfinfo, 'rtimestamp'),1))
                            load(refsetpath, 'ref_coord_tbl', 'rtimestamp');
                            obj.ref_last_modified = rtimestamp;
                        else
                            load(refsetpath, 'ref_coord_tbl');
                            obj.ref_last_modified = datetime;
                        end
                    else
                        load(refsetpath, 'ref_coord_tbl');
                        obj.ref_last_modified = datetime;
                    end
                else
                    ref_coord_tbl = [];
                    obj.ref_last_modified = datetime;
                end
            else
                obj.ref_last_modified = datetime;
            end
            if ~isempty(ref_coord_tbl)
                good_rows = RNA_Threshold_SpotSelector.findNonzeroCoords(ref_coord_tbl, bool3d);
                if isempty(good_rows)
                    ref_coord_tbl = [];
                else
                    ref_coord_tbl = ref_coord_tbl(good_rows,:);
                end
                obj.ref_last_modified = datetime;
                flag_fscores_dirty = true;
            end
            obj.ref_coords = ref_coord_tbl;
            
            s_img = istructs(1).image;
            obj.alloc_size = size(s_img, 1) * size(s_img,2);
            obj.n_alloc = 0;
            obj.n_used = 0;
            obj.n_table = [];
            obj.slice_drawn = 0;
            
            obj.ztrim = z_trim;
            obj.selmcoords = mask_selection;
            
            if save_ver >= 4
                obj.toggle_allz = toggle_az;
                obj.toggle_3dcount = toggle_3dc;
                obj.toggle_clr_local = toggle_cl;
                obj.toggle_singleSlice = toggle_ss;
                obj.toggle_del_unsnapped = toggle_du;
            else
                obj.toggle_allz = false;
                obj.toggle_3dcount = false;
                obj.toggle_clr_local = false;
                obj.toggle_singleSlice = false;
                obj.toggle_del_unsnapped = false;
            end
            
            if save_ver >= 5
                obj.toggle_cscale_max = toggle_cs;
            else
                obj.toggle_cscale_max = false;
            end
            
            if save_ver >= 6
                obj.f_scores = ftable;
            else
                t_count = size(obj.threshold_table,1);
                obj.f_scores = NaN(t_count,4);
            end
            
            if save_ver >= 8
                obj.z_min = zmin;
                obj.z_max = zmax;
            else
                obj.z_min = obj.ztrim + 1;
                obj.z_max = obj.max_slice - obj.ztrim;
            end
            
            if save_ver >= 9
                if ~isempty(find(ismember(rfinfo, 'rtimestamp'),1))
                    obj.f_scores_dirty = flag_fscores_dirty;
                else
                    obj.f_scores_dirty = true;
                end
            else
                obj.f_scores_dirty = true;
            end
            
            if save_ver >= 10
                obj.crosshair_color = crossclr;
            else
                obj.crosshair_color = [0.0, 0.0, 0.0];
            end
            
            if save_ver >= 12
                obj.toggle_raw_view = toggle_rvr;
            else
                obj.toggle_raw_view = false;
            end
 
            obj.toggle_change_zminmax = false;
            
            %obj.positives = pos_tbl(:,1);
            %obj.false_negs = neg_tbl(:,1);
            %obj.positives

            %DEBUG
            %t_count = size(obj.threshold_table,1);
            %obj.false_negs = cell(t_count);
            
        end
        
        %%
        function not_idxs = reverseIndices(vec_size, idxs)
            
            mtx = 1:1:vec_size;
            im = ismember(mtx,idxs); %Bool matrix - true means val is in original idx set
            %mtx = mtx(~im); 
            mtx = immultiply(mtx, ~im); % Set vals in orig set to 0
            [~,~,not_idxs] = find(mtx); %Get array of remaining nonzero values
            
        end
        
        %%
        function [src, targ] = moveCoordPairs(src, targ, src_idxs)
            
            %This is grossly inefficient since it requires a lot of memory
            %reallocation. If the lag ends up being a perpetual annoyance,
            %this should probably be updated to something better suited.
            %(Like storing coords in a linked list instead of a bloody array)
            
            %Move from src to targ
            diff_sz = size(src_idxs,2);
            t_sz = size(targ,1);
            new_t_sz = t_sz + diff_sz;
            t_new = NaN([new_t_sz 3]);
            t_new(1:t_sz,:) = targ(:,:);
            t_new(t_sz+1:new_t_sz,:) = src(src_idxs,:);
            
            %Delete from src
            s_sz = size(src,1);
            new_s_sz = s_sz - diff_sz;
            s_new = NaN([new_s_sz 3]);
            not_idxs = reverseIndices(s_sz, src_idxs);
            s_new(1:new_s_sz, :) = src(not_idxs,:);
            
            targ = t_new;
            src = s_new;
            
        end
        
        %%
        %TODO
        function idxs = findSelectedSpots(x_adj, y_adj, z, c_table, rad)
            %Update to only find those on current slice!
            
            %Any points within +- rad pixels of click (already scale
            %adjusted)
            
            %This is very clumsy, but it should work for now
            x_min = x_adj - rad;
            x_max = x_adj + rad;
            y_min = y_adj - rad;
            y_max = y_adj + rad;
            %fprintf("X Range: %f - %f\n", x_min, x_max);
            %fprintf("Y Range: %f - %f\n", y_min, y_max);
            
            x_bool = (c_table(:,1) >= x_min);
            x_bool = immultiply(x_bool, (c_table(:,1) <= x_max));
            
            y_bool = (c_table(:,2) >= y_min);
            y_bool = immultiply(y_bool, (c_table(:,2) <= y_max));
            
            z_bool = ones(size(c_table,1),1); %Defaults to pass all z
            if z > 0
                z_bool = (c_table(:,3) >= z);
            end
            
            net_bool = and(and(x_bool,y_bool), z_bool);
            
            %[x_rows,~,~] = find(x_bool);
            %[y_rows,~,~] = find(y_bool);
            %[z_rows,~,~] = find(z_bool);
            %fprintf("X matches: %d\n", size(x_rows,1));
            %fprintf("Y matches: %d\n", size(y_rows,1));
            
            %both_mtx = ismember(x_rows, y_rows);
            %[raw_idxs,~,~] = find(both_mtx);
            %fprintf("Total matches: %d\n", size(raw_idxs,1));
            %idxs = x_rows(raw_idxs);
 
            idxs = find(net_bool);
        end
        
        %%
        function bool = valueInMatrix(mtx, value)
            
            mtx2 = (mtx == value);
            m_max = max(max(mtx2, [], 2), [], 1);
            bool = (m_max > 0);
            
        end
        
        %%
        function xy_table = filterByZ(my_tbl, z_value)
            if isempty(my_tbl)
                xy_table = [];
                return;
            end
            row_match = (my_tbl(:,3) == z_value);
            row_idx = find(row_match);
            
            xy_table = my_tbl(row_idx, 1:2);
        end
        
        %%
        function keep_idx = findNonzeroCoords(my_tbl, check_z)
            
            bad_mtx = my_tbl(:,1) == 0;
            bad_mtx = bad_mtx | (my_tbl(:,2) == 0);
            if check_z
                bad_mtx = bad_mtx | (my_tbl(:,3) == 0);
            end
            [keep_idx,~,~] = find(~bad_mtx);
            
        end
        
        %%
        function [true_idx, false_idx] = isolateTFIndices(my_tbl)
            
            c_three = my_tbl(:,4);
            [true_idx,~,~] = find(c_three);
            [false_idx,~,~] = find(~c_three);
            
        end
       
        %%
        function [true_idx, false_idx, exclude_idx] = isolateTFEIndices(my_tbl)
            
            if isempty(my_tbl)
                true_idx = [];
                false_idx = [];
                exclude_idx = [];
                return;
            end
            
            c_three = my_tbl(:,4);
            %[pos_idx,~,~] = find(c_three);
            [false_idx,~,~] = find(~c_three);
            
            %Split the pos_idx rows into true(1) and exclude(2)
            t_1 = c_three == 1;
            t_2 = c_three == 2;
            [true_idx,~,~] = find(t_1);
            [exclude_idx,~,~] = find(t_2);
            
        end
        
        %%
        function ref_idx = findAReferenceIndex(th_tbl, target_idx)

            %Look down...
            ref_idx = 0;
            for i = 1:(target_idx-1)
                if th_tbl(i, 2)
                    ref_idx = i;
                end
            end
            
            %If nothing found, look up
            if ~ref_idx
                t_count = size(th_tbl,1);
                for i = (target_idx+1):t_count
                    if th_tbl(i, 2)
                        ref_idx = i;
                        break;
                    end
                end
            end
            
        end
        
        %%
        function radius = calculateClickRadius(fig_handle, fullX, fullY)
            
            rad_factor = 10;
            
            my_axes = findobj(fig_handle, 'type', 'axes');
            
            x_limits = get(my_axes,'XLim');
            y_limits = get(my_axes,'YLim');
            
            x_min = x_limits(1,1);
            x_max = x_limits(1,2);
            y_min = y_limits(1,1);
            y_max = y_limits(1,2);
            
            %fprintf("X Range: %f - %f\n", x_min, x_max);
            %fprintf("Y Range: %f - %f\n", y_min, y_max);
            
            x_range = x_max - x_min;
            y_range = y_max - y_min;
            
            x_scale = x_range./fullX;
            y_scale = y_range./fullY;
            
            %Use the larger one
            scale = x_scale;
            if y_scale > x_scale
                scale = y_scale;
            end
            
            
            radius = rad_factor .* scale;
            
            %fprintf("Scale Factor: %f\n", scale);
            %fprintf("Radius: %f\n", radius);
            
            
        end
        
        %%
        function [match_idxs, not_idxs] = findRowMatches(m_ref, m_other)
            
            %For differently sized matrices
            r_size = size(m_ref, 1);
            o_size = size(m_other, 1);
            cols = size(m_ref, 2);
            
            if o_size > r_size
                c_size = o_size;
                %Realloc ref
                m_temp = zeros([c_size cols]);
                m_temp(1:r_size,:) = m_ref(:,:);
                m_ref = m_temp;
            else
                c_size = r_size;
                %Realloc other
                m_temp = zeros([c_size cols]);
                m_temp(1:o_size,:) = m_other(:,:);
                m_other = m_temp;
            end
            
            %Do is member
            m_matcher = ismember(m_ref, m_other, 'rows');
            [idx_m,~,~] = find(m_matcher);
            [idx_n,~,~] = find(~m_matcher);
            
            %Dump any indices that are above the size of ref mtx
            copy_mtx = idx_m <= r_size;
            [copy_idx,~,~] = find(copy_mtx);
            match_idxs = idx_m(copy_idx,:);
            
            copy_mtx = idx_n <= r_size;
            [copy_idx,~,~] = find(copy_mtx);
            not_idxs = idx_n(copy_idx,:);
            
        end
        
        %%
        function tblout = maskFilterSpotTable(tblin, mask)

            if isempty(mask)
                tblout = tblin;
                return;
            end
            
            R = size(tblin,1);
            L = size(tblin,2);
            count = 0;
            tbltmp = zeros(R,L);
            %size(mask)
            
            for r = 1:R
                x = tblin(r,1);
                y = tblin(r,2);
                z = tblin(r,3);
                if x < 1 | y < 1 | z < 1
                    fprintf("Warning: Spot table contains invalid coordinate: %d,%d,%d\n",x,y,z);
                    continue;
                end
                if mask(y,x,z)
                    count = count+1;
                    tbltmp(count,:) = tblin(r,:);
                end
            end
            
            tblout = tbltmp(1:count,:);
        end
       
        %%
        function bool_res = selectorExists(save_stem)
            objpath = [save_stem 'spotAnnoObj.mat'];
            bool_res = isfile(objpath);
        end
        
        %%
        function bool_res = refsetExists(save_stem)
            bool_res = false;
            objpath = [save_stem 'spotAnnoObj_refset.mat'];
            if isfile(objpath)
                %Check if it is empty.
                load(objpath, 'ref_coord_tbl');
                bool_res = ~isempty(ref_coord_tbl);
            end
        end
        
        %%
        function fscores = loadFScores(save_stem, zmin, zmax, all_fields)
            if nargin < 4; all_fields = false; end
            fscores = [];
            
            if nargin < 3
                zmin = -1;
                zmax = -1;
            end
            
            %Actually, allow it to use manual selections as well...
            %if ~refsetExists(save_stem); return; end
            if ~RNA_Threshold_SpotSelector.selectorExists(save_stem); return; end
            selector = RNA_Threshold_SpotSelector.openSelector(save_stem, true);
            if zmin >= 1
                if selector.z_min ~= zmin
                    selector.z_min = zmin;
                    selector.f_scores_dirty = true;
                end
            end
            
            if zmax >= 1 & zmax <= selector.max_slice
                if selector.z_max ~= zmax
                    selector.z_max = zmax;
                    selector.f_scores_dirty = true;
                end
            end
            
            if selector.z_min > selector.z_max
                selector.z_min = selector.z_max;
                selector.f_scores_dirty = true;
            end
            
            if selector.f_scores_dirty
                selector = selector.updateFTable();
                selector = selector.saveMe();
            end
           
            if all_fields
                fscores = selector.f_scores;
            else
                fscores = selector.f_scores(:,1);
            end
            clear selector;
        end
        
        %%
        function timestamp = getRefsetTimestamp(save_stem)
            timestamp = [];
            if ~RNA_Threshold_SpotSelector.selectorExists(save_stem); return; end
            
            save_path = [save_stem 'spotAnnoObj'];
            load(save_path, 'save_ver');
            
            if save_ver >= 11
                refsetpath = [save_path '_refset.mat'];
                if isfile(refsetpath)
                    load(refsetpath, 'rtimestamp');
                    timestamp = rtimestamp;
                end
            end
            
        end
        
        %%
        function touchRefset(save_stem)
            if ~RNA_Threshold_SpotSelector.selectorExists(save_stem); return; end
            
            save_path = [save_stem 'spotAnnoObj'];
            load(save_path, 'save_ver');
            
            if save_ver >= 11
                refsetpath = [save_path '_refset.mat'];
                if isfile(refsetpath)
                    rtimestamp = datetime;
                    load(refsetpath, 'ref_coord_tbl');
                    save(refsetpath, 'ref_coord_tbl', 'rtimestamp');
                end
            end
        end
        
        %%
        function fixPosTable(save_stem)
            if ~RNA_Threshold_SpotSelector.selectorExists(save_stem); return; end
            
            %Check for pos tbl. Substitute dummy if file is empty.
            save_path = [save_stem 'spotAnnoObj'];
            pt_path = [save_path '_ptbl.mat'];
            finfo = who('-file', pt_path);
            if isempty(find(ismember(finfo, 'pos_tbl'),1))
                pos_tbl = [];
                save(pt_path, 'pos_tbl');
            else
                load(pt_path, 'pos_tbl');
                if ~isempty(pos_tbl)
                    fprintf('Pos table is intact. To replace callset, use loadNewSpotset.\n');
                    return;
                end
            end
            
            coord_path = [save_stem '_coordTable.mat'];
            spot_table_path = [save_stem '_spotTable.mat'];
            
            load(coord_path, 'coord_table');
            load(spot_table_path, 'spot_table');
            
            myanno = RNA_Threshold_SpotSelector.openSelector(save_stem);
            myanno = myanno.loadNewSpotset(spot_table, coord_table);
            myanno.saveMe();
        end
        
    end
    
    
end