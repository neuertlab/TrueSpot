%
%%
% Tool Codes:
%   'homebrew'
%   'bigfish'
%   'bigfish_nor'
%   'rsfish'
%   'deepblink'
%   'all' - Not always applicable

%%
classdef ImageResults

    %TODO Import BigFISH gauss fits
    
    properties
        image_name = '';
        mask_bounds; %Coordinates of area of image actually evaluated
        
        threshold_value_hb = NaN;
        threshold_index_hb = NaN;
        threshold_results;
        
        threshold_value_bf = NaN;
        threshold_index_bf = NaN;
        zmin_bf = 1;
        zmax_bf = 1;

        threshold_value_bfnr = NaN;
        threshold_index_bfnr = NaN;
        zmin_bfnr = 1;
        zmax_bfnr = 1;
        
        %From our tool - either fscore peak, or auto detected
        %This is to mark which th callset was imported
        threshold_value_rs = NaN;
        threshold_index_rs = NaN;
        
        %Vectors of SpotCalls
        truthset; %Table of true spots.
        callset_homebrew;
        callset_bigfish;
        callset_bigfish_nr;
        callset_rsfish;
        callset_deepblink; % Calls above .95 prob
        
        %Cells of tables for each truthset
        %Tables of threshold values vs.
        %   Th Value
        %   Spot count
        %   Sensitivity
        %   Precision
        %   F-Score
        res_homebrew;
        res_bigfish;
        res_bigfish_nr;
        res_rsfish;
        res_deepblink;
        
        %Single cell counting stats
        %Table of ON prop, avg rna per cell, avg rna per ON cell,
        %std rna per cell, std rna per ON cell, (previous 4 repeated with
        %nuc)
        cell_rna_counts;
        on_rna_count = 2; %Min number of spots in a cell for it to be ON
    end
    
    methods

        function [obj, index] = addTruthset(obj, save_stem)
            index = 0;
            if ~RNA_Threshold_SpotSelector.refsetExists(save_stem)
                return;
            end
            spotanno = RNA_Threshold_SpotSelector.openSelector(save_stem); 
            ref_tbl = spotanno.ref_coords;
            clear spotanno;
            
            spot_count = size(ref_tbl,1);
            ref_spots = ImageResults.initializeTrueCallTable(spot_count);
            ref_spots(:,1:3) = array2table(ref_tbl(:,1:3));
            
            if isempty(obj.truthset)
                obj.truthset = cell(1,1);
                obj.truthset{1} = ref_spots;
                index = 1;
            else
                existing_ts = size(obj.truthset,2);
                index = existing_ts + 1;
                obj.truthset{index} = ref_spots; %I assume MATLAB takes care of reallocation this way?
            end
            
            obj.res_homebrew{index} = [];
            obj.res_bigfish{index} = [];
            obj.res_bigfish_nr{index} = [];
            obj.res_rsfish{index} = [];
            obj.res_deepblink{index} = [];
        end
        
        function obj = setTruthset(obj, save_stem, ts_index)
            if ~RNA_Threshold_SpotSelector.refsetExists(save_stem)
                return;
            end
            spotanno = RNA_Threshold_SpotSelector.openSelector(save_stem); 
            ref_tbl = spotanno.ref_coords;
            clear spotanno;
            
            spot_count = size(ref_tbl,1);
            ref_spots = ImageResults.initializeTrueCallTable(spot_count);
            ref_spots(:,1:3) = array2table(ref_tbl(:,1:3));
            
            if isempty(obj.truthset)
                obj.truthset = cell(1,ts_index);
            end
            
            obj.truthset{ts_index} = ref_spots;
        end
        
        function [obj, success] = importHBResults(obj, save_stem, fixed_th, truthset_index)
            if nargin < 3
                fixed_th = 0;
            end
            if nargin < 4
                truthset_index = 1;
            end

            %Open up the spotsrun. If not there, return false.
            spotsrun = RNASpotsRun.loadFrom(save_stem);
            if isempty(spotsrun)
                success = false;
                return;
            end
            spotsrun.out_stem = save_stem;
            spotsrun.saveMe();

            %Update mask dims from spotsrun (overwrite with spotanno later, if it
            %is present)
            obj.mask_bounds.z0 = spotsrun.z_min_apply;
            obj.mask_bounds.z1 = spotsrun.z_max_apply;
            obj.mask_bounds.x0 = 1;
            obj.mask_bounds.x1 = spotsrun.idims_sample.x;
            obj.mask_bounds.y0 = 1;
            obj.mask_bounds.y1 = spotsrun.idims_sample.y;

            %Import thresholding results
            if ~isempty(spotsrun.threshold_results)
                obj.threshold_results = spotsrun.threshold_results;
                if fixed_th > 0
                    obj.threshold_value_hb = fixed_th;
                else
                    obj.threshold_value_hb = spotsrun.threshold_results.threshold;
                end
                obj.threshold_index_hb = RNAUtils.findThresholdIndex(obj.threshold_value_hb, transpose(spotsrun.threshold_results.x));
            else
                obj.threshold_results = struct.empty();
                if fixed_th > 0
                    obj.threshold_value_hb = fixed_th;
                else
                    obj.threshold_value_hb = spotsrun.intensity_threshold;
                end
                obj.threshold_index_hb = spotsrun.intensity_threshold - spotsrun.t_min + 1;
            end
            
            if obj.threshold_index_hb < 1
                fprintf('ERROR: HB run failed to threshold properly! Import failed.\n');
                success = false;
                return;
            end

            %Import callset at picked threshold
            [~, coord_table] = spotsrun.loadCoordinateTable();
            th_coords = coord_table{obj.threshold_index_hb, 1};
            clear coord_table;
            scount = size(th_coords,1);
            callset = ImageResults.initializeSpotCallTable(scount);
            callset(:,1:3) = array2table(th_coords(:,1:3));

            %Import spot counts
            [~, spot_table] = spotsrun.loadSpotsTable();
            T = size(spot_table,1);
            res_tbl = NaN(T,8);
            res_tbl(:,1:2) = spot_table(:,1:2);

            %Check for annoobj
            %If there, update fscores etc.
            if RNA_Threshold_SpotSelector.refsetExists(save_stem)
                spotanno = RNA_Threshold_SpotSelector.openSelector(save_stem, true);
                obj.mask_bounds.z0 = spotanno.z_min;
                obj.mask_bounds.z1 = spotanno.z_max;
                if ~isempty(spotanno.selmcoords)
                    obj.mask_bounds.x0 = spotanno.selmcoords(1,1);
                    obj.mask_bounds.x1 = spotanno.selmcoords(2,1);
                    obj.mask_bounds.y0 = spotanno.selmcoords(3,1);
                    obj.mask_bounds.y1 = spotanno.selmcoords(4,1);
                end
                spotanno.f_scores_dirty = true;
                spotanno = spotanno.updateFTable();
                spotanno.saveMe();
                res_tbl(:,5) = spotanno.f_scores(:,1);
                res_tbl(:,6) = spotanno.f_scores(:,2);
                res_tbl(:,7) = spotanno.f_scores(:,3);
                res_tbl(:,8) = spotanno.f_scores(:,4);
                res_tbl(:,3) = spotanno.f_scores(:,2) ./ (spotanno.f_scores(:,2) + spotanno.f_scores(:,4));
                res_tbl(:,4) = spotanno.f_scores(:,2) ./ (spotanno.f_scores(:,2) + spotanno.f_scores(:,3));

                %Match spots to truthset spots?
                thpos = spotanno.positives{obj.threshold_index_hb,1};
                callset(:,'tfcall') = array2table(~thpos(:,4));
%                 for i = 1:scount
%                     if thpos(i,4) == 1
%                         callset(i,'tfcall') = table(0);
%                         callset(i,'tfcall_desc') = table("true pos");
%                     else
%                         callset(i,'tfcall') = table(1);
%                         callset(i,'tfcall_desc') = table("false pos");
%                     end
%                 end
            end

            %Save res_tbl to obj
            varnames = ImageResults.getResTableVarNames();
            if (truthset_index > 1) | iscell(obj.res_homebrew)
                obj.res_homebrew{truthset_index} = array2table(res_tbl, 'VariableNames', varnames);
            else
                obj.res_homebrew = array2table(res_tbl, 'VariableNames', varnames);
            end

            %Check for quant. Import if present.
            [hb_dir, ~, ~] = fileparts(save_stem);
            dir_contents = dir(hb_dir);
            fcount = size(dir_contents,1);
            match_idx = 0;
            for i = 1:fcount
                if endsWith(dir_contents(i,1).name, '_quantData.mat')
                    match_idx = i;
                    break;
                end
            end
            if match_idx > 0
                %Quant found, load
                load([hb_dir filesep dir_contents(match_idx,1).name], 'quant_results');
%                 qcells = quant_results.cell_rna_data;
%                 cell_count = size(qcells,2);
%                 for c = 1:cell_count
%                     this_cell = qcells(c);
%                     if isempty(this_cell.spots); continue; end
%                     cell_spot_count = size(this_cell.spots, 2);
%                     cell_x = this_cell.cell_loc.left;
%                     cell_y = this_cell.cell_loc.top;
%                     cell_z = this_cell.cell_loc.z_bottom;
%                     for j = 1:cell_spot_count
%                         this_spot = this_cell.spots(j);
%                         spot_x = this_spot.x + cell_x;
%                         spot_y = this_spot.y + cell_y;
%                         spot_z = this_spot.z + cell_z;
%                         [~, match_idx] = ImageResults.findMatchingCall(callset, spot_x, spot_y, spot_z);
%                         if match_idx > 0
%                             callset(match_idx) = callset(match_idx).importFromQuantSpot(this_spot, cell_x, cell_y, cell_z);
%                         end
%                     end
%                 end
            end

            obj.callset_homebrew = callset;

            success = true;
        end
        
        function [obj, success] = importBFResults(obj, save_stem, is_rescaled, truthset_index)
            if nargin < 3; is_rescaled = true; end
            if nargin < 4; truthset_index = 1; end
            success = false;
            
            %Check for presence of run
            [bf_dir, bf_pfx, ~] = fileparts(save_stem);
            summary_path = [bf_dir filesep 'summary.txt'];
            if ~isfile(summary_path)
                fprintf('BigFISH run could not be found in %s!\n', bf_dir);
                return;
            end
            
            %Import threshold values first
            [zmin, zmax, bfthresh] = BigfishCompare.readSummaryTxt(summary_path);
            if is_rescaled
                obj.zmin_bf = zmin;
                obj.zmax_bf = zmax;
                obj.threshold_value_bf = bfthresh;
            else
                obj.zmin_bfnr = zmin;
                obj.zmax_bfnr = zmax;
                obj.threshold_value_bfnr = bfthresh;
            end
            
            %See if need to import...
            spot_table_path = [bf_dir filesep bf_pfx '_spotTable.mat'];
            coord_table_path = [bf_dir filesep bf_pfx '_coordTable.mat'];
            if ~isfile(coord_table_path)
                [coord_table, spot_table] = BigfishCompare.importBigFishCsvs(bf_dir, save_stem, zmin);
            else
                load(coord_table_path, 'coord_table');
                load(spot_table_path, 'spot_table');
            end
            
            %Get threshold index
            x_tbl = transpose(spot_table(:,1));
            thresh_idx = RNAUtils.findThresholdIndex(bfthresh, x_tbl);
            if is_rescaled
                obj.threshold_index_bf = thresh_idx;
                obj.threshold_value_bf = spot_table(thresh_idx,1);
            else
                obj.threshold_index_bfnr = thresh_idx;
                obj.threshold_value_bfnr = spot_table(thresh_idx,1);
            end
            
            %Import callset
            th_coords = coord_table{thresh_idx,1};
            while isempty(th_coords)
                if thresh_idx < 2
                    fprintf('BigFISH callset for %s is empty!\n', bf_dir);
                    return;
                end
                thresh_idx = thresh_idx - 1;
                th_coords = coord_table{thresh_idx,1};
                
                if is_rescaled
                    obj.threshold_index_bf = thresh_idx;
                else
                    obj.threshold_index_bfnr = thresh_idx;
                end
            end
            
            clear coord_table;
            spot_count = size(th_coords,1);
            callset = ImageResults.initializeSpotCallTable(spot_count);
            callset(:,1:3) = array2table(th_coords(:,1:3));
            
            %Import fit table, if present
            fit_table_path = [bf_dir filesep 'fitspots.csv'];
            if isfile(fit_table_path)
                fit_table = csvread(fit_table_path); %z,y,x, 0 based coords?
                fit_count = size(fit_table,1);
                
                if fit_count ~= spot_count
                    %Manually match :)
                    cs = 1;
                    crem = spot_count;
                    frem = fit_count;
                    for s = 1:fit_count
                        fprintf('WARNING: Fit table is not the same size as call table! Will attempt to match...\n');
                        
                        if cs > spot_count
                            fprintf('End of call table hit. Stopping fit matching...\n');
                            break;
                        end
                        
                        fx = fit_table(s,3) + 1;
                        fy = fit_table(s,2) + 1;
                        fz = fit_table(s,1) + 1;
                        
                        if fit_count > spot_count
                            if frem <= crem
                                %Just copy.
                                callset{cs,'fit_x'} = fx;
                                callset{cs,'fit_y'} = fy;
                                callset{cs,'fit_z'} = fz;
                                frem = frem - 1;
                                crem = crem - 1;
                                cs = cs + 1;
                                continue;
                            end
                        else
                            if frem >= crem
                                callset{cs,'fit_x'} = fx;
                                callset{cs,'fit_y'} = fy;
                                callset{cs,'fit_z'} = fz;
                                frem = frem - 1;
                                crem = crem - 1;
                                cs = cs + 1;
                                continue;
                            end
                        end
                        
                        smatch = false;
                        while ~smatch
                            cx = callset{cs,'init_x'};
                            cy = callset{cs,'init_y'};
                            cz = callset{cs,'init_z'};
                            
                            dxq = (fx - cx) ^ 2;
                            dyq = (fy - cy) ^ 2;
                            dzq = (fz - cz) ^ 2;
                            dist3 = sqrt(dxq + dyq + dzq);
                            
                            if dist3 <= 4
                                smatch = true;
                                callset{cs,'fit_x'} = fx;
                                callset{cs,'fit_y'} = fy;
                                callset{cs,'fit_z'} = fz;
                            else
                                if fit_count > spot_count
                                    if frem <= crem
                                        callset{cs,'fit_x'} = fx;
                                        callset{cs,'fit_y'} = fy;
                                        callset{cs,'fit_z'} = fz;
                                        smatch = true;
                                    end
                                else
                                    if frem >= crem
                                        callset{cs,'fit_x'} = fx;
                                        callset{cs,'fit_y'} = fy;
                                        callset{cs,'fit_z'} = fz;
                                        smatch = true;
                                    end
                                end
                            end
                            
                            crem = crem - 1;
                            cs = cs + 1;
                        end
                        
                        frem = frem - 1;
                    end
                else
                    for s = 1:fit_count
                        callset{s,'fit_x'} = fit_table(s,3) + 1;
                        callset{s,'fit_y'} = fit_table(s,2) + 1;
                        callset{s,'fit_z'} = fit_table(s,1) + 1;
                    end
                end
            end
            
            %Start restbl
            T = size(spot_table,1);
            res_tbl = NaN(T,8);
            res_tbl(:,1:2) = spot_table(:,1:2);
            
            %Check for spotanno, import if present
            if RNA_Threshold_SpotSelector.refsetExists(save_stem)
                spotanno = RNA_Threshold_SpotSelector.openSelector(save_stem, true);
                spotanno.f_scores_dirty = true;
                spotanno = spotanno.updateFTable();
                res_tbl(:,5) = spotanno.f_scores(:,1);
                res_tbl(:,6) = spotanno.f_scores(:,2);
                res_tbl(:,7) = spotanno.f_scores(:,3);
                res_tbl(:,8) = spotanno.f_scores(:,4);
                res_tbl(:,3) = spotanno.f_scores(:,2) ./ (spotanno.f_scores(:,2) + spotanno.f_scores(:,4));
                res_tbl(:,4) = spotanno.f_scores(:,2) ./ (spotanno.f_scores(:,2) + spotanno.f_scores(:,3));

                %Match spots to truthset spots?
                thpos = spotanno.positives{thresh_idx,1};
                callset(:,'tfcall') = array2table(~thpos(:,4));
%                 for i = 1:scount
%                     if thpos(i,4) == 1
%                         spot_list(i).tfcall = 0;
%                     else
%                         spot_list(i).tfcall = 1;
%                     end
%                 end
            end
            
            %Save
            varnames = ImageResults.getResTableVarNames();
            res_tbl_table = array2table(res_tbl, 'VariableNames', varnames);
            if is_rescaled
                obj.callset_bigfish = callset;
                if (truthset_index > 1) | iscell(obj.res_bigfish)
                    obj.res_bigfish{truthset_index} = res_tbl_table;
                else
                    obj.res_bigfish = res_tbl_table;
                end
            else
                obj.callset_bigfish_nr = callset;
                if (truthset_index > 1) | iscell(obj.res_bigfish_nr)
                    obj.res_bigfish_nr{truthset_index} = res_tbl_table;
                else
                    obj.res_bigfish_nr = res_tbl_table;
                end
            end
            
            
            success = true;
        end
        
        function [obj, success] = importRSResults(obj, save_stem, th_interval, truthset_index)
            %Most of the runs I had it delete the fit files. Maybe don't.
            if nargin < 3
                th_interval = 0.1/250;
            end
            if nargin < 4
                truthset_index = 1;
            end
            success = false;
            
            %Check for run. Return if not found.
            coord_table_path = [save_stem '_coordTable.mat'];
            if ~isfile(coord_table_path)
                fprintf('Coord table for RS run %s not found. Check to see if import is required...\n', save_stem);
                return;
            end
            
            %Load spot table and coord table
            load([save_stem '_spotTable'], 'spot_table');
            
            %Fix tables
            if nnz(isnan(spot_table(:,1))) > 0
                %Trim NaN rows, coord tables, and fit tables
                fit_table_old = [];
                fit_table_path = [save_stem '_fitTable.mat'];
                if isfile(fit_table_path)
                    load(fit_table_path, 'fit_table');
                    fit_table_old = fit_table;
                end
                
                load(coord_table_path, 'coord_table');
                
                realloc = nnz(~isnan(spot_table(:,1)));
                spot_table_old = spot_table;
                coord_table_old = coord_table;
                spot_table = NaN(realloc, 2);
                coord_table = cell(realloc,1);
                fit_table = [];
                if ~isempty(fit_table_old); fit_table = cell(realloc,1); end
                T = size(spot_table_old,1);
                
                j = 1;
                for i = 1:T
                    if isnan(spot_table_old(i,1)); continue; end
                    spot_table(j,:) = spot_table_old(i,:);
                    coord_table{j,1} = coord_table_old{i,1};
                    if ~isempty(fit_table_old)
                        fit_table{j,1} = fit_table_old{i,1};
                    end
                    j = j+1;
                end
                
                
                %Save
                writerver = 23022100;
                writer_ver_str = 'v 23.02.21.0*';
                save([save_stem '_coordTable.mat'], 'coord_table', 'writerver', 'writer_ver_str');
                save([save_stem '_spotTable.mat'], 'spot_table', 'writerver', 'writer_ver_str');
                save([save_stem '_fitTable.mat'], 'fit_table', 'writerver', 'writer_ver_str');
                
                clear coord_table;
                clear coord_table_old;
                clear spot_table_old;
                if ~isempty(fit_table)
                    clear fit_table;
                    clear fit_table_old;
                end
            end
            
            %Update spot table threshold values
            spot_table_int = spot_table;
            spot_table(:,1) = spot_table(:,1) .* th_interval;
            
            %Initialize restbl
            T = size(spot_table,1);
            res_tbl = NaN(T,8);
            res_tbl(:,1:2) = spot_table(:,1:2);
            
            %Check for spotsanno.
            th_idx = 0;
            call_values = [];
            if RNA_Threshold_SpotSelector.refsetExists(save_stem)
                spotanno = RNA_Threshold_SpotSelector.openSelector(save_stem, true);
                spotanno = spotanno.updateFTable();
                res_tbl(:,5) = spotanno.f_scores(:,1);
                res_tbl(:,6) = spotanno.f_scores(:,2);
                res_tbl(:,7) = spotanno.f_scores(:,3);
                res_tbl(:,8) = spotanno.f_scores(:,4);
                res_tbl(:,3) = spotanno.f_scores(:,2) ./ (spotanno.f_scores(:,2) + spotanno.f_scores(:,4));
                res_tbl(:,4) = spotanno.f_scores(:,2) ./ (spotanno.f_scores(:,2) + spotanno.f_scores(:,3));
                
                [~, th_idx] = max(res_tbl(:,5));
                postbl = spotanno.positives{th_idx,1};
                clear spotanno;
                call_values = transpose(postbl(:,4));
                clear postbl;
            end
            
            %Threshold if no fscore curve available
            if th_idx < 1
                thresh_res = RNAThreshold.runWithPreset(spot_table_int, [], 3);
                if isempty(thresh_res)
                    fprintf('Thresholding failed. Using lowest value...\n');
                    th_idx = 1;
                else
                    th_idx = thresh_res.threshold - spot_table_int(1,1) + 1;
                end
            end
            if isnan(th_idx)
                %Still having issues with the thresholder malfunctioning
                %Look into this AGAIN when have more time.
                fprintf('Thresholding failed. Using lowest value...\n');
                th_idx = 1;
            end
            obj.threshold_index_rs = th_idx;
            obj.threshold_value_rs = spot_table(th_idx,1);
            
             %Import callset
            load(coord_table_path, 'coord_table');
            th_coords = coord_table{th_idx,1};
            clear coord_table;
            spot_count = size(th_coords,1);
            callset = ImageResults.initializeSpotCallTable(spot_count);
            callset(:,1:3) = array2table(th_coords(:,1:3));
            if ~isempty(call_values)
                callset(:,'tfcall') = array2table(~call_values(:,1));
            end

            %Check for fit data and import if present
            fit_table_path = [save_stem '_fitTable.mat'];
            if isfile(fit_table_path)
                load(fit_table_path, 'fit_table');
                th_fits = fit_table{th_idx,1};
                callset(:,'fit_x') = array2table(th_fits(:,1));
                callset(:,'fit_y') = array2table(th_fits(:,2));
                callset(:,'fit_z') = array2table(th_fits(:,3));
                callset(:,'fit_peak') = array2table(th_fits(:,4));
            end
            
            %Save callset and restable
            varnames = ImageResults.getResTableVarNames();
            if (truthset_index > 1) | iscell(obj.res_rsfish)
                obj.res_rsfish{truthset_index} = array2table(res_tbl, 'VariableNames', varnames);
            else
                obj.res_rsfish = array2table(res_tbl, 'VariableNames', varnames);
            end
            
            if (truthset_index > 1) | iscell(obj.callset_rsfish)
                obj.callset_rsfish{truthset_index} = callset;
            else
                obj.callset_rsfish = callset;
            end
            
            success = true;
        end
        
        function [obj, success] = importDBResults(obj, save_stem, truthset_index)
            if nargin < 3
                truthset_index = 1;
            end
            success = false;
            coord_table_path = [save_stem '_coordTable.mat'];
            if ~isfile(coord_table_path)
                fprintf('Coord table for DeepBlink run %s not found. Check to see if import is required...\n', save_stem);
                return;
            end
            
            spot_table_path = [save_stem '_spotTable.mat'];
            if isfile(spot_table_path)
                load(spot_table_path, 'spot_table');
                prob_cutoffs = transpose(spot_table(:,1));
                prob_count = size(prob_cutoffs,2);
                th_idx = RNAUtils.findThresholdIndex(0.95, prob_cutoffs);
            else
                prob_cutoffs = [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];
                prob_count = size(prob_cutoffs,2);
                th_idx = prob_count - 1;
            end

            %Import callset
            load(coord_table_path, 'coord_table');
            th_coords = coord_table{th_idx,1};
            %Slide down until there is one that is not empty...
            while isempty(th_coords)
                if th_idx < 2
                    fprintf('Coord table for DeepBlink run %s is empty! Check to see if import is required...\n', save_stem);
                    return;
                end
                th_idx = th_idx - 1;
                th_coords = coord_table{th_idx,1};
            end
            
            spot_count = size(th_coords,1);
            callset = ImageResults.initializeSpotCallTable(spot_count);
            callset(:,1:3) = array2table(th_coords(:,1:3));

            %Init restbl
            res_tbl = NaN(prob_count,8);
            res_tbl(:,1) = prob_cutoffs(:);
            for i = 1:prob_count
                res_tbl(i,2) = size(coord_table{i,1}, 1);
            end
            clear coord_table;
            
            %Check for spotsanno.
            if RNA_Threshold_SpotSelector.refsetExists(save_stem)
                spotanno = RNA_Threshold_SpotSelector.openSelector(save_stem, true);
                spotanno = spotanno.updateFTable();
                res_tbl(:,5) = spotanno.f_scores(:,1);
                res_tbl(:,6) = spotanno.f_scores(:,2);
                res_tbl(:,7) = spotanno.f_scores(:,3);
                res_tbl(:,8) = spotanno.f_scores(:,4);
                res_tbl(:,3) = spotanno.f_scores(:,2) ./ (spotanno.f_scores(:,2) + spotanno.f_scores(:,4));
                res_tbl(:,4) = spotanno.f_scores(:,2) ./ (spotanno.f_scores(:,2) + spotanno.f_scores(:,3));
                
                postbl = spotanno.positives{th_idx,1};
                clear spotanno;
                
                callset(:,'tfcall') = array2table(~postbl(:,4));
            end

            %Check for original table
            [db_dir, ~, ~] = fileparts(save_stem);
            dir_contents = dir(db_dir);

            input_file = [];
            content_count = size(dir_contents,1);
            for i = 1:content_count
                if endsWith(dir_contents(i,1).name, '.csv')
                    input_file = [db_dir filesep dir_contents(i,1).name];
                    fprintf('Found csv in DeepBlink directory: "%s"...\n', input_file);
                    break;
                end
            end
            if ~isempty(input_file)
                import_table = readtable(input_file,'Delimiter',',','ReadVariableNames',true,'Format',...
                    '%f%f%f%f');
                import_mtx = table2array(import_table);
                import_keep_idxs = find(import_mtx(:,3) >= 0.95);
                
                if ~isempty(import_keep_idxs)
                    callset(:,'fit_x') = array2table(import_mtx(import_keep_idxs,1));
                    callset(:,'fit_y') = array2table(import_mtx(import_keep_idxs,2));
                    callset(:,'fit_z') = array2table(import_mtx(import_keep_idxs,4));
                end
            end

            %Save to obj
            varnames = ImageResults.getResTableVarNames();
            if (truthset_index > 1) | iscell(obj.res_deepblink)
                obj.res_deepblink{truthset_index} = array2table(res_tbl, 'VariableNames', varnames);
            else
                obj.res_deepblink = array2table(res_tbl, 'VariableNames', varnames);
            end
            
            if (truthset_index > 1) | iscell(obj.callset_deepblink)
                obj.callset_deepblink{truthset_index} = callset;
            else
                obj.callset_deepblink = callset;
            end
            
            success = true;
        end
        
        function [obj, success] = importCellRNACounts(obj, quant_file, toolcode)
            success = false;
            if ~isfile(quant_file)
                fprintf('Quant file %s does not exist!\n', quant_file);
                return;
            end

            %Get row correlating to tool code...
            row_idx = ImageResults.getToolQTableRowIndex(toolcode);
            if row_idx < 1
                fprintf('Tool code "%s" not recognized!\n', toolcode);
                return;
            end

            finfo = who('-file', quant_file);
            if ~isempty(find(ismember(finfo, 'quant_results'),1))
                load(quant_file, 'quant_results');
                qcells = quant_results.cell_rna_data;
                cell_count = size(qcells, 2);

                %Dump to vectors...
                cell_rna_count = NaN(1,cell_count);
                nuc_rna_count = NaN(1,cell_count);
                for c = 1:cell_count
                    this_cell = qcells(c);
                    if ~isempty(this_cell.spots)
                        scount = size(this_cell.spots,2);
                        cell_rna_count(c) = scount;
                        nuc_rna_count(c) = this_cell.spotcount_nuc;
                    else
                        cell_rna_count(c) = 0;
                        nuc_rna_count(c) = 0;
                    end
                end

                %Calculate stats...
                on_idxs = find(cell_rna_count >= obj.on_rna_count);
                obj.cell_rna_counts(row_idx, 1) = cell_count;
                obj.cell_rna_counts(row_idx, 3) = nanmean(cell_rna_count, 'all');
                obj.cell_rna_counts(row_idx, 4) = nanstd(cell_rna_count, 0, 'all');
                obj.cell_rna_counts(row_idx, 7) = nanmean(nuc_rna_count, 'all');
                obj.cell_rna_counts(row_idx, 8) = nanstd(nuc_rna_count, 0, 'all');

                if ~isempty(on_idxs)
                    on_count = size(on_idxs,1);
                    obj.cell_rna_counts(row_idx, 2) = on_count/cell_count;

                    on_cell_count = cell_rna_count(on_idxs);
                    obj.cell_rna_counts(row_idx, 5) = nanmean(on_cell_count, 'all');
                    obj.cell_rna_counts(row_idx, 6) = nanstd(on_cell_count, 0, 'all');

                    on_nuc_count = nuc_rna_count(on_idxs);
                    obj.cell_rna_counts(row_idx, 9) = nanmean(on_nuc_count, 'all');
                    obj.cell_rna_counts(row_idx, 10) = nanstd(on_nuc_count, 0, 'all');
                else
                    obj.cell_rna_counts(row_idx, 2) = 0.0;
                    obj.cell_rna_counts(row_idx, 5) = 0.0;
                    obj.cell_rna_counts(row_idx, 6) = 0.0;
                    obj.cell_rna_counts(row_idx, 9) = 0.0;
                    obj.cell_rna_counts(row_idx, 10) = 0.0;
                end

            elseif ~isempty(find(ismember(finfo, 'cell_spot_counts'),1))
                %Quick quant results
                load(quant_file, 'cell_spot_counts');
                cell_count = size(cell_spot_counts,1);
                on_idxs = find(cell_spot_counts(:,2) >= obj.on_rna_count);

                obj.cell_rna_counts(row_idx, 1) = cell_count;
                obj.cell_rna_counts(row_idx, 3) = nanmean(cell_spot_counts(:,2), 'all');
                obj.cell_rna_counts(row_idx, 4) = nanstd(cell_spot_counts(:,2), 0, 'all');

                if ~isempty(on_idxs)
                    on_count = size(on_idxs,1);
                    obj.cell_rna_counts(row_idx, 2) = on_count/cell_count;

                    on_cell_count = cell_spot_counts(on_idxs,2);
                    obj.cell_rna_counts(row_idx, 5) = nanmean(on_cell_count, 'all');
                    obj.cell_rna_counts(row_idx, 6) = nanstd(on_cell_count, 0, 'all');
                else
                    obj.cell_rna_counts(row_idx, 2) = 0.0;
                    obj.cell_rna_counts(row_idx, 5) = 0.0;
                    obj.cell_rna_counts(row_idx, 6) = 0.0;
                end
            else
                fprintf('Valid quant table in %s not found!\n', quant_file);
                return;
            end

            success = true;
        end
        
        function fig_handle = renderSpotCountPlot(obj, tool_code, color, figno, existing_fig)
            if nargin < 4
                figno = round(rand() * 10000);
            end
            if nargin < 3
                color = [0 0 0];
            end
            if nargin < 5
                existing_fig = [];
            end

            %Get data
            thval = 0;
            res_tbl = []; fig_handle = [];
            if isempty(tool_code); return; end
            if strcmp(tool_code, 'homebrew')
                if isempty(obj.res_homebrew); return; end
                res_tbl = obj.res_homebrew{1};
                thval = obj.threshold_value_hb;
            elseif strcmp(tool_code, 'bigfish')
                if isempty(obj.res_bigfish); return; end
                res_tbl = obj.res_bigfish{1};
                thval = obj.threshold_value_bf;
            elseif strcmp(tool_code, 'bigfish_nor')
                if isempty(obj.res_bigfish_nr); return; end
                res_tbl = obj.res_bigfish_nr{1};
                thval = obj.threshold_value_bfnr;
            elseif strcmp(tool_code, 'rsfish')
                if isempty(obj.res_rsfish); return; end
                res_tbl = obj.res_rsfish{1};
            elseif strcmp(tool_code, 'deepblink')
                if isempty(obj.res_deepblink); return; end
                res_tbl = obj.res_deepblink{1};
            end
            if isempty(res_tbl); return; end

            color_light = color + ((1.0 - color) .* 0.5);

            if isempty(existing_fig)
                fig_handle = figure(figno);
                clf;
            else
                fig_handle = existing_fig;
                figure(fig_handle);
            end
            if isempty(res_tbl); return; end

            x = table2array(res_tbl(:,1));
            y = log10(double(table2array(res_tbl(:,2))));
            plot(x,y,'LineWidth',2,'Color',color);
            hold on;

            if thval > 0
                xline(thval, '--', 'Threshold', 'Color', color_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end

            if strcmp(tool_code, 'homebrew')
                %Draw threshold range
                if ~isempty(obj.threshold_results)
                    tlist = RNAThreshold.getAllThresholdSuggestions(obj.threshold_results);
                    if ~isempty(tlist)
                        tavg = nanmean(tlist, 'all');
                        tstdev = nanstd(tlist, 0, 'all');
                        xline(tavg + tstdev, ':', 'Color', color_light,'LineWidth',1);
                        xline(tavg - tstdev, ':', 'Color', color_light,'LineWidth',1);
                    end
                end
            end

            set(gca,'XTickLabel',[]);
            xlabel('Threshold (a.u.)');
            ylabel('log10(# Spots Detected)');
        end
        
        function fig_handle = renderFScorePlot(obj, tool_code, color, truthset_index, figno, existing_fig)
            if nargin < 5
                figno = round(rand() * 10000);
            end
            if nargin < 3
                color = [0 0 0];
            end
            if nargin < 4
                truthset_index = 1;
            end
            if nargin < 6
                existing_fig = [];
            end

            %Get data
            thval = 0;
            res_tbl = []; fig_handle = [];
            if isempty(tool_code); return; end
            if strcmp(tool_code, 'homebrew')
                if isempty(obj.res_homebrew); return; end
                res_tbl = obj.res_homebrew{truthset_index};
                thval = obj.threshold_value_hb;
            elseif strcmp(tool_code, 'bigfish')
                if isempty(obj.res_bigfish); return; end
                res_tbl = obj.res_bigfish{truthset_index};
                thval = obj.threshold_value_bf;
            elseif strcmp(tool_code, 'bigfish_nor')
                if isempty(obj.res_bigfish_nr); return; end
                res_tbl = obj.res_bigfish_nr{truthset_index};
                thval = obj.threshold_value_bfnr;
            elseif strcmp(tool_code, 'rsfish')
                if isempty(obj.res_rsfish); return; end
                res_tbl = obj.res_rsfish{truthset_index};
            elseif strcmp(tool_code, 'deepblink')
                if isempty(obj.res_deepblink); return; end
                res_tbl = obj.res_deepblink{truthset_index};
            end

            if isempty(res_tbl); return; end

            color_light = color + ((1.0 - color) .* 0.5);

            if isempty(existing_fig)
                fig_handle = figure(figno);
                clf;
            else
                fig_handle = existing_fig;
                figure(fig_handle);
            end

            x = table2array(res_tbl(:,1));
            y = table2array(res_tbl(:,5));
            plot(x,y,'LineWidth',2,'Color',color);
            hold on;

            if thval > 0
                xline(thval, '--', 'Threshold', 'Color', color_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end

            if strcmp(tool_code, 'homebrew')
                %Draw threshold range
                if ~isempty(obj.threshold_results)
                    tlist = RNAThreshold.getAllThresholdSuggestions(obj.threshold_results);
                    if ~isempty(tlist)
                        tavg = nanmean(tlist, 'all');
                        tstdev = nanstd(tlist, 0, 'all');
                        xline(tavg + tstdev, ':', 'Color', color_light,'LineWidth',1);
                        xline(tavg - tstdev, ':', 'Color', color_light,'LineWidth',1);
                    end
                end
            end

            ylim([0.0 1.0]);
            xticks([]);
            xlabel('Threshold (a.u.)');
            ylabel('F-Score');
        end
        
        function fig_handle = renderROCPlot(obj, tool_code, color, truthset_index, figno, existing_fig)
            if nargin < 5
                figno = round(rand() * 10000);
            end
            if nargin < 3
                color = [0 0 0];
            end
            if nargin < 4
                truthset_index = 1;
            end
            if nargin < 6
                existing_fig = [];
            end

            %Get data
            thidx = 0;
            res_tbl = []; fig_handle = [];
            if isempty(tool_code); return; end
            if strcmp(tool_code, 'homebrew')
                if isempty(obj.res_homebrew); return; end
                res_tbl = obj.res_homebrew{truthset_index};
                thidx = obj.threshold_index_hb;
            elseif strcmp(tool_code, 'bigfish')
                if isempty(obj.res_bigfish); return; end
                res_tbl = obj.res_bigfish{truthset_index};
                thidx = obj.threshold_index_bf;
            elseif strcmp(tool_code, 'bigfish_nor')
                if isempty(obj.res_bigfish_nr); return; end
                res_tbl = obj.res_bigfish_nr{truthset_index};
                thidx = obj.threshold_index_bfnr;
            elseif strcmp(tool_code, 'rsfish')
                if isempty(obj.res_rsfish); return; end
                res_tbl = obj.res_rsfish{truthset_index};
            elseif strcmp(tool_code, 'deepblink')
                if isempty(obj.res_deepblink); return; end
                res_tbl = obj.res_deepblink{truthset_index};
            end
            if isempty(res_tbl); return; end

            color_dark = max(color - ((1.0 - color) .* 0.5),[0,0,0]);

            if isempty(existing_fig)
                fig_handle = figure(figno);
                clf;
            else
                fig_handle = existing_fig;
                figure(fig_handle);
            end

%             ax = axes();
%             ax.XLim = [0.0 1.0];
%             ax.YLim = [0.0 1.0];
%             ax.XDir = 'reverse';
%             hold on;

            x = table2array(res_tbl(:,3));
            y = table2array(res_tbl(:,4));
            plot(x, y, 'LineStyle', 'none', 'Marker', '.', 'MarkerEdgeColor', color);
            hold on;

            if thidx > 0
                plot(x(thidx), y(thidx), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 15, 'MarkerEdgeColor', color_dark);
            end

            ylim([0.0 1.0]);
            xlim([0.0 1.0]);
            ax = gca;
            ax.XDir = 'reverse';

            xlabel('Sensitivity');
            ylabel('Precision');
            %axes(ax);
        end
        
        function [obj, success] = updateCallsetTruthsetCompares(obj)
            %TODO
        end

        function obj = clearResultsForTool(obj, tool_code)
            if isempty(tool_code); return; end
            if strcmp(tool_code, 'homebrew') | strcmp(tool_code, 'all')
                obj.threshold_value_hb = NaN;
                obj.threshold_index_hb = NaN;
                obj.callset_homebrew = [];
                tscount = size(obj.res_homebrew,2);
                for i = 1:tscount
                    obj.res_homebrew{i} = [];
                end
                obj.threshold_results = [];
            end
            if strcmp(tool_code, 'bigfish') | strcmp(tool_code, 'all')
                obj.threshold_value_bf = NaN;
                obj.threshold_index_bf = NaN;
                obj.callset_bigfish = [];
                tscount = size(obj.res_bigfish,2);
                for i = 1:tscount
                    obj.res_bigfish{i} = [];
                end
                obj.zmin_bf = 1;
                obj.zmax_bf = 1;
            end
            if strcmp(tool_code, 'bigfish_nor') | strcmp(tool_code, 'all')
                obj.threshold_value_bfnr = NaN;
                obj.threshold_index_bfnr = NaN;
                obj.callset_bigfish_nr = [];
                tscount = size(obj.res_bigfish_nr,2);
                for i = 1:tscount
                    obj.res_bigfish_nr{i} = [];
                end
                obj.zmin_bfnr = 1;
                obj.zmax_bfnr = 1;
            end
            if strcmp(tool_code, 'rsfish') | strcmp(tool_code, 'all')
                obj.threshold_value_rs = NaN;
                obj.threshold_index_rs = NaN;
                obj.callset_rsfish = [];
                tscount = size(obj.res_rsfish,2);
                for i = 1:tscount
                    obj.res_rsfish{i} = [];
                end
            end
            if strcmp(tool_code, 'deepblink') | strcmp(tool_code, 'all')
                obj.callset_deepblink = [];
                tscount = size(obj.res_deepblink,2);
                for i = 1:tscount
                    obj.res_deepblink{i} = [];
                end
            end
        end
        
    end
    
    methods (Static)
        
        function row = getToolQTableRowIndex(tool_code)
            row = 0;
            if isempty(tool_code); return; end
            if strcmp(tool_code, 'homebrew'); row = 1; end
            if strcmp(tool_code, 'bigfish'); row = 2; end
            if strcmp(tool_code, 'bigfish_nor'); row = 3; end
            if strcmp(tool_code, 'rsfish'); row = 4; end
            if strcmp(tool_code, 'deepblink'); row = 5; end
        end

        function res_table_var_names = getResTableVarNames()
            res_table_var_names = {'thresholdValue' 'spotCount' 'sensitivity' 'precision' 'fScore'...
                'true_pos', 'false_pos', 'false_neg'};
        end
        
        function call_table_var_names = getCallTableVarNames()
            call_table_var_names = {'init_x' 'init_y' 'init_z' 'init_peak' 'init_total'...
                'fit_x' 'fit_y' 'fit_z' 'fit_peak' 'fit_total'...
                'fit_xFWHM' 'fit_yFWHM' 'tfcall'...
                'x_nearest_true' 'y_nearest_true' 'z_nearest_true' 'dist_nearest_true'};
        end
        
        function res_table = initializeResTable(th_count)
            varTypes = {'double' 'uint32' 'double' 'double' 'double', 'uint32', 'uint32', 'uint32'};
            varNames = ImageResults.getResTableVarNames();
            res_table = table('Size', [th_count, 5], 'VariableTypes',varTypes, 'VariableNames',varNames);
        end
        
        function call_table = initializeTrueCallTable(alloc)
            varNames = {'isnap_x' 'isnap_y' 'isnap_z' 'isnap_peak' 'isnap_total'...
                'exact_x' 'exact_y' 'exact_z' 'exact_peak' 'exact_total'...
                'fit_xFWHM' 'fit_yFWHM'};
            varTypes = {'uint16' 'uint16' 'uint16' 'double' 'double'...
                'double' 'double' 'double' 'double' 'double'...
                'double' 'double'};
            table_size = [alloc size(varNames,2)];
            call_table = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);
            nanvec = NaN(alloc,1);
            nantbl = array2table(nanvec);
            call_table(:,'isnap_peak') = nantbl;
            call_table(:,'isnap_total') = nantbl;
            call_table(:,'exact_x') = nantbl;
            call_table(:,'exact_y') = nantbl;
            call_table(:,'exact_z') = nantbl;
            call_table(:,'exact_peak') = nantbl;
            call_table(:,'exact_total') = nantbl;
            call_table(:,'fit_xFWHM') = nantbl;
            call_table(:,'fit_yFWHM') = nantbl;
        end
        
        function call_table = initializeSpotCallTable(alloc)
            varNames = ImageResults.getCallTableVarNames();
            varTypes = {'uint16' 'uint16' 'uint16' 'double' 'double'...
                'double' 'double' 'double' 'double' 'double'...
                'double' 'double' 'uint8' ...
                'double' 'double' 'double' 'double'};
            table_size = [alloc size(varNames,2)];
            call_table = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);
            nanvec = NaN(alloc,1);
            nantbl = array2table(nanvec);
            call_table(:,'init_peak') = nantbl;
            call_table(:,'init_total') = nantbl;
            call_table(:,'fit_x') = nantbl;
            call_table(:,'fit_y') = nantbl;
            call_table(:,'fit_z') = nantbl;
            call_table(:,'fit_peak') = nantbl;
            call_table(:,'fit_total') = nantbl;
            call_table(:,'fit_xFWHM') = nantbl;
            call_table(:,'fit_yFWHM') = nantbl;
            call_table(:,'x_nearest_true') = nantbl;
            call_table(:,'y_nearest_true') = nantbl;
            call_table(:,'z_nearest_true') = nantbl;
            call_table(:,'dist_nearest_true') = nantbl;
            
            call_table(:,'tfcall') = array2table(repmat(2,alloc,1));
            
%             tfcall_vec = repmat("none", alloc, 1);
%             tfcall_tbl = array2table(tfcall_vec);
%             call_table(:,'tfcall_desc') = tfcall_tbl;
        end
        
        function image_results = initializeNew()
            image_results = ImageResults;
            image_results.threshold_results = struct.empty();
            image_results.mask_bounds = struct('x0', 1, 'y0', 1, 'z0', 1, 'x1', 1, 'y1', 1, 'z1', 1);
            
            image_results.truthset = cell.empty();
            image_results.callset_homebrew = table.empty();
            image_results.callset_bigfish = table.empty();
            image_results.callset_bigfish_nr = table.empty();
            image_results.callset_rsfish = table.empty();
            image_results.callset_deepblink = table.empty();
            
            image_results.res_homebrew = cell.empty();
            image_results.res_bigfish = cell.empty();
            image_results.res_bigfish_nr = cell.empty();
            image_results.res_rsfish = cell.empty();
            image_results.res_deepblink = cell.empty();
            
            varTypes = {'uint16' 'double' 'double' 'double' 'double'...
                'double' 'double' 'double' 'double' 'double'};
            varNames = {'cellCount' 'onProp' 'meanRNACell'...
                'stdevRNACell' 'meanRNAOnCell' 'stdevRNAOnCell'...
                'meanRNANuc' 'stdevRNANuc'...
                'meanRNAOnNuc' 'stdevRNAOnNuc'};
            rowNames = {'Homebrew' 'BigFISH' 'BigFISH_NR'...
                'RSFISH' 'DeepBlink'};
            tblsize = [size(rowNames,2) size(varNames,2)];
            image_results.cell_rna_counts = table('Size',tblsize, ...
                'VariableTypes',varTypes,...
                'VariableNames',varNames,...
                'RowNames',rowNames);
        end
        
        function [spot_call, index] = findMatchingCall(calls, x, y, z)
            spot_call = SpotCall.empty();
            index = 0;
            call_count = size(calls,2);
            for i = 1:call_count
                this_spot = calls(i);
                if(this_spot.init_x ~= x); continue; end
                if(this_spot.init_y ~= y); continue; end
                if(this_spot.init_z ~= z); continue; end
                spot_call = this_spot;
                index = i;
                break;
            end
        end

    end
end