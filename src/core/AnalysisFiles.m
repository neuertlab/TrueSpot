classdef AnalysisFiles
   
    methods (Static)

        %%
        function spot_table = analysisCallset2SpotcountTable(analysis, toolid)
            spot_table = [];
            rfieldname = ['results_' toolid];
            if ~isfield(analysis, rfieldname); return; end
            if ~isfield(analysis.(rfieldname), 'callset'); return; end

            spot_table = AnalysisFiles.callset2SpotcountTable(analysis.(rfieldname).callset);
        end

        %%
        function spot_table = callset2SpotcountTable(call_table)
            %If all found th are integers, try to fill in gaps.
            th_vals = unique(call_table{:, 'dropout_thresh'});
            th_vals = th_vals(th_vals > 0);

            %Guess interval by finding smallest adjacent gap.
            mingap = min(diff(th_vals), [], 'all', 'omitnan');
            min_val = min(th_vals, [], 'all', 'omitnan');
            max_val = max(th_vals, [], 'all', 'omitnan');

            vals = min_val:mingap:max_val;
            T = size(vals, 2);
            spot_table = NaN(T,2);

            spot_table(:,1) = vals;
            for t = 1:T
                spot_table(t,2) = nnz(call_table{:, 'dropout_thresh'} >= spot_table(t,1));
            end
        end

        %%
        function [analysis, okay] = addExpRefSet(analysis, ref_coords, ref_region, name, snaprad_3, snaprad_z, snap_minth)
            if nargin < 5; snaprad_3 = 4; end
            if nargin < 6; snaprad_z = 2; end
            if nargin < 7; snap_minth = 0.00001; end

            okay = false;
            if isempty(analysis); return; end
            if isempty(ref_coords); return; end
            if isempty(name); return; end

            %Add to main analysis
            if ~isfield(analysis, 'refsets')
                analysis.refsets = struct();
            end
            analysis.refsets.(name) = struct('name', name);
            analysis.refsets.(name).exprefset = ref_coords;
            if ~isempty(ref_region)
                analysis.refsets.(name).truthset_region = ref_region;
            else
                analysis.refsets.(name).truthset_region = struct();
                analysis.refsets.(name).truthset_region.x0 = 1;
                analysis.refsets.(name).truthset_region.x1 = analysis.image_dims.x;
                analysis.refsets.(name).truthset_region.y0 = 1;
                analysis.refsets.(name).truthset_region.y1 = analysis.image_dims.y;
                analysis.refsets.(name).truthset_region.z0 = 1;
                analysis.refsets.(name).truthset_region.z1 = analysis.image_dims.z;
            end
            analysis.refsets.(name).timestamp = datetime;

            %Add to results structs
            allfields = fieldnames(analysis);
            fieldcount = size(allfields, 1);
            for i = 1:fieldcount
                myfield = allfields{i,1};
                if startsWith(myfield, 'results_')
                    rstruct = analysis.(myfield);
                    if ~isfield(rstruct, 'callset'); continue; end
                    [~, ref_call_map] = RNACoords.updateTFCalls(...
                        rstruct.callset, ref_coords, snaprad_3, snaprad_z, snap_minth);

                    if ~isfield(rstruct, 'benchmarks')
                        rstruct.benchmarks = struct();
                    end

                    pstruct = struct();
                    pstruct.name = name;
                    pstruct.ref_call_map = ref_call_map;
                    rstruct.benchmarks.(name) = pstruct; %Activation will add the rest

                    analysis.(myfield) = rstruct;
                end
            end

            [analysis, ~] = AnalysisFiles.activateExpRefSet(analysis, name);
            okay = true;
        end

        %%
        function restable = initializeResTable(th_count)
            varNames = {'thresholdValue' 'spotCount' 'sensitivity' 'precision' 'fScore' 'true_pos' 'false_pos' 'false_neg'};
            varTypes = {'single' 'uint32' 'double' 'double' 'double' 'uint32' 'uint32' 'uint32'};

            table_size = [th_count size(varNames,2)];
            restable = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);
        end

        %%
        function pmetrics = calculatePerformanceMetrics(call_table, th_val, pmetrics, emptyTSBehavior, minTh, maxTh)
            %This method ASSUMES that call_table flag columns are all up to
            %date!

            if nargin < 4
                emptyTSBehavior = -1; %Set NaN
            end
            if nargin < 5; minTh = 0; end
            if nargin < 6; maxTh = 0; end

            if isempty(pmetrics)
                pmetrics = struct();
            end

            vec_istrimmed = call_table{:,'is_trimmed_out'};
            vec_intsreg = call_table{:,'in_truth_region'};
            vec_isreal = call_table{:,'is_true'};
            vec_dropth = call_table{:,'dropout_thresh'};

            vec_real = vec_isreal & vec_intsreg & ~vec_istrimmed;
            empty_ts = nnz(vec_real) <= 0;

            any_trimmed = nnz(vec_istrimmed) > 0;
            spot_table = RNAUtils.spotTableFromCallTable(call_table, true, minTh, maxTh);
            %spot_table = AnalysisFiles.callset2SpotcountTable(call_table);

            th_count = size(spot_table,1);
            res_untrimmed = AnalysisFiles.initializeResTable(th_count);
            res_trimmed = table.empty();

            res_untrimmed{:,'thresholdValue'} = double(spot_table(:,1));

            if any_trimmed
                res_trimmed = AnalysisFiles.initializeResTable(th_count);
                res_trimmed{:,'thresholdValue'} = double(spot_table(:,1));
            else
                %Clean trimmed struct if it is present
                if isfield(pmetrics, 'performance_trimmed')
                    pmetrics = rmfield(pmetrics, 'performance_trimmed');
                end
                if isfield(pmetrics, 'pr_auc_trimmed')
                    pmetrics = rmfield(pmetrics, 'pr_auc_trimmed');
                end
                if isfield(pmetrics, 'fscore_peak_trimmed')
                    pmetrics = rmfield(pmetrics, 'fscore_peak_trimmed');
                end
                if isfield(pmetrics, 'fscore_autoth_trimmed')
                    pmetrics = rmfield(pmetrics, 'fscore_autoth_trimmed');
                end
            end

            sc_all = NaN(th_count,2);
            tp_all = NaN(th_count,2);
            fp_all = NaN(th_count,2);
            fn_all = NaN(th_count,2);
            for t = 1:th_count
                th = spot_table(t,1);
                pos_vec = (vec_dropth >= th) & vec_intsreg;
                tp_vec = pos_vec & vec_isreal;
                fp_vec = pos_vec & ~vec_isreal;
                fn_vec = (vec_dropth < th) & vec_intsreg & vec_isreal;

                tp_all(t,1) = nnz(tp_vec);
                fp_all(t,1) = nnz(fp_vec);
                fn_all(t,1) = nnz(fn_vec);
                sc_all(t,1) = nnz(pos_vec);

                %Repeat for trimmed, if applicable
                if any_trimmed
                    pos_vec = (vec_dropth >= th) & vec_intsreg & ~vec_istrimmed;
                    neg_vec = (vec_dropth < th) & vec_intsreg & ~vec_istrimmed;
                    tp_vec = pos_vec & vec_isreal;
                    fp_vec = pos_vec & ~vec_isreal;
                    fn_vec = neg_vec & vec_isreal;

                    tp_all(t,2) = nnz(tp_vec);
                    fp_all(t,2) = nnz(fp_vec);
                    fn_all(t,2) = nnz(fn_vec);
                    sc_all(t,2) = nnz(pos_vec);
                end
            end

            %Let's speed up the easy calculations...
            res_untrimmed{:, 'spotCount'} = uint32(sc_all(:,1));
            res_untrimmed{:, 'true_pos'} = uint32(tp_all(:,1));
            res_untrimmed{:, 'false_pos'} = uint32(fp_all(:,1));
            res_untrimmed{:, 'false_neg'} = uint32(fn_all(:,1));

            %NaN-ify when total denominator is 0.
            denom = tp_all(:,1) + fn_all(:,1);
            denom(denom == 0) = NaN;
            recall = tp_all(:,1) ./ denom;
            if empty_ts
                if emptyTSBehavior == 0
                    recall(:) = 0;
                elseif emptyTSBehavior == 1
                    recall(:) = 1;
                end
            end

            denom = tp_all(:,1) + fp_all(:,1);
            denom(denom == 0) = NaN;
            precision = tp_all(:,1) ./ denom;

            fscores = (2 .* precision .* recall) ./ (precision + recall);
            pr_auc = RNAUtils.calculateAUC(recall, precision);
            peak_fscore = max(fscores, [], 'all', 'omitnan');
            peak_recall = max(recall, [], 'all', 'omitnan');
            res_untrimmed{:, 'sensitivity'} = recall;
            res_untrimmed{:, 'precision'} = precision;
            res_untrimmed{:, 'fScore'} = fscores;
            if any_trimmed
                res_trimmed{:, 'spotCount'} = uint32(sc_all(:,2));
                res_trimmed{:, 'true_pos'} = uint32(tp_all(:,2));
                res_trimmed{:, 'false_pos'} = uint32(fp_all(:,2));
                res_trimmed{:, 'false_neg'} = uint32(fn_all(:,2));

                denom = tp_all(:,2) + fn_all(:,2);
                denom(denom == 0) = NaN;
                recall = tp_all(:,2) ./ denom;
                if empty_ts
                    if emptyTSBehavior == 0
                        recall(:) = 0;
                    elseif emptyTSBehavior == 1
                        recall(:) = 1;
                    end
                end

                denom = tp_all(:,2) + fp_all(:,2);
                denom(denom == 0) = NaN;
                precision = tp_all(:,2) ./ denom;

                fscores = (2 .* precision .* recall) ./ (precision + recall);
                pr_auc_trim = RNAUtils.calculateAUC(recall, precision);
                peak_fscore_trim = max(fscores, [], 'all', 'omitnan');
                peak_recall_trim = max(recall, [], 'all', 'omitnan');
                res_trimmed{:, 'sensitivity'} = recall;
                res_trimmed{:, 'precision'} = precision;
                res_trimmed{:, 'fScore'} = fscores;
            end

            th_idx = 0;
            if th_val > 0
                th_idx = RNAUtils.findThresholdIndex(th_val, spot_table(:,1).');
            end

            %Save to output struct
            pmetrics.performance = res_untrimmed;
            pmetrics.pr_auc = pr_auc;
            pmetrics.fscore_peak = peak_fscore;
            pmetrics.max_recall = peak_recall;
            if th_idx > 0
                pmetrics.fscore_autoth = res_untrimmed{th_idx, 'fScore'};
            end
            if any_trimmed
                pmetrics.performance_trimmed = res_trimmed;
                pmetrics.pr_auc_trimmed = pr_auc_trim;
                pmetrics.fscore_peak_trimmed = peak_fscore_trim;
                pmetrics.max_recall_trimmed = peak_recall_trim;
                if th_idx > 0
                    pmetrics.fscore_autoth_trimmed = res_trimmed{th_idx, 'fScore'};
                end
            end

        end

        %%
        function [rstruct, okay] = activateExpRefSetForTool(analysis, rstruct, refset, name)
            okay = false;
            if isempty(rstruct); return; end
            if isempty(analysis); return; end
            if ~isfield(rstruct, 'callset'); return; end

            [rstruct, okay] = AnalysisFiles.updateCallsetTrimRes(rstruct, analysis.image_dims);
            if ~okay; return; end
            
            callset = rstruct.callset;
            callset{:, 'is_true'} = 0;
            callset{:, 'in_truth_region'} = 1;
            callset{:, 'xdist_ref'} = NaN;
            callset{:, 'ydist_ref'} = NaN;
            callset{:, 'zdist_ref'} = NaN;
            callset{:, 'xydist_ref'} = NaN;
            callset{:, 'xyzdist_ref'} = NaN;
            if ~isempty(refset)
                %Check for an existing ref/call match set
                ref_call_map = [];
                if isfield(rstruct, 'benchmarks')
                    if isfield(rstruct.benchmarks, name)
                        ref_call_map = rstruct.benchmarks.(name).ref_call_map;
                    end
                end

                if isempty(ref_call_map)
                    %Generate.
                    [callset, ref_call_map] = RNACoords.updateTFCalls(...
                        callset, refset.exprefset, 4, 2, 0.0001);
                end

                pstruct = struct();
                pstruct.name = name;
                pstruct.ref_call_map = ref_call_map;

                %Update callset table.
                ref_matched_rows = find(ref_call_map > 0);
                ref_call_map_i = ref_call_map(ref_matched_rows);
                refset_matched = refset.exprefset(ref_matched_rows, :);
                callset{ref_call_map_i, 'is_true'} = 1;

                %False negatives
                ref_rows_fn = find(ref_call_map == 0);
                if ~isempty(ref_rows_fn)
                    idims = analysis.image_dims;
                    sz = [idims.y idims.x idims.z];
                    xx = refset.exprefset(ref_rows_fn,1);
                    yy = refset.exprefset(ref_rows_fn,2);
                    zz = refset.exprefset(ref_rows_fn,3);
                    fn_1d = sub2ind(sz, yy, xx, zz);

                    fn_match = ismember(callset{:, 'coord_1d'}, fn_1d);
                    if ~isempty(fn_match)
                        callset{fn_match, 'is_true'} = 1;
                    end

                    clear idims sz xx yy zz fn_match fn_1d
                end

                if isfield(refset, 'truthset_region')
                    rows = find(callset{:, 'isnap_x'} < refset.truthset_region.x0);
                    if ~isempty(rows); callset{rows, 'in_truth_region'} = 0; end
                    rows = find(callset{:, 'isnap_x'} > refset.truthset_region.x1);
                    if ~isempty(rows); callset{rows, 'in_truth_region'} = 0; end

                    rows = find(callset{:, 'isnap_y'} < refset.truthset_region.y0);
                    if ~isempty(rows); callset{rows, 'in_truth_region'} = 0; end
                    rows = find(callset{:, 'isnap_y'} > refset.truthset_region.y1);
                    if ~isempty(rows); callset{rows, 'in_truth_region'} = 0; end

                    rows = find(callset{:, 'isnap_z'} < refset.truthset_region.z0);
                    if ~isempty(rows); callset{rows, 'in_truth_region'} = 0; end
                    rows = find(callset{:, 'isnap_z'} > refset.truthset_region.z1);
                    if ~isempty(rows); callset{rows, 'in_truth_region'} = 0; end
                end

                callset = RNACoords.updateRefDistances(callset, refset_matched, ref_call_map_i);

                thval = 0;
                if isfield(rstruct, 'threshold'); thval = rstruct.threshold; end
                pstruct = AnalysisFiles.calculatePerformanceMetrics(callset, thval, pstruct);
                pstruct.timestamp = datetime;
                if ~isfield(rstruct, 'benchmarks')
                    rstruct.benchmarks = struct();
                end
                rstruct.benchmarks.(name) = pstruct;
            end

            rstruct.callset = callset;
            rstruct.timestamp = datetime;
        end

        %%
        function [analysis, okay] = activateExpRefSet(analysis, name)
            %If name is empty, then deactivate.
            %Check to make sure requested refset exists...
            okay = false;
            if ~isfield(analysis, 'refsets'); return; end

            %Check to make sure requested refset exists... (If not, then do nothing)
            if ~isempty(name)
                if ~isfield(analysis.refsets, name); return; end
                myref = analysis.refsets.(name);
            else
                myref = [];
            end

            allfields = fieldnames(analysis);
            fieldcount = size(allfields, 1);
            for i = 1:fieldcount
                myfield = allfields{i,1};
                if startsWith(myfield, 'results_') 
                    %toolid = replace(myfield, 'results_', '');
                    %[analysis, ~] = AnalysisFiles.updateCallsetTrim(analysis, toolid);
                    [analysis.(myfield), ~] = AnalysisFiles.activateExpRefSetForTool(analysis, analysis.(myfield), myref, name);
                end
            end

            if ~isempty(myref)
                analysis.activeRefId = name;
            else
                analysis.activeRefId = '';
            end

            okay = true;
        end

        %%
        function analysis = fixExpRefsetOrganization(analysis)
            if ~isfield(analysis, 'exprefset'); return; end

            aremove = {'exprefset' 'truthset_region' 'last_ts'};
            rremove = {'performance' 'pr_auc' 'fscore_peak' 'fscore_autoth'...
                'ref_call_map' 'performance_trimmed' 'pr_auc_trimmed'...
                'fscore_peak_trimmed' 'fscore_autoth_trimmed' 'last_ts'};

            analysis.refsets = struct();
            analysis.activeRefId = '';

            allfields = fieldnames(analysis);
            fieldcount = size(allfields, 1);
            for i = 1:fieldcount
                myfield = allfields{i,1};
                if startsWith(myfield, 'truthset_')
                    if ~strcmp(myfield, 'truthset_region')
                        %Move.
                        rsname = replace(myfield, 'truthset_', '');
                        analysis.refsets.(rsname) = analysis.(myfield);
                        analysis.refsets.(rsname).name = rsname;
                        analysis = rmfield(analysis, myfield);
                    end
                elseif startsWith(myfield, 'results_') 
                    %toolid = replace(myfield, 'results_', '');
                    rstruct = analysis.(myfield);
                    rstruct.benchmarks = struct();

                    rfields = fieldnames(rstruct);
                    J = size(rfields,1);
                    for j = 1:J
                        rfield = rfields{j,1};
                        if startsWith(rfield, 'truthset_')
                            rsname = replace(rfield, 'truthset_', '');
                            rstruct.benchmarks.(rsname) = rstruct.(rfield);
                            rstruct.benchmarks.(rsname).name = rsname;
                            rstruct = rmfield(rstruct, rfield);
                        end
                    end

                    %Remove columns in callset table.
                    tvars = rstruct.callset.Properties.VariableNames;
                    V = size(tvars, 2);
                    for v = 1:V
                        if startsWith(tvars{v}, 'is_true_')
                            rstruct.callset = removevars(rstruct.callset, tvars{v});
                        elseif startsWith(tvars{v}, 'in_truth_region_')
                            rstruct.callset = removevars(rstruct.callset, tvars{v});
                        end
                    end

                    R = size(rremove, 2);
                    for r = 1:R
                        if isfield(rstruct, rremove{r})
                            rstruct = rmfield(rstruct, rremove{r});
                        end
                    end
                    rstruct.timestamp = datetime;

                    analysis.(myfield) = rstruct;
                end
            end

            R = size(aremove, 2);
            for r = 1:R
                if isfield(analysis, aremove{r})
                    analysis = rmfield(analysis, aremove{r});
                end
            end

            %Activate truthset BH if it exists.
            if isfield(analysis.refsets, 'BH')
                [analysis, ~] = AnalysisFiles.activateExpRefSet(analysis, 'BH');
            else
                [analysis, ~] = AnalysisFiles.activateExpRefSet(analysis, []);
            end
        end

        %%
        function [rstruct, okay] = updateCallsetTrimRes(rstruct, idims, ignoreZ)

            if nargin < 3; ignoreZ = false; end

            okay = false;
            if isempty(rstruct); return; end
            if isempty(idims); return; end

            xmin = 1; ymin = 1; zmin = 1;
            xmax = idims.x;
            ymax = idims.y;
            zmax = idims.z;

            if isfield(rstruct, 'x_min'); xmin = rstruct.x_min; end
            if isfield(rstruct, 'x_max'); xmax = rstruct.x_max; end
            if isfield(rstruct, 'y_min'); ymin = rstruct.y_min; end
            if isfield(rstruct, 'y_max'); ymax = rstruct.y_max; end
            if isfield(rstruct, 'z_min'); zmin = rstruct.z_min; end
            if isfield(rstruct, 'z_max'); zmax = rstruct.z_max; end

            calls_full = rstruct.callset;
            x_okay = and(calls_full{:, 'isnap_x'} >= xmin, calls_full{:, 'isnap_x'} <= xmax);
            y_okay = and(calls_full{:, 'isnap_y'} >= ymin, calls_full{:, 'isnap_y'} <= ymax);

            if ignoreZ
                z_okay = true(size(calls_full, 1), 1);
            else
                z_okay = and(calls_full{:, 'isnap_z'} >= zmin, calls_full{:, 'isnap_z'} <= zmax);
            end
            keep_rows = find(x_okay & y_okay & z_okay);

            if ~isempty(keep_rows)
                rstruct.callset{:,'is_trimmed_out'} = 1;
                rstruct.callset{keep_rows,'is_trimmed_out'} = 0;
            else
                rstruct.callset{:,'is_trimmed_out'} = 1;
            end

            okay = true;
        end

        %%
        function [callset, okay] = applyTruthRegionMask(truthdims, callset, idims, ignoreZ)

            if nargin < 4; ignoreZ = false; end

            okay = false;
            if isempty(callset); return; end
            if isempty(idims); return; end

            xmin = 1; ymin = 1; zmin = 1;
            xmax = idims.x;
            ymax = idims.y;
            zmax = idims.z;

            if ~isempty(truthdims)
                if isfield(truthdims, 'x0'); xmin = truthdims.x0; end
                if isfield(truthdims, 'x1'); xmax = truthdims.x1; end
                if isfield(truthdims, 'y0'); ymin = truthdims.y0; end
                if isfield(truthdims, 'y1'); ymax = truthdims.y1; end
                if isfield(truthdims, 'z0'); zmin = truthdims.z0; end
                if isfield(truthdims, 'z1'); zmax = truthdims.z1; end
            end

            x_okay = and(callset{:, 'isnap_x'} >= xmin, callset{:, 'isnap_x'} <= xmax);
            y_okay = and(callset{:, 'isnap_y'} >= ymin, callset{:, 'isnap_y'} <= ymax);
            if ignoreZ
                z_okay = true(size(callset, 1), 1);
            else
                z_okay = and(callset{:, 'isnap_z'} >= zmin, callset{:, 'isnap_z'} <= zmax);
            end
            keep_rows = find(x_okay & y_okay & z_okay);

            if ~isempty(keep_rows)
                callset{:,'in_truth_region'} = 0;
                callset{keep_rows,'in_truth_region'} = 1;
            else
                callset{:,'in_truth_region'} = 0;
            end

            okay = true;
        end

        %%
        function [analysis, okay] = updateCallsetTrim(analysis, toolid)
            okay = false;
            resfieldname = ['results_' toolid];
            if ~isfield(analysis, resfieldname); return; end

            idims = analysis.image_dims;
            rstruct = analysis.(resfieldname);
            
            [rstruct, okay] = AnalysisFiles.updateCallsetTrimRes(rstruct, idims);

            analysis.(resfieldname) = rstruct;
        end
        
        %%
        function analysis = rethresholdExp(analysis, preset, verbose, updatePerf)
            if nargin < 4; updatePerf = false; end
            if ~isfield(analysis, 'results_hb'); return; end
            rstruct = analysis.results_hb;

            %Regen spot table.
            if verbose; fprintf('>> Extracting spot count table...\n'); end
            spot_table = RNAUtils.spotTableFromCallTable(rstruct.callset, false, rstruct.th_scan_min, rstruct.th_scan_max);

            if verbose; fprintf('>> Rethresholding with preset %d...\n', preset); end
            th_res = RNAThreshold.runWithPreset(spot_table, [], preset);
            th_res.lowNoiseFlag = RNAThreshold.checkLowNoise(analysis.image_dims, spot_table);
            rstruct.threshold_details = th_res;
            if th_res.lowNoiseFlag
                rstruct.threshold = 1;
            else
                rstruct.threshold = th_res.threshold;
            end

            if updatePerf
                if isfield(rstruct, 'benchmarks')
                    refsetnames = fieldnames(rstruct.benchmarks);
                    setcount = size(refsetnames, 1);
                    for i = 1:setcount
                        setname = refsetnames{i};
                        if verbose; fprintf('>> Updating performance calculations for refset %s...\n', setname); end
                        bstruct = rstruct.benchmarks.(setname);
                        bstruct = AnalysisFiles.calculatePerformanceMetrics(...
                            rstruct.callset, rstruct.threshold, bstruct, -1, rstruct.th_scan_min, rstruct.th_scan_max);
                        rstruct.benchmarks.(setname) = bstruct;
                    end
                end
            else
                thidx = 0;
                if th_res.threshold > 0
                    thidx = RNAUtils.findThresholdIndex(th_res.threshold, spot_table(:,1)');
                end

                if isfield(rstruct, 'benchmarks')
                    refsetnames = fieldnames(rstruct.benchmarks);
                    setcount = size(refsetnames, 1);
                    for i = 1:setcount
                        setname = refsetnames{i};
                        if verbose; fprintf('>> Updating threshold FScores for refset %s...\n', setname); end
                        bstruct = rstruct.benchmarks.(setname);
                        if isfield(bstruct, 'performance_trimmed')
                            if thidx > 0
                                bstruct.fscore_autoth_trimmed = bstruct.performance_trimmed{thidx, 'fScore'};
                            else
                                bstruct.fscore_autoth_trimmed = NaN;
                            end
                        end
                        if isfield(bstruct, 'performance')
                            if thidx > 0
                                bstruct.fscore_autoth = bstruct.performance{thidx, 'fScore'};
                            else
                                bstruct.fscore_autoth = NaN;
                            end
                        end
                        rstruct.benchmarks.(setname) = bstruct;
                    end
                end
            end

            analysis.results_hb = rstruct; 
        end

        %%
        function analysis = rethresholdSim(analysis, preset, verbose, updatePerf)
            if nargin < 4; updatePerf = false; end
            if ~isfield(analysis, 'results_hb'); return; end
            rstruct = analysis.results_hb;

            %Regen spot table.
            if verbose; fprintf('>> Extracting spot count table...\n'); end
            spot_table = RNAUtils.spotTableFromCallTable(rstruct.callset, false, rstruct.th_scan_min, rstruct.th_scan_max);

            if verbose; fprintf('>> Rethresholding with preset %d...\n', preset); end
            th_res = RNAThreshold.runWithPreset(spot_table, [], preset);
            th_res.lowNoiseFlag = RNAThreshold.checkLowNoise(analysis.image_dims, spot_table);
            rstruct.threshold_details = th_res;
            if th_res.lowNoiseFlag
                rstruct.threshold = 1;
            else
                rstruct.threshold = th_res.threshold;
            end

            if updatePerf
                rstruct = AnalysisFiles.calculatePerformanceMetrics(...
                    rstruct.callset, rstruct.threshold, rstruct, -1, rstruct.th_scan_min, rstruct.th_scan_max);
            else
                thidx = 0;
                if th_res.threshold > 0
                    thidx = RNAUtils.findThresholdIndex(rstruct.threshold, spot_table(:,1)');
                end

                if isfield(rstruct, 'performance_trimmed')
                    if thidx > 0
                        rstruct.fscore_autoth_trimmed = rstruct.performance_trimmed{thidx, 'fScore'};
                    else
                        rstruct.fscore_autoth_trimmed = NaN;
                    end
                end
                if isfield(rstruct, 'performance')
                    if thidx > 0
                        rstruct.fscore_autoth = rstruct.performance{thidx, 'fScore'};
                    else
                        rstruct.fscore_autoth = NaN;
                    end
                end
            end

            analysis.results_hb = rstruct; 
        end

    end
end