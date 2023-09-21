%
%%

classdef RNAThreshold
    
    methods (Static)

        %------------ Internal (Moved from RNA_Threshold_Common) ------------

        %%
        %
        function param_struct = genEmptyThresholdParamStruct()
            param_struct = struct("window_pos", 0.5);
            param_struct.window_sizes = [3:3:21];
            param_struct.suggestion = [0 0]; %Suggested thresh range from outside data
            param_struct.mad_factor_min = -1.0;
            param_struct.mad_factor_max = 1.0;
            param_struct.spline_iterations = 3;
            param_struct.verbosity = 0;
            param_struct.sample_spot_table = [];
            param_struct.control_spot_table = [];
            param_struct.test_data = false;
            param_struct.test_diff = false;
            param_struct.test_winsc = true;
            param_struct.reweight_fit = true;
            param_struct.fit_to_log = true;
            param_struct.fit_strat = 'default';
            param_struct.fit_ri_weight = 0.0;
            param_struct.madth_weight = 0.0;
            param_struct.fit_weight = 1.0;
            param_struct.std_factor = 0.0;
            param_struct.min_log_diff = 0.15;
        end
        
        %%
        %
        function sub_struct = genEmptyThresholdInfoStruct()
            sub_struct = struct("spline_fit", []);
            sub_struct.spline_knot_x = 0;
            sub_struct.spline_knot_y = 0;
            sub_struct.max_index = 0; %Index in data of maximum y
            sub_struct.median = 0;
            sub_struct.mad = 0;
            sub_struct.med_suggested_threshold = [];
            sub_struct.scan_value_threshold = [];
            sub_struct.medth_min = 0;
            sub_struct.medth_max = 0;
            sub_struct.medth_avg = 0;
            sub_struct.medth_std = 0;
            sub_struct.spline_okay = false;
            sub_struct.med_okay = false;
            sub_struct.log_used_spline = false;
        end
        
        %%
        %
        function tres_struct = genEmptyThresholdResultStruct()
            res = struct("window_pos", 0.0);
            res.window_sizes = [];
            res.suggestion = [0 0];
            res.mad_factor_min = -1.0;
            res.mad_factor_max = 1.5;
            res.spline_iterations = 0;
            res.verbosity = 0;
            res.control_floor = 0;
            res.x = []; %x values (thresholds) n x 1 mtx
            res.window_scores = []; %Win score plots. Rows are values, columns are window sizes
            res.window_scores_ctrl = [];
            res.threshold = 0; %Overall suggestion
            res.test_data = [];
            res.test_diff = [];
            res.struct_ver = 3;
            res.reweight_fit = true;
            res.fit_to_log = true;
            res.fit_strat = 'default';
            res.fit_ri_weight = 0.0;
            res.madth_weight = 0.0;
            res.fit_weight = 1.0;
            res.std_factor = 0.0;
            res.min_log_diff = 0.5;
            tres_struct = res;
        end
        
        %%
        function thresh_info = thresholdTestCurveDoFit(data, data_trimmed, thresh_info, params, min_points)
            thresh_info.spline_okay = true;
            
            thresh_info.spline_fit = [];
            thresh_info.spline_knot_x = 0;
            thresh_info.spline_knot_y = 0;
            
            fitdata = data;
            point_count = size(data,1);
            
            if thresh_info.log_used_spline
                fitdata(:,2) = log10(fitdata(:,2));
                
                %Remove invalid values...
                [okay_rows,~] = find(isfinite(fitdata(:,2)));
                if ~isempty(okay_rows)
                    keep_row_count = size(okay_rows,1);
                    if keep_row_count < min_points
                        thresh_info.spline_okay = false;
                        return;
                    end
                else
                    thresh_info.spline_okay = false;
                    return;
                end
                
                fitdata = fitdata(okay_rows,:);
                [okay_rows,~] = find(fitdata(:,2) > -5); %Remove extremes
                if ~isempty(okay_rows)
                    keep_row_count = size(okay_rows,1);
                    if keep_row_count < min_points
                        thresh_info.spline_okay = false;
                        return;
                    end
                else
                    thresh_info.spline_okay = false;
                    return;
                end
                
                try
                    fitdata = fitdata(okay_rows,:);
                    
                    local_count = size(fitdata,1);
                    [~,maxidx] = max(fitdata(:,2),[],'omitnan');
                    thresh_info.max_index = maxidx;
                    fitdata = fitdata(maxidx:local_count,:);
                    
                    %Re-interpolate the removed values
                    min_x = fitdata(1,1);
                    max_x = fitdata(size(fitdata,1),1);
                    xx = [min_x:1:max_x].';
                    yy = pchip(fitdata(:,1), fitdata(:,2), xx);
                    fitdata = NaN(size(xx,1),2);
                    fitdata(:,1) = xx;
                    fitdata(:,2) = yy;
                catch ME
                    thresh_info.spline_okay = false;
                    return;
                end
            else
                fitdata = data_trimmed;
                [okay_rows,~] = find(isfinite(fitdata(:,2)));
                if ~isempty(okay_rows)
                    keep_row_count = size(okay_rows,1);
                    if keep_row_count < min_points
                        thresh_info.spline_okay = false;
                        return;
                    end
                else
                    thresh_info.spline_okay = false;
                    return;
                end

                fitdata = fitdata(okay_rows,:);
            end
            
            if params.verbosity > 0
                fprintf("Fitting two-piece linear spline...\n");
            end
            %thresh_info.spline_fit = Seglr2.fitTo(data_trimmed, thresh_info.medth_min, thresh_info.medth_max, params.spline_iterations, verbosity);
            if params.reweight_fit
                bpstrat = 'weighted_avg';
            else
                bpstrat = 'segrsq';
            end
            
            %Trycatch, try again as non log if tried log
            try
                thresh_info.spline_fit = Seglr2.fitToSpotCurve(fitdata, thresh_info.medth_min, thresh_info.medth_max, params.verbosity, params.fit_strat, bpstrat, params.spline_iterations);
                if ~isempty(thresh_info.spline_fit)
                    thresh_info.spline_knot_x = fitdata(thresh_info.spline_fit.break_index,1);
                    
                    %Adjust break index for original data...
                    for i = 1:point_count
                        if data(i,1) >= thresh_info.spline_knot_x
                            thresh_info.spline_fit.break_index = i;
                            break;
                        end
                    end
                    
                    %Calculate the y coord of knot.
                    thresh_info.spline_knot_y = (thresh_info.spline_fit.right.slope * thresh_info.spline_knot_x) + thresh_info.spline_fit.right.yintr;
                end
            catch ME
                thresh_info.spline_okay = false;
                return;
            end
        end
        
        %%
        function thresh_info = thresholdTestCurve(data, thresh_info, params)
            
            thresh_info.median = median(data(:,2), 'omitnan');
            thresh_info.mad = mad(data(:,2),1);
            thresh_info.scan_value_threshold = thresh_info.median + (thresh_info.mad .* [params.mad_factor_min:0.25:params.mad_factor_max]);
            thresh_info.log_used_spline = params.fit_to_log;
            
            %Trim everything before the initial peak
            point_count = size(data,1);
            [~,maxidx] = max(data(:,2),[],'omitnan');
            thresh_info.max_index = maxidx; 
            if (point_count - maxidx) < 5
                %Not enough to work with. Return.
                thresh_info.med_okay = false;
                thresh_info.spline_okay = false;
                return;
            end
            data_trimmed = data(maxidx:point_count,:);
            
            %Run median-based threshold determination
            %This could probably be fully vectorized?
            data_smoothed(:,1) = smooth(data_trimmed(:,2));
            mfcount = size(thresh_info.scan_value_threshold,2);
            thresh_info.med_suggested_threshold = NaN(1,mfcount);
            bool_has_nonnan = false;
            for i = 1:mfcount
                findval = find(data_smoothed(:,1) <= thresh_info.scan_value_threshold(1,i), 1);
                if ~isempty(findval)
                    thresh_info.med_suggested_threshold(1,i) = data(findval+maxidx-1,1);
                    bool_has_nonnan = true;
                end
            end

            %Evaluate medth values.
            if bool_has_nonnan
                thresh_info.medth_avg = mean(thresh_info.med_suggested_threshold,'all','omitnan');
                thresh_info.medth_std = std(thresh_info.med_suggested_threshold,0,'all','omitnan');
                medstd_factor = 2.0;
                medstd_amt = thresh_info.medth_std * medstd_factor;
                thresh_info.medth_min = round(thresh_info.medth_avg - medstd_amt);
                if thresh_info.medth_min < 1
                    thresh_info.medth_min = 1;
                end
                maxt = data(point_count,1);
                thresh_info.medth_max = round(thresh_info.medth_avg + medstd_amt);
                if thresh_info.medth_max > maxt
                    thresh_info.medth_max = maxt;
                end
                thresh_info.med_okay = true;
            else
                thresh_info.medth_avg = NaN;
                thresh_info.medth_std = NaN;
                thresh_info.medth_min = data(1,1);
                thresh_info.medth_max = data(point_count,1);
            end
            
            %If user specified (spline itr > 0), fit spline
            if params.spline_iterations > 0
                thresh_info = RNA_Threshold_Common.thresholdTestCurveDoFit(data, data_trimmed, thresh_info, params, 5);
                
                if params.fit_to_log
                    if ~thresh_info.spline_okay
                        thresh_info.log_used_spline = false;
                        thresh_info = RNA_Threshold_Common.thresholdTestCurveDoFit(data, data_trimmed, thresh_info, params, 5);
                    end
                end
                
                if ~thresh_info.spline_okay
                    thresh_info.spline_fit = [];
                    thresh_info.spline_knot_x = 0;
                    thresh_info.spline_knot_y = 0;
                end
            end
            
        end

        %%
        % (Description)
        %
        % ARGS
        %
        %
        % RETURN
        %
        function [window_scores, adj_window_size, adj_window_pos] = calculateWindowScores(diff_curve, window_size, window_pos)
            %Take some counts
            P_count = size(diff_curve,1);
            
            %Calculate window position shift.
            winshift = 0;
            if window_size < 1
                window_size = 1;
            elseif window_size > (P_count - 2)
                window_size = (P_count - 2);
            end
            adj_window_size = window_size;
            
            if window_pos > 0.0
                if window_pos >= 1.0
                    winshift = window_size;
                    adj_window_pos = 1.0;
                else
                    winshift = round(window_size * window_pos);
                    adj_window_pos = window_pos;
                end
            else
                adj_window_pos = 0.0;
            end
            
            %Calculate window scores for sample
            winout = NaN(P_count, 1);
            winmax = P_count;
            for i = 1:winmax
                w_back = i;
                w_front = i + window_size - 1;
                if w_back < 1
                    w_back = 1;
                end
                if w_front > P_count
                    w_front = P_count;
                end
                wmean = mean(diff_curve(w_back:w_front,1));
                winout(i) = var(diff_curve(w_back:w_front,1)) / wmean;
            end
            
            %Shift back
            win_adj = NaN(P_count, 1);
            for i = 1:(P_count - winshift)
                win_adj(i+winshift) = winout(i);
            end
            window_scores = win_adj;
        end
        
        %%
        % (Description)
        %
        % ARGS
        %
        %
        % RETURN
        %
        function threshold_results = estimateThreshold(parameter_info)
            %Generate return struct
            threshold_results = RNA_Threshold_Common.genEmptyThresholdResultStruct();
            %threshold_results.window_size = parameter_info.window_size;
            threshold_results.window_sizes = parameter_info.window_sizes;
            threshold_results.window_pos = parameter_info.window_pos;
            threshold_results.suggestion = parameter_info.suggestion;
            %threshold_results.mad_factor = parameter_info.mad_factor;
            threshold_results.mad_factor_min = parameter_info.mad_factor_min;
            threshold_results.mad_factor_max = parameter_info.mad_factor_max;
            threshold_results.spline_iterations = parameter_info.spline_iterations;
            threshold_results.verbosity = parameter_info.verbosity;
            threshold_results.reweight_fit = parameter_info.reweight_fit;
            threshold_results.fit_to_log = parameter_info.fit_to_log;
            threshold_results.fit_strat = parameter_info.fit_strat;
            threshold_results.fit_ri_weight = parameter_info.fit_ri_weight;
            threshold_results.madth_weight = parameter_info.madth_weight;
            threshold_results.fit_weight = parameter_info.fit_weight;
            threshold_results.std_factor = parameter_info.std_factor;
            threshold_results.min_log_diff = parameter_info.min_log_diff;
            
            mweight = threshold_results.madth_weight;
            fweight = threshold_results.fit_weight;
            iweight = threshold_results.fit_ri_weight;
            if mweight < 0.0; mweight = 0.0; end
            if fweight < 0.0; fweight = 0.0; end
            if iweight < 0.0; iweight = 0.0; end
            if strcmp(threshold_results.fit_strat, 'three_piece')
                iweight = 0.0;
            end
            if threshold_results.spline_iterations < 1
                fweight = 0.0; iweight = 0.0; mweight = 1.0;
            end
            wsum = mweight + fweight + iweight;
            if wsum ~= 1.0
                if mweight == 0.0 & fweight == 0.0 & iweight == 0.0
                    fweight = 1.0;
                else
                    fweight = fweight/wsum;
                    mweight = mweight/wsum;
                    iweight = iweight/wsum;
                end
            end
            threshold_results.fit_weight = fweight;
            threshold_results.madth_weight = mweight;
            threshold_results.fit_ri_weight = iweight;
            
            %Start by calculating the diff
            spotcount_table = double(parameter_info.sample_spot_table);
            threshold_results.x = spotcount_table(:,1);
            deriv1 = diff(spotcount_table(:,2));
            deriv1 = smooth(deriv1);
            deriv1 = abs(deriv1);
            T_sample = size(spotcount_table,1);
            P_sample = T_sample - 1;
            
            %Repeat for control, if present
            ctrl_spotcount_table = double(parameter_info.control_spot_table);
            if ~isempty(ctrl_spotcount_table)
                ctrlderiv = diff(ctrl_spotcount_table(:,2));
                ctrlderiv = smooth(ctrlderiv);
                ctrlderiv = abs(ctrlderiv);
                T_control = size(ctrl_spotcount_table,1);
                P_control = T_control - 1;
            else
                ctrlderiv = [];
                T_control = 0;
                P_control = 0;
            end
       
            %Find control floor.
            if ~isempty(ctrlderiv)
                if parameter_info.verbosity > 0
                    fprintf("Finding noise floor...\n");
                end
                [~, cdmaxidx] = max(ctrlderiv(:,1), [], 'omitnan');
                cdtest = ctrlderiv(cdmaxidx:P_control,1);
                
                [findres, ~] = find(cdtest < 1,1);
                if isempty(findres)
                    [findres, ~] = find(cdtest < 2,1);
                end
                if ~isempty(findres)
                    threshold_results.control_floor = spotcount_table(findres+cdmaxidx,1);
                end
            end
            
            %Alloc window score mtxs in return struct
            wincount = 0;
            if ~isempty(parameter_info.window_sizes)
                wincount = size(parameter_info.window_sizes,2);
                threshold_results.window_scores = NaN(P_sample,wincount);
                if ~isempty(ctrlderiv)
                    threshold_results.window_scores_ctrl = NaN(P_control,wincount);
                end
                threshold_results.test_winsc(1,wincount) = RNA_Threshold_Common.genEmptyThresholdInfoStruct();
            end
           
            %For each window size...
            for i = 1:wincount
                wsz = parameter_info.window_sizes(1,i);
                [threshold_results.window_scores(:,i), threshold_results.window_sizes(1,i), threshold_results.window_pos] =...
                    RNA_Threshold_Common.calculateWindowScores(deriv1, wsz, threshold_results.window_pos);
                
                %Repeat for control
                if ~isempty(ctrlderiv)
                    [threshold_results.window_scores_ctrl(:,i), ~, ~] =...
                        RNA_Threshold_Common.calculateWindowScores(ctrlderiv, wsz, threshold_results.window_pos);
                end
            end
           
            %Test the spot curve, the diff, and the win score curve (as
            %   wanted)
            %If window scores are not useful, we need to turn on curve and
            %diff use
            med_okay = 0;
            fit_okay = 0;
            if parameter_info.test_winsc
                if parameter_info.verbosity > 0
                    fprintf("Testing winscore curve...\n");
                end
                for i = 1:wincount
                    thresh_info = RNA_Threshold_Common.genEmptyThresholdInfoStruct();
                    threshold_results.test_winsc(1,i) = thresh_info;
                    %threshold_results.test_winsc(1,i) = thresh_info;
                    winscores_test = NaN(P_sample,2);
                    winscores_test(:,1) = threshold_results.x(1:P_sample,1);
                    winscores_test(:,2) = threshold_results.window_scores(:,i);
                    
                    ws_results = ...
                        RNA_Threshold_Common.thresholdTestCurve(winscores_test, thresh_info, parameter_info);
                    threshold_results.test_winsc(1,i) = ws_results;
                    
                    if ~isempty(ws_results)
                        if ws_results.med_okay; med_okay = med_okay + 1; end
                        if ws_results.spline_okay; fit_okay = fit_okay + 1; end
                    end
                end
            else 
                threshold_results.test_winsc = [];
            end
            
            %Check med_okay and fit_okay to determine whether to auto
            %turn on data and diff checks
            if threshold_results.madth_weight > 0.0
                if med_okay < 1
                    parameter_info.test_data = true;
                    parameter_info.test_diff = true;
                end
            end

            if threshold_results.fit_weight > 0.0 | threshold_results.fit_ri_weight > 0.0
                if fit_okay < 1
                    parameter_info.test_data = true;
                    parameter_info.test_diff = true;
                end
            end
            
            if parameter_info.test_data
                if parameter_info.verbosity > 0
                    fprintf("Testing sample spot curve...\n");
                end
                threshold_results.test_data = ...
                    RNA_Threshold_Common.thresholdTestCurve(spotcount_table, RNA_Threshold_Common.genEmptyThresholdInfoStruct(), parameter_info);
            end
            
            if parameter_info.test_diff
                if parameter_info.verbosity > 0
                    fprintf("Testing diff curve...\n");
                end
                diff_curve = NaN(P_sample,1);
                diff_curve(:,1) = spotcount_table(1:P_sample,1);
                diff_curve(:,2) = deriv1(:,1);
                threshold_results.test_diff = ...
                    RNA_Threshold_Common.thresholdTestCurve(diff_curve, RNA_Threshold_Common.genEmptyThresholdInfoStruct(), parameter_info);
            end
            
            
            %Determine a best threshold from test info.
            mscores = RNAThreshold.getAllMedThresholds(threshold_results);
            fscores = RNAThreshold.getAllFitThresholds(threshold_results);
            iscores = RNAThreshold.getAllRightISectThresholds(threshold_results);
            allscores = [];

            if ~isempty(mscores)
                mfac = mean(mscores, 'all', 'omitnan') + (std(mscores, 0, 'all', 'omitnan') * threshold_results.std_factor);
                allscores = cat(2, allscores, mscores);
            else
                mfac = 0.0;
            end

            if ~isempty(fscores)
                ffac = mean(fscores, 'all', 'omitnan') + (std(fscores, 0, 'all', 'omitnan') * threshold_results.std_factor);
                allscores = cat(2, allscores, fscores);
            else
                ffac = 0.0;
            end

            if ~isempty(iscores)
                ifac = mean(iscores, 'all', 'omitnan') + (std(iscores, 0, 'all', 'omitnan') * threshold_results.std_factor);
                allscores = cat(2, allscores, iscores);
            else
                ifac = 0.0;
            end
            threshold_results.threshold = round((mfac * mweight) + (ffac * fweight) + (ifac * iweight));

            sugg_min = threshold_results.suggestion(1);
            sugg_max = threshold_results.suggestion(2);
            if sugg_min > 0
                if sugg_max > 0
                    %Find std of all suggested thresholds to generate a
                    %range.
                    %Find the overlap (if there is one) between suggested
                    %range and this image's range.
                    allstd = std(allscores, 0, 'all', 'omitnan');
                    det_min = (threshold_results.threshold - allstd);
                    det_max = (threshold_results.threshold + allstd);
                    ovl_min = max(det_min, sugg_min);
                    ovl_max = min(det_max, sugg_max);
                    if ovl_min > ovl_max
                        %No overlap. Just pull up or down thresh.
                        if sugg_min >= det_max
                            %Suggestion is higher
                            threshold_results.threshold = round((sugg_min + det_max) / 2);
                        else
                            %Suggestion is lower.
                            threshold_results.threshold = round((sugg_max + det_min) / 2);
                        end
                    else
                        %Overlap. Set thresh to middle of overlap.
                        threshold_results.threshold = round((ovl_min + ovl_max) / 2);
                    end
                end
            end
            
            if threshold_results.threshold < threshold_results.control_floor
                if parameter_info.verbosity > 0
                    fprintf("Warning: Auto selected threshold below noise floor. Adjusting.\n");
                end
                threshold_results.threshold = threshold_results.control_floor;
            end
        end

        %------------ Thresholder Run Wrappers ------------
        
        function threshold_results = runDefaultParameters(spot_count_table, ctrl_count_table)
            param_struct = RNA_Threshold_Common.genEmptyThresholdParamStruct();
            param_struct.sample_spot_table = spot_count_table;
            
            if nargin > 1
                param_struct.control_spot_table = ctrl_count_table;
            end
            
            threshold_results = RNA_Threshold_Common.estimateThreshold(param_struct);
        end
        
        function threshold_results = runDefaultParametersSpotsRun(rnaspots_run, verbosity)
            if nargin < 2
                verbosity = 0;
            end
            
            param_struct = RNA_Threshold_Common.genEmptyThresholdParamStruct();
            [rnaspots_run, param_struct.sample_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
            [~, param_struct.control_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Control();
            param_struct.verbosity = verbosity;
            threshold_results = RNA_Threshold_Common.estimateThreshold(param_struct);
        end
        
        function threshold_results = runSavedParameters(rnaspots_run, verbosity, spot_table, ctrl_table)
            if nargin < 2
                verbosity = 0;
            end
            
            param_struct = RNA_Threshold_Common.genEmptyThresholdParamStruct();
            
            if nargin >= 3 & ~isempty(spot_table)
                param_struct.sample_spot_table = double(spot_table);
            else
                [rnaspots_run, spot_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
                param_struct.sample_spot_table = double(spot_table);
            end
            
            
            if nargin >= 4 & ~isempty(ctrl_table)
                param_struct.control_spot_table = ctrl_table;
            else
                [rnaspots_run, param_struct.control_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Control();
            end
            
            param_struct.verbosity = verbosity;
            
            winmin = rnaspots_run.ttune_winsz_min;
            winmax = rnaspots_run.ttune_winsz_max;
            if winmin == 0 & winmax == 0
                param_struct.window_sizes = [];
            else
                if winmin < 2; winmin = 2; end
                winincr = rnaspots_run.ttune_winsz_incr;
                if winincr < 1; winincr = 1; end
                if winmax < winmin; winmax = winmin+winincr; end
                param_struct.window_sizes = [winmin:winincr:winmax];
            end

            madmin = rnaspots_run.ttune_madf_min;
            if isnan(madmin); madmin = -1.0; end
            madmax = rnaspots_run.ttune_madf_max;
            if isnan(madmax); madmax = madmin+0.5; end
            param_struct.mad_factor_min = madmin;
            param_struct.mad_factor_max = madmax;
            
            param_struct.test_data = rnaspots_run.ttune_use_rawcurve;
            param_struct.test_diff = rnaspots_run.ttune_use_diffcurve;
            
            spline_itr = rnaspots_run.ttune_spline_itr;
            if spline_itr < 0; spline_itr = 0; end
            param_struct.spline_iterations = spline_itr;
            
            param_struct.reweight_fit = rnaspots_run.ttune_reweight_fit;
            param_struct.fit_to_log = rnaspots_run.ttune_fit_to_log;
            param_struct.fit_ri_weight = rnaspots_run.ttune_thweight_fisect;
            param_struct.madth_weight = rnaspots_run.ttune_thweight_med;
            param_struct.fit_weight = rnaspots_run.ttune_thweight_fit;
            param_struct.std_factor = rnaspots_run.ttune_std_factor;
            
            if isempty(rnaspots_run.ttune_fit_strat) | (rnaspots_run.ttune_fit_strat == 0)
                param_struct.fit_strat = 'default';
            elseif rnaspots_run.ttune_fit_strat == 1
                param_struct.fit_strat = 'slow';
            elseif rnaspots_run.ttune_fit_strat == 2
                param_struct.fit_strat = 'section_fit';
            elseif rnaspots_run.ttune_fit_strat == 3
                param_struct.fit_strat = 'three_piece';
            else
                param_struct.fit_strat = 'default';
            end
            
            threshold_results = RNA_Threshold_Common.estimateThreshold(param_struct);
        end
        
        function threshold_results = runWithPreset(spot_count_table, ctrl_count_table, preset_index)
            presets = RNAThreshold.loadPresets();
            if isempty(presets); return; end
            if preset_index < 1; return; end
            if preset_index > size(presets,2); return; end
            preset_struct = presets(preset_index);
            
            param_struct = RNA_Threshold_Common.genEmptyThresholdParamStruct();
            param_struct.sample_spot_table = spot_count_table;
            param_struct.control_spot_table = ctrl_count_table;
            
            win_min = preset_struct.ttune_winsz_min;
            win_max = preset_struct.ttune_winsz_max;
            win_incr = preset_struct.ttune_winsz_incr;
            
            param_struct.window_sizes = [win_min:win_incr:win_max];
            param_struct.reweight_fit = preset_struct.ttune_reweight_fit;
            param_struct.fit_to_log = preset_struct.ttune_fit_to_log;
            param_struct.madth_weight = preset_struct.ttune_thweight_med;
            param_struct.fit_weight = preset_struct.ttune_thweight_fit;
            param_struct.fit_ri_weight = preset_struct.ttune_thweight_fisect;
            param_struct.std_factor = preset_struct.ttune_std_factor;
            param_struct.spline_iterations = preset_struct.ttune_spline_itr;
            param_struct.test_data = preset_struct.ttune_use_rawcurve;
            param_struct.test_diff = preset_struct.ttune_use_diffcurve;
            param_struct.mad_factor_min = preset_struct.ttune_madf_min;
            param_struct.mad_factor_max = preset_struct.ttune_madf_max;
            
            if isempty(preset_struct.ttune_fit_strat) | (preset_struct.ttune_fit_strat == 0)
                param_struct.fit_strat = 'default';
            elseif preset_struct.ttune_fit_strat == 1
                param_struct.fit_strat = 'slow';
            elseif preset_struct.ttune_fit_strat == 2
                param_struct.fit_strat = 'section_fit';
            elseif preset_struct.ttune_fit_strat == 3
                param_struct.fit_strat = 'three_piece';
            else
                param_struct.fit_strat = 'default';
            end
            
            threshold_results = RNA_Threshold_Common.estimateThreshold(param_struct);
        end
        
        %------------ Parameters ------------

        function param_info = paramsFromSpotsrun(rnaspots_run)
            param_info = RNA_Threshold_Common.genEmptyThresholdParamStruct();
            
            winmin = rnaspots_run.ttune_winsz_min;
            winmax = rnaspots_run.ttune_winsz_max;
            if winmin == 0 & winmax == 0
                param_info.window_sizes = [];
            else
                if winmin < 2; winmin = 2; end
                winincr = rnaspots_run.ttune_winsz_incr;
                if winincr < 1; winincr = 1; end
                if winmax < winmin; winmax = winmin+winincr; end
                param_info.window_sizes = [winmin:winincr:winmax];
            end
            
            madmin = rnaspots_run.ttune_madf_min;
            if isnan(madmin); madmin = -1.0; end
            madmax = rnaspots_run.ttune_madf_max;
            if isnan(madmax); madmax = madmin+0.5; end
            param_info.mad_factor_min = madmin;
            param_info.mad_factor_max = madmax;
            
            param_info.test_data = rnaspots_run.ttune_use_rawcurve;
            param_info.test_diff = rnaspots_run.ttune_use_diffcurve;
            
            spline_itr = rnaspots_run.ttune_spline_itr;
            if spline_itr < 0; spline_itr = 0; end
            param_info.spline_iterations = spline_itr;
            
            if isempty(rnaspots_run.ttune_fit_strat) | (rnaspots_run.ttune_fit_strat == 0)
                param_info.fit_strat = 'default';
            elseif rnaspots_run.ttune_fit_strat == 1
                param_info.fit_strat = 'slow';
            elseif rnaspots_run.ttune_fit_strat == 2
                param_info.fit_strat = 'section_fit';
            elseif rnaspots_run.ttune_fit_strat == 3
                param_info.fit_strat = 'three_piece';
            else
                param_info.fit_strat = 'default';
            end
            
            param_info.reweight_fit = rnaspots_run.ttune_reweight_fit;
            param_info.fit_to_log = rnaspots_run.ttune_fit_to_log;
            param_info.fit_ri_weight = rnaspots_run.ttune_thweight_fisect;
            param_info.madth_weight = rnaspots_run.ttune_thweight_med;
            param_info.fit_weight = rnaspots_run.ttune_thweight_fit;
            param_info.std_factor = rnaspots_run.ttune_std_factor;
        end
        
        %------------ Presets ------------

        function preset_struct = genPresetStruct()
            preset_struct = struct("ttune_fit_strat", 0);
        	preset_struct.ttune_winsz_min = 3;
            preset_struct.ttune_winsz_max = 21;
            preset_struct.ttune_winsz_incr = 3;
            preset_struct.ttune_reweight_fit = false;
            preset_struct.ttune_fit_to_log = true;
            preset_struct.ttune_thweight_med = 0.0;
            preset_struct.ttune_thweight_fit = 1.0;
            preset_struct.ttune_thweight_fisect = 0.0;
            preset_struct.ttune_std_factor = 0.0;
            preset_struct.ttune_spline_itr = 3;
            preset_struct.ttune_use_rawcurve = false;
            preset_struct.ttune_use_diffcurve = false;
            preset_struct.ttune_madf_min = -1.0;
            preset_struct.ttune_madf_max = 1.0;
        end
        
        function presets = loadPresets()
            filepath = ['.' filesep 'core' filesep 'ths_presets.mat'];
            if isfile(filepath)
                load(filepath, 'presets');
            else
                presets = [];
            end
        end
        
        function presets = setPreset(preset_struct, preset_index)
            presets = RNAThreshold.loadPresets();
            presets(preset_index) = preset_struct;
            filepath = ['.' filesep 'core' filesep 'ths_presets.mat'];
            save(filepath, 'presets');
        end
        
        function spotsrun = applyPreset(spotsrun, preset_index)
            presets = RNAThreshold.loadPresets();
            if isempty(presets); return; end
            if preset_index < 1; return; end
            if preset_index > size(presets,2); return; end
            preset_struct = presets(preset_index);
            
            spotsrun.ttune_fit_strat = preset_struct.ttune_fit_strat;
            spotsrun.ttune_winsz_min = preset_struct.ttune_winsz_min;
            spotsrun.ttune_winsz_max = preset_struct.ttune_winsz_max;
            spotsrun.ttune_winsz_incr = preset_struct.ttune_winsz_incr;
            spotsrun.ttune_reweight_fit = preset_struct.ttune_reweight_fit;
            spotsrun.ttune_fit_to_log = preset_struct.ttune_fit_to_log;
            spotsrun.ttune_thweight_med = preset_struct.ttune_thweight_med;
            spotsrun.ttune_thweight_fit = preset_struct.ttune_thweight_fit;
            spotsrun.ttune_thweight_fisect = preset_struct.ttune_thweight_fisect;
            spotsrun.ttune_std_factor = preset_struct.ttune_std_factor;
            spotsrun.ttune_spline_itr = preset_struct.ttune_spline_itr;
            spotsrun.ttune_use_rawcurve = preset_struct.ttune_use_rawcurve;
            spotsrun.ttune_use_diffcurve = preset_struct.ttune_use_diffcurve;
            spotsrun.ttune_madf_min = preset_struct.ttune_madf_min;
            spotsrun.ttune_madf_max = preset_struct.ttune_madf_max;
        end
        
        %------------ Special Case ------------

        function [peak_x, min_x, max_x] = threshSuggestionFromFScores(fscores)
            if isempty(fscores)
                peak_x = 0; min_x = 0; max_x = 0;
                return;
            end

            [maxval, idx] = max(fscores(:,2), [], 'omitnan');
            minval = min(fscores(:,2), [], 'omitnan');
            roff = (maxval - minval) * 0.05;
            peak_x = f_score(idx, 1);
            min_x = peak_x; max_x = peak_x;
            
            scount = size(fscores,1);
            t = maxval - roff;
            for i = idx:-1:1
                if fscores(i,2) <= t
                    min_x = fscore(i,1);
                    break;
                end
            end
            for i = idx:scount
                if fscores(i,2) <= t
                    max_x = fscore(i,1);
                    break;
                end
            end
        end

        %------------ Getters ------------

        function score_list = getAllMedThresholds(threshold_results)
            curve_info_list = RNAThreshold.getAllCurveResultStructs(threshold_results);
            if isempty(curve_info_list); return; end
            
            curve_count = size(curve_info_list,2);
            madf_count = 0;
            
            for i = 1:curve_count
                if isempty(curve_info_list(i)); continue; end
                if ~curve_info_list(i).med_okay; continue; end
                if ~isempty(curve_info_list(i).med_suggested_threshold)
                    madf_count = size(curve_info_list(i).med_suggested_threshold,2);
                    break;
                end
            end
            
            alloc = madf_count * curve_count;
            score_list = NaN(1, alloc);
            
            i = 1;
            for c = 1:curve_count
                if isempty(curve_info_list(c)); continue; end
                if ~curve_info_list(c).med_okay; continue; end
                if ~isempty(curve_info_list(c).med_suggested_threshold)
                    score_list(1,i:i+madf_count-1) = curve_info_list(c).med_suggested_threshold(1,:);
                    i = i + madf_count;
                end
            end
        end
        
        function score_list = getAllFitThresholds(threshold_results)
            curve_info_list = RNAThreshold.getAllCurveResultStructs(threshold_results);
            if isempty(curve_info_list); return; end
            
            curve_count = size(curve_info_list,2);
            score_list = NaN(1, curve_count);
            
            i = 1;
            for c = 1:curve_count
                if isempty(curve_info_list(c)); continue; end
                if ~curve_info_list(c).spline_okay; continue; end
                if ~isempty(curve_info_list(c).spline_fit)
                    score_list(i) = curve_info_list(c).spline_knot_x;
                    i = i + 1;
                end
            end
        end

        function score_list = getAllRightISectThresholds(threshold_results)
            curve_info_list = RNAThreshold.getAllCurveResultStructs(threshold_results);
            if isempty(curve_info_list); return; end
            
            curve_count = size(curve_info_list,2);
            score_list = NaN(1, curve_count);
            
            i = 1;
            for c = 1:curve_count
                if isempty(curve_info_list(c)); continue; end
                if ~curve_info_list(c).spline_okay; continue; end
                if ~isempty(curve_info_list(c).spline_fit)
                    score_list(i) = curve_info_list(c).spline_fit.rcurve_intr_x;
                    i = i + 1;
                end
            end
        end
        
        function score_list = getAllThresholdSuggestions(threshold_results)
            curve_info_list = RNAThreshold.getAllCurveResultStructs(threshold_results);
            if isempty(curve_info_list); return; end
            
            curve_count = size(curve_info_list,2);
            madf_count = 0;
            
            for c = 1:curve_count
                if isempty(curve_info_list(c)); continue; end
                if ~curve_info_list(c).med_okay; continue; end
                if ~isempty(curve_info_list(c).med_suggested_threshold)
                    madf_count = size(curve_info_list(c).med_suggested_threshold,2);
                    break;
                end
            end
            
            alloc = (madf_count * curve_count) + (curve_count * 2);
            score_list = NaN(1, alloc);
            
            i = 1;
            for c = 1:curve_count
                if isempty(curve_info_list(c)); continue; end
                if ~isempty(curve_info_list(c).med_suggested_threshold)
                    if ~curve_info_list(c).med_okay; continue; end
                    score_list(1,i:i+madf_count-1) = curve_info_list(c).med_suggested_threshold(1,:);
                    i = i + madf_count;
                end
                if ~isempty(curve_info_list(c).spline_fit)
                    if ~curve_info_list(c).spline_okay; continue; end
                    score_list(i) = curve_info_list(c).spline_knot_x;
                    i = i + 1;
                    
                    score_list(i) = curve_info_list(c).spline_fit.rcurve_intr_x;
                    i = i + 1;
                end
            end
        end

        function curve_info_list = getAllCurveResultStructs(threshold_results)
            curve_count = 0;
            wincount = 0;
            
            if ~isempty(threshold_results.test_data)
                curve_count = curve_count + 1; 
            end
            if ~isempty(threshold_results.test_diff)
                curve_count = curve_count + 1; 
            end
            
            if isfield(threshold_results, 'test_winsc') & ~isempty(threshold_results.test_winsc)
                wincount = size(threshold_results.test_winsc, 2); 
                curve_count = curve_count + wincount;
            end
            
            curve_info_list(curve_count) = RNA_Threshold_Common.genEmptyThresholdInfoStruct();
            i = 1;
            if ~isempty(threshold_results.test_data)
                curve_info_list(i) = threshold_results.test_data;
                i = i + 1;
            end
            if ~isempty(threshold_results.test_diff)
                curve_info_list(i) = threshold_results.test_diff;
                i = i + 1;
            end
            if wincount > 0
                for j = 1:wincount
                    curve_info_list(i) = threshold_results.test_winsc(1,j);
                    i = i + 1;
                end
            end
        end
        
        %------------ Plots ------------

        function [plot_handle, line_x] = draw_thres_plot(ax, data_x, data_y, curveres, include_madrange, color_base, color_dark, color_light, is_multi)
            %Plot curve
            point_count = size(data_y,1);
            plot_handle = plot(data_x(1:point_count,1), data_y,'LineWidth',2,'Color',color_base);
            hold on;
                
            %Plot fit
            thx = 0;
            if ~isempty(curveres.spline_fit)
                 line_x = data_x(1:curveres.spline_fit.break_index,1);
                 line_y = (line_x .* curveres.spline_fit.left.slope) + curveres.spline_fit.left.yintr;
                 plot(line_x, line_y,'LineStyle',':','LineWidth',1,'Color',color_light);
                    
                 line_x = data_x(curveres.spline_fit.break_index:point_count,1);
                 line_y = (line_x .* curveres.spline_fit.right.slope) + curveres.spline_fit.right.yintr;
                 plot(line_x, line_y,'LineStyle',':','LineWidth',1,'Color',color_light);
                 
                 thx = curveres.spline_knot_x;
            end
            line_x = [0 0 0];
            
            line_x(1,2) = curveres.medth_min;
            line_x(1,3) = curveres.medth_max;
            
            if include_madrange
                %Draw mad range lines
                linex = curveres.medth_min;
                if linex > 0 & ~is_multi
                    line([linex linex], get(ax,'YLim'),'Color',color_dark,'LineStyle',':','LineWidth',1);
                end
                
                linex = curveres.medth_max;
                if linex > 0 & ~is_multi
                    line([linex linex], get(ax,'YLim'),'Color',color_dark,'LineStyle',':','LineWidth',1);
                end
            end
            
            if isnan(thx) | thx <= 0
                %Grab threshold from mad scan avg?
                thx = round(curveres.medth_avg);
            end
                
            %Draw threshold line
            if thx > 0 & ~is_multi
                line([thx thx], get(ax,'YLim'),'Color',color_dark,'LineStyle','--','LineWidth',1);
            end
            line_x(1,1) = thx;
        end
        
        function fig_handle = resultPlotCombine(rnaspots_run, figno)
            %Plots all curves used for thresholding together on one plot.
            %Log10 scale
            %Also shows linear fits and chosen thresholds.
            %Does not show MADf ranges (too crowded)
            %Shows up to 7 curves.
            %If more than that, it will just return (too crowded).
            if isempty(rnaspots_run)
                fig_handle = [];
                return;
            end
            if isempty(rnaspots_run.threshold_results)
                fig_handle = [];
                return;
            end
            
            curve_count = 0;
            if rnaspots_run.ttune_use_rawcurve; curve_count = curve_count+1; end
            if rnaspots_run.ttune_use_diffcurve; curve_count = curve_count+1; end
            
            thres = rnaspots_run.threshold_results;
            if ~isempty(thres.window_sizes)
                curve_count = curve_count + size(thres.window_sizes, 2);
            end
            if curve_count < 1 | curve_count > 7
                fprintf("Too many or too few curves to plot! Returning...\n");
                fig_handle = [];
                return;
            end
            
            %Colors: Red, Orange/Yellow, Green, Aqua, Blue, Purple, Magenta
            colors_base = NaN(7,3);
            colors_light = NaN(7,3);
            colors_dark = NaN(7,3);
            
            colors_base(1,:) = [0.678 0.243 0.243]; %#ad3e3e Red
            colors_base(2,:) = [0.831 0.718 0.165]; %#d4b72a Yellow-Orange
            colors_base(3,:) = [0.243 0.678 0.243]; %#3ead3e Green
            colors_base(4,:) = [0.165 0.831 0.753]; %#2ad4c0 Aqua
            colors_base(5,:) = [0.243 0.243 0.678]; %#3e3ead Blue
            colors_base(6,:) = [0.490 0.129 0.749]; %#7d21bf Purple
            colors_base(7,:) = [0.678 0.243 0.678]; %#ad3ead Magenta
            
            for i = 1:7
                %This should output some interesting results...
                for j = 1:3
                    colors_light(i,j) = colors_base(i,j) + ((1.0 - colors_base(i,j))/2);
                    colors_dark(i,j) = colors_base(i,j) * 0.75;
                end
            end
            
            %Decide which color curve to use.
            usecolors = [1:1:7];
            remove_colors = 7 - curve_count;
            if remove_colors > 0
                for i = 1:remove_colors
                    ccount = size(usecolors,2);
                    mididx = round(ccount/2);
                    if mod(ccount,2) == 0
                        rmidx = mididx + 1;
                    else
                        rmidx = mididx - 1;
                    end
                    tempcolors = usecolors;
                    usecolors = zeros(1,ccount-1);
                    usecolors(:,1:rmidx-1) = tempcolors(:,1:rmidx-1);
                    usecolors(:,rmidx:ccount-1) = tempcolors(:,rmidx+1:ccount);
                end
            end
            clear tempcolors;
            
            %Initialize plot
            line_x = NaN(curve_count,3);
            c_idx = 1;
            legend_names = cell(1, curve_count);
            fig_handle = figure(figno);
            clf;
            ax = axes;
            
            %Raw curve (if applicable)
            spots_table = [];
            if rnaspots_run.ttune_use_rawcurve
                legend_names{1,c_idx} = 'Spot Count';
                
                %Plot curve
                [~, spots_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
                color_idx = usecolors(1,c_idx);
                log_y = log10(spots_table(:,2));
                log_y(log_y <= -5)= NaN;
                [plotlist(1,c_idx), lines] = RNAThreshold.draw_thres_plot(ax, spots_table(:,1), log_y, thres.test_data,...
                    false, colors_base(color_idx,:), colors_dark(color_idx,:), colors_light(color_idx,:), true);
                line_x(c_idx,:) = lines(1,:);
                c_idx = c_idx + 1;
            end
            
            if rnaspots_run.ttune_use_diffcurve
                legend_names{1,c_idx} = 'Smoothed Diff';
                
                %Plot curve
                if isempty(spots_table)
                    [~, spots_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
                end
                diff_curve = diff(spots_table(:,2));
                diff_curve = smooth(diff_curve);
                diff_curve = abs(diff_curve);
                
                color_idx = usecolors(1,c_idx);
                log_y = log10(diff_curve(:,1));
                log_y(log_y <= -5)= NaN;
                [plotlist(1,c_idx), lines] = RNAThreshold.draw_thres_plot(ax, spots_table(:,1), log_y, thres.test_diff,...
                    false, colors_base(color_idx,:), colors_dark(color_idx,:), colors_light(color_idx,:), true);
                line_x(c_idx,:) = lines(1,:);
                c_idx = c_idx + 1;
                clear diff_curve;
            end
            if ~isempty(spots_table)
                clear spots_table;
            end
            
            if ~isempty(thres.window_sizes)
                wincount = size(thres.window_sizes,2);
                for j = 1:wincount
                    legend_names{1,c_idx} = ['Window Size ' num2str(thres.window_sizes(1,j))];
                    win_curve = thres.window_scores(:,j);

                    color_idx = usecolors(1,c_idx);
                    log_y = log10(win_curve(:,1));
                    log_y(log_y <= -5)= NaN;
                    [plotlist(1,c_idx), lines] = RNAThreshold.draw_thres_plot(ax, thres.x, log_y, thres.test_winsc(1,j),...
                    false, colors_base(color_idx,:), colors_dark(color_idx,:), colors_light(color_idx,:), true);
                    line_x(c_idx,:) = lines(1,:);
                    c_idx = c_idx + 1;
                end
            end
            
            %Draw vertical lines.
            for i = 1:curve_count
                xval = line_x(i,1);
                color_idx = usecolors(1,i);
                if ~isnan(xval)
                    line([xval xval], get(ax,'YLim'),'Color',colors_dark(color_idx,:),'LineStyle','--','LineWidth',1);
                end
            end
            
            xlabel('Threshold');
            ylabel('log10(Score)');
            legend(plotlist, legend_names);
        end
        
        function fig_handles = resultPlotsIndiv(rnaspots_run, base_figno)
            %Plots curves used for thresholding each getting an individual
            %plot.
            %Log10 scale.
            %Shows linear fits, chosen thresholds, and scan ranges.
            if isempty(rnaspots_run)
                fig_handles = [];
                return;
            end
            if isempty(rnaspots_run.threshold_results)
                fig_handles = [];
                return;
            end
            
            curve_count = 0;
            if rnaspots_run.ttune_use_rawcurve; curve_count = curve_count+1; end
            if rnaspots_run.ttune_use_diffcurve; curve_count = curve_count+1; end
            
            thres = rnaspots_run.threshold_results;
            if ~isempty(thres.window_sizes)
                curve_count = curve_count + size(thres.window_sizes, 2);
            end
            
            if curve_count < 1
                %Nothing to draw.
                fig_handles = [];
                return;
            end
            
            c_idx = 1; figno = base_figno;
            color_base = [0.290, 0.290, 0.290];
            color_dark = [0.100, 0.100, 0.100];
            color_light = [0.500, 0.500, 0.500];
            
            spots_table = [];
            if rnaspots_run.ttune_use_rawcurve
                fig_handles(1,c_idx) = figure(figno);
                clf;
                ax = axes;  
                
                %Plot curve
                [~, spots_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
                log_y = log10(spots_table(:,2));
                log_y(log_y <= -5)= NaN;
                RNAThreshold.draw_thres_plot(ax, spots_table(:,1), log_y, thres.test_data,...
                    true, color_base, color_dark, color_light, false);
                
                xlabel('Threshold');
                ylabel('log10(Spot Count)');
                title('Spot Count');
                
                c_idx = c_idx + 1;
                figno = figno+1;
            end
            
            if rnaspots_run.ttune_use_diffcurve
                fig_handles(1,c_idx) = figure(figno);
                clf;
                ax = axes;
                
                %Plot curve
                if isempty(spots_table)
                    [~, spots_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
                end
                diff_curve = diff(spots_table(:,2));
                diff_curve = smooth(diff_curve);
                diff_curve = abs(diff_curve);
                
                log_y = log10(diff_curve(:,1));
                log_y(log_y <= -5)= NaN;
                RNAThreshold.draw_thres_plot(ax, spots_table(:,1), log_y, thres.test_diff,...
                    true, color_base, color_dark, color_light, false);
                
                xlabel('Threshold');
                ylabel('log10(|diff(Spot Count)|)');
                title('Smoothed Diff');
                
                c_idx = c_idx + 1;
                figno = figno+1;
                clear diff_curve;
            end
            if ~isempty(spots_table)
                clear spots_table;
            end
            
            if ~isempty(thres.window_sizes)
                wincount = size(thres.window_sizes,2);
                for j = 1:wincount
                    fig_handles(1,c_idx) = figure(figno);
                    clf;
                    ax = axes;
                    
                    win_curve = thres.window_scores(:,j);

                    log_y = log10(win_curve(:,1));
                    log_y(log_y <= -5)= NaN;
                    RNAThreshold.draw_thres_plot(ax, thres.x, log_y, thres.test_winsc(1,j),...
                        true, color_base, color_dark, color_light, false);

                    xlabel('Threshold');
                    ylabel('log10(Window Score)');
                    title(['Window Size: ' num2str(thres.window_sizes(1,j))]);
                    
                    c_idx = c_idx + 1;
                    figno = figno+1;
                end
            end
            
        end
        
        function rec = drawVerticalBox(ax, x_min, x_max, color, label, label_color)
            %Should go under plot since I'm not sure how to change alpha.
            %get(ax,'YLim')
            y_lim = get(ax,'YLim');
            x_lim = get(ax,'XLim');
            y_max = y_lim(1,2);
            y_min = y_lim(1,1);
            r_height = y_max - y_min;
            r_width = x_max - x_min;
            xlim_hi = x_lim(1,2);
            
            rec = rectangle('Position', [x_min y_min r_width r_height],...
                'FaceColor', color, 'LineStyle', 'none');
            %dim = [x_min/x_lim(1,2) 1.0 r_width/x_lim(1,2) 1.0];
            if nargin > 4
%                 rec = annotation('textbox', dim, 'FitBoxToText', 'off', ...
%                     'String', label, 'Color', label_color, 'LineStyle', 'none', 'HorizontalAlignment', 'center', ...
%                     'BackgroundColor', color, 'FaceAlpha', 0.5);
                xctr = round((x_max + x_min) / 2);
                xpos = xctr/xlim_hi;
                if xpos > 0.99; xpos = 0.99; end
                %y_shift = r_height / 40;
                %text(xctr, y_max - y_shift, label, 'Color', label_color, 'HorizontalAlignment', 'center');
                annotation('textbox', [xpos 0.82 0.1 0.1], 'String', label, 'HorizontalAlignment', 'center',...
                    'FontWeight', 'bold', 'BackgroundColor', 'w', 'FaceAlpha', 0.7,...
                    'LineWidth', 1, 'FitBoxToText', 'on');
            else
%                 rec = annotation('rectangle', dim, ...
%                     'LineStyle', 'none', ...
%                     'FaceColor', color, 'FaceAlpha', 0.5);
            end
            %fprintf("break;\n");
        end
        
        function fig_handle = plotThreshRanges(rnaspots_run, curve, y_lbl, yrange, figno)
            if isempty(rnaspots_run); fig_handle = []; return; end
            if isempty(curve); fig_handle = []; return; end
            if isempty(rnaspots_run.threshold_results); fig_handle = []; return; end

            thres = rnaspots_run.threshold_results;
            main_th = thres.threshold;
            curve_info_list = RNAThreshold.getAllCurveResultStructs(thres);
            curve_count = size(curve_info_list,2);
            mad_count = size(curve_info_list(1).med_suggested_threshold,2);
            alloc = 0;

            if thres.fit_ri_weight > 0.0
                alloc = alloc + curve_count;
            end
            if thres.madth_weight > 0.0
                alloc = alloc + (curve_count * mad_count);
            end
            if thres.fit_weight > 0.0
                alloc = alloc + curve_count;
            end

            all_thresh = NaN(alloc,1);
            i = 1;
            if thres.fit_ri_weight > 0.0
                for j = 1:curve_count
                    if ~isempty(curve_info_list(j).spline_fit)
                        all_thresh(i,1) = curve_info_list(j).spline_fit.rcurve_intr_x;
                        i = i + 1;
                    end
                end
            end
            if thres.madth_weight > 0.0
                for j = 1:curve_count
                    for k = 1:mad_count
                        all_thresh(i,1) = curve_info_list(j).med_suggested_threshold(k);
                        i = i + 1;
                    end
                end
            end
            if thres.fit_weight > 0.0
                for j = 1:curve_count
                    if ~isempty(curve_info_list(j).spline_fit)
                        all_thresh(i,1) = curve_info_list(j).spline_knot_x;
                        i = i + 1;
                    end
                end
            end

            th_avg = mean(all_thresh, 'all', 'omitnan');
            th_std = std(all_thresh, 0, 'all', 'omitnan');
            th_max = max(all_thresh,[],'all','omitnan');
            th_min = min(all_thresh,[],'all','omitnan');
            
            color_curve = [0.290, 0.290, 0.290];
            color_A = [0.567, 0.133, 0.133]; %#912222
            %color_B = [0.780, 0.118, 0.118]; %#
            color_B = [1.000, 0.718, 0.718]; %#ffb7b7
            color_C = [0.392, 0.020, 0.020]; %#640505

            fig_handle = figure(figno);
            clf;
            ax = axes;
            if ~isempty(yrange)
                ylim([0.0 1.0]);
            else
                y_hi = max(curve(:,2),[],'all','omitnan');
                ylim([0.0 y_hi]);
            end
            x_hi = max(curve(:,1), [], 'all', 'omitnan');
            xlim([0.0 x_hi]);
            hold on;
            RNAThreshold.drawVerticalBox(ax, th_avg - th_std, th_avg + th_std, color_B, 'Mean +- 1 StDev', color_C);
            plot(curve(:,1), curve(:,2), 'LineWidth',2,'Color',color_curve);
            xline(main_th, '--', 'Selected Threshold', 'Color', color_A,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            xline(th_min, ':', 'Minimum', 'Color', color_A,'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            xline(th_max, ':', 'Maximum', 'Color', color_A,'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            %RNAThreshold.drawVerticalLines(ax, main_th, [th_min th_max], color_A, color_A, 2);

            xlabel('Threshold');
            ylabel(y_lbl);
        end
        
    end
    
end