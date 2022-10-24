%
%%

classdef RNAThreshold
    
    methods (Static)
        
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
                param_struct.sample_spot_table = spot_table;
            else
                [rnaspots_run, param_struct.sample_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
            end
            
            if nargin >= 4 & ~isempty(ctrl_table)
                param_struct.control_spot_table = ctrl_table;
            else
                [rnaspots_run, param_struct.control_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Control();
            end
            
            param_struct.verbosity = verbosity;
            
            winmin = rnaspots_run.ttune_winsz_min;
            if winmin < 2; winmin = 2; end
            winincr = rnaspots_run.ttune_winsz_incr;
            if winincr < 1; winincr = 1; end
            winmax = rnaspots_run.ttune_winsz_max;
            if winmax < winmin; winmax = winmin+winincr; end
            param_struct.window_sizes = [winmin:winincr:winmax];
            
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
        
        function param_info = paramsFromSpotsrun(rnaspots_run)
            param_info = RNA_Threshold_Common.genEmptyThresholdParamStruct();
            
            winmin = rnaspots_run.ttune_winsz_min;
            if winmin < 2; winmin = 2; end
            winincr = rnaspots_run.ttune_winsz_incr;
            if winincr < 1; winincr = 1; end
            winmax = rnaspots_run.ttune_winsz_max;
            if winmax < winmin; winmax = winmin+winincr; end
            param_info.window_sizes = [winmin:winincr:winmax];
            
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
        
        function score_list = getAllMedThresholds(threshold_results)
            madf_count = 0;
            curve_count = 0;
            wincount = 0;
            
            if ~isempty(threshold_results.test_data)
                curve_count = curve_count + 1; 
                cres = threshold_results.test_data;
                if ~isempty(cres.med_suggested_threshold)
                	madf_count = size(cres.med_suggested_threshold,2);
                end
            end
            if ~isempty(threshold_results.test_diff)
                curve_count = curve_count + 1; 
                cres = threshold_results.test_diff;
                if madf_count < 1 & ~isempty(cres.med_suggested_threshold)
                	madf_count = size(cres.med_suggested_threshold,2);
                end
            end
            
            if isfield(threshold_results, 'test_winsc') & ~isempty(threshold_results.test_winsc)
                wincount = size(threshold_results.test_winsc, 2); 
                curve_count = curve_count + wincount;
                if madf_count < 1
                    cres = threshold_results.test_winsc(1,1);
                    if ~isempty(cres.med_suggested_threshold)
                        madf_count = size(cres.med_suggested_threshold,2);
                    end
                end
            end
            
            alloc = madf_count * curve_count;
            score_list = NaN(1, alloc);
            
            i = 1;
            if ~isempty(threshold_results.test_data)
                score_list(1,i:i+madf_count-1) = threshold_results.test_data.med_suggested_threshold(1,:);
                i = i + madf_count;
            end
            if ~isempty(threshold_results.test_diff)
                score_list(1,i:i+madf_count-1) = threshold_results.test_diff.med_suggested_threshold(1,:);
                i = i + madf_count;
            end
            if wincount > 0
                for j = 1:wincount
                    cres = threshold_results.test_winsc(1,j);
                    score_list(1,i:i+madf_count-1) = cres.med_suggested_threshold(1,:);
                    i = i + madf_count;
                end
            end
        end
        
        function score_list = getAllFitThresholds(threshold_results)
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
            
            score_list = NaN(1, curve_count);
            i = 1;
            if ~isempty(threshold_results.test_data)
                score_list(1,i) = threshold_results.test_data.spline_knot_x;
                i = i + 1;
            end
            if ~isempty(threshold_results.test_diff)
                score_list(1,i) = threshold_results.test_diff.spline_knot_x;
                i = i + 1;
            end
            if wincount > 0
                for j = 1:wincount
                    cres = threshold_results.test_winsc(1,j);
                    score_list(1,i) = cres.spline_knot_x;
                    i = i + 1;
                end
            end
        end

        function score_list = getAllRightISectThresholds(threshold_results)
            curve_count = 0;
            wincount = 0;
            
            if ~isempty(threshold_results.test_data)
                if ~isempty(threshold_results.test_data.spline_fit)
                    curve_count = curve_count + 1; 
                end
            end
            if ~isempty(threshold_results.test_diff)
                if ~isempty(threshold_results.test_diff.spline_fit)
                    curve_count = curve_count + 1; 
                end
            end
            
            if isfield(threshold_results, 'test_winsc') & ~isempty(threshold_results.test_winsc)
                wincount = size(threshold_results.test_winsc, 2); 
                curve_count = curve_count + wincount;
            end
            
            score_list = NaN(1, curve_count);
            i = 1;
            if ~isempty(threshold_results.test_data)
                if ~isempty(threshold_results.test_data.spline_fit)
                    score_list(1,i) = threshold_results.test_data.spline_fit.rcurve_intr_x;
                    i = i + 1;
                end
            end
            if ~isempty(threshold_results.test_diff)
                if ~isempty(threshold_results.test_diff.spline_fit)
                    score_list(1,i) = threshold_results.test_diff.spline_fit.rcurve_intr_x;
                    i = i + 1;
                end
            end
            if wincount > 0
                for j = 1:wincount
                    cres = threshold_results.test_winsc(1,j);
                    if ~isempty(cres.spline_fit)
                        score_list(1,i) = cres.spline_fit.rcurve_intr_x;
                        i = i + 1;
                    end
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
        
        function [plot_handle, line_x] = draw_thres_plot(ax, data_x, data_y, curveres, include_madrange, color_base, color_dark, color_light, is_multi)
            %Plot curve
            point_count = size(data_y,1);
            plot_handle = plot(data_x(1:point_count,1), data_y,'LineWidth',2,'Color',color_base);
            hold on;
                
            %Plot fit
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
            
            if thx <= 0
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
        
        function rec = drawVerticalBox(ax, x_min, x_max, color)
            %Should go under plot since I'm not sure how to change alpha.
            %get(ax,'YLim')
            y_lim = get(ax,'YLim');
            r_height = y_lim(1,2) - y_lim(1,1);
            r_width = x_max - x_min;
            
            rec = rectangle('Position', [x_min y_lim(1,1) r_width r_height],...
                'FaceColor', color, 'LineStyle', 'none');
            %fprintf("break;\n");
        end
        
        function drawVerticalLines(ax, primary_x, secondary_x, primary_color, secondary_color, primary_thickness)
            if ~isempty(primary_x)
                if isvector(primary_x)
                    vsize = size(primary_x,2);
                    for i = 1:vsize
                        line_x = primary_x(1,i);
                        line([line_x line_x], get(ax,'YLim'),'Color',primary_color,'LineStyle','--','LineWidth',primary_thickness);
                    end
                else
                    line_x = primary_x;
                    line([line_x line_x], get(ax,'YLim'),'Color',primary_color,'LineStyle','--','LineWidth',primary_thickness);
                end
            end
            
            if ~isempty(secondary_x)
                if isvector(secondary_x)
                    vsize = size(secondary_x,2);
                    for i = 1:vsize
                        line_x = secondary_x(1,i);
                        line([line_x line_x], get(ax,'YLim'),'Color',secondary_color,'LineStyle',':','LineWidth',1);
                    end
                else
                    line_x = secondary_x;
                    line([line_x line_x], get(ax,'YLim'),'Color',secondary_color,'LineStyle',':','LineWidth',1);
                end
            end
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
            color_B = [0.937, 0.627, 0.627]; %#efa0a0

            fig_handle = figure(figno);
            clf;
            ax = axes;
            if ~isempty(yrange)
                ylim([0.0 1.0]);
            else
                y_hi = max(curve(:,2),[],'all','omitnan');
                ylim([0.0 y_hi]);
            end
            hold on;
            RNAThreshold.drawVerticalBox(ax, th_avg - th_std, th_avg + th_std, color_B);
            plot(curve(:,1), curve(:,2), 'LineWidth',2,'Color',color_curve);
            RNAThreshold.drawVerticalLines(ax, main_th, [th_min th_max], color_A, color_A, 2);

            xlabel('Threshold');
            ylabel(y_lbl);
        end

        function fig_handle = resultPlotFScore(rnaspots_run, f_scores, include_madth, include_splth, figno)
            %Tertiary (filled rectangle) represents stdev of all threshold
            %candidates
            %Secondary (light outer lines) represents full range.
            %Primary (main darker line) represents called threshold
            
            %Stdev and range are recalculated based on the two bool params
            if isempty(rnaspots_run); fig_handle = []; return; end
            if isempty(f_scores); fig_handle = []; return; end
            if isempty(rnaspots_run.threshold_results); fig_handle = []; return; end
            
            thres = rnaspots_run.threshold_results;
            main_th = thres.threshold;
            curve_count = 0;
            wincount = 0;
            madf_count = -1;
            if ~isempty(thres.test_data)
                curve_count = curve_count + 1; 
                cres = thres.test_data;
                if ~isempty(cres.med_suggested_threshold)
                	madf_count = size(cres.med_suggested_threshold,2);
                end
            end
            if ~isempty(thres.test_diff) 
                curve_count = curve_count + 1; 
                if madf_count < 0
                    cres = thres.test_diff;
                    if ~isempty(cres.med_suggested_threshold)
                        madf_count = size(cres.med_suggested_threshold,2);
                    end
                end
            end
            if isfield(thres, 'test_winsc') & ~isempty(thres.test_winsc)
                wincount = size(thres.test_winsc, 2); 
                curve_count = curve_count + wincount;
                if madf_count < 0
                    cres = thres.test_winsc(1,1);
                    if ~isempty(cres.med_suggested_threshold)
                        madf_count = size(cres.med_suggested_threshold,2);
                    end
                end
            end
            if curve_count < 1; fig_handle = []; return; end
            
            alloc = 0;
            if include_madth
                if madf_count > 0
                    alloc = madf_count * curve_count;
                end
            end
            if include_splth
                alloc = alloc + curve_count;
            end
            
            i = 1;
            thvals = NaN(1, alloc);
            if ~isempty(thres.test_data)
                if include_madth & (madf_count > 0)
                    thvals(1,i:i+madf_count-1) = thres.test_data.med_suggested_threshold(1,:);
                    i = i + madf_count;
                end
                if include_splth
                    thvals(1,i) = thres.test_data.spline_knot_x;
                    i = i + 1;
                end
            end
            if ~isempty(thres.test_diff)
                if include_madth & (madf_count > 0)
                    thvals(1,i:i+madf_count-1) = thres.test_diff.med_suggested_threshold(1,:);
                    i = i + madf_count;
                end
                if include_splth
                    thvals(1,i) = thres.test_diff.spline_knot_x;
                    i = i + 1;
                end
            end
            if wincount > 0
                for j = 1:wincount
                    cres = thres.test_winsc(1,j);
                    if include_madth & (madf_count > 0)
                        thvals(1,i:i+madf_count-1) = cres.med_suggested_threshold(1,:);
                        i = i + madf_count;
                    end
                    if include_splth
                        thvals(1,i) = cres.spline_knot_x;
                        i = i + 1;
                    end
                end
            end
            th_avg = mean(thvals, 'all', 'omitnan');
            th_std = std(thvals, 0, 'all', 'omitnan');
            th_max = max(thvals,[],'all','omitnan');
            th_min = min(thvals,[],'all','omitnan');
            
            color_curve = [0.290, 0.290, 0.290];
            color_A = [0.567, 0.133, 0.133]; %#912222
            color_B = [0.937, 0.627, 0.627]; %#efa0a0
            
            fig_handle = figure(figno);
            clf;
            ax = axes;
            ylim([0.0 1.0]);
            hold on;
            RNAThreshold.drawVerticalBox(ax, th_avg - th_std, th_avg + th_std, color_B);
            plot(thres.x(:,1), f_scores(:,1), 'LineWidth',2,'Color',color_curve);
            RNAThreshold.drawVerticalLines(ax, main_th, [th_min th_max], color_A, color_A, 2);
%             if th_avg ~= main_th
%                 RNAThreshold.drawVerticalLines(ax, th_avg, [], color_A, color_A, 1);
%             end

            xlabel('Threshold');
            ylabel('F-Score');
            
        end
        
    end
    
end