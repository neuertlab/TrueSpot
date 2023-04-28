%
%%

classdef SpotPlots
    methods (Static)

        function fig_handle = renderThreshVersusPlot(x, y, y_label, y_max, base_color, threshold, thresh_res, figno, existing_fig)

            if nargin < 3; y_label = ''; end
            if nargin < 4; y_max = 1.0; end
            if nargin < 5; base_color = [0.25 0.25 0.25]; end
            if nargin < 6; threshold = 0; end
            if nargin < 7; thresh_res = []; end
            if nargin < 8; figno = round(rand() * 10000); end
            if nargin < 9; existing_fig = []; end

            fig_handle = [];
            if isempty(x); return; end
            if isempty(y); return; end

            color_light = min(base_color + ((1.0 - base_color) .* 0.6),[1.0,1.0,1.0]);
            color_dark = max((base_color .* 0.5),[0,0,0]);

            if isempty(existing_fig)
                fig_handle = figure(figno);
                clf;
            else
                fig_handle = existing_fig;
                figure(fig_handle);
            end

            thr_min = 0;
            thr_max = 0;
            if ~isempty(thresh_res)
                score_list = double(RNAThreshold.getAllThresholdSuggestions(thresh_res));
                thr_mean = mean(score_list, 'all', 'omitnan');
                thr_std = std(score_list, 0, 'all', 'omitnan');
                thr_min = min(score_list, [], 'all', 'omitnan');
                thr_max = max(score_list, [], 'all', 'omitnan');
                x_min = thr_mean - thr_std;
                x_max = thr_mean + thr_std;

                r_width = x_max - x_min;
                rectangle('Position', [x_min 0 r_width y_max],...
                    'FaceColor', color_light, 'LineStyle', 'none');
                hold on;
            end

            plot(x,y,'LineWidth',2,'Color',base_color);
            hold on;

            if threshold > 0
                xline(threshold, '--', 'Selected Threshold', 'Color', color_dark,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end

            if ~isempty(thresh_res)
                %Min and max
                xline(thr_min, ':', 'Minimum', 'Color', color_dark,'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
                xline(thr_max, ':', 'Maximum', 'Color', color_dark,'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end

            ylim([0.0 y_max]);
            set(gca,'XTickLabel',[]);
            xlabel('Threshold (a.u.)');
            ylabel(y_label);
        end

        function fig_handle = renderLogSpotCountPlot(x, y, base_color, threshold, thresh_res, figno, existing_fig)
            y_max = max(y,[],'all','omitnan') + 0.5;
            fig_handle = SpotPlots.renderThreshVersusPlot(x, y, 'log10(# Spots Detected)', y_max, base_color, threshold, thresh_res, figno, existing_fig);
        end

        function fig_handle = renderFScorePlot(x, y, base_color, threshold, thresh_res, figno, existing_fig)
            fig_handle = SpotPlots.renderThreshVersusPlot(x, y, 'F-Score', 1.0, base_color, threshold, thresh_res, figno, existing_fig);
        end

        function fig_handle = renderPRPlot(recall, precision, base_color, figno, existing_fig)

            if nargin < 3; base_color = [0.25 0.25 0.25]; end
            if nargin < 4; figno = round(rand() * 10000); end
            if nargin < 5; existing_fig = []; end

            %color_light = base_color + ((1.0 - base_color) .* 0.5);
            color_dark = max((base_color .* 0.5),[0,0,0]);

            if isempty(existing_fig)
                fig_handle = figure(figno);
                clf;
            else
                fig_handle = existing_fig;
                figure(fig_handle);
            end

            plot(recall, precision, 'LineStyle', 'none', 'Marker', '.', 'MarkerEdgeColor', base_color);
            hold on;

            xline(0.95, ':', '95% Recall', 'Color', color_dark,'LineWidth',1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');

            ylim([0.0 1.0]);
            xlim([0.0 1.0]);

            xlabel('Recall');
            ylabel('Precision');
        end

    end
end