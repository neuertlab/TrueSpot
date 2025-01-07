%%
classdef ProbDistroPlots

    %%
    properties
        timePoints = [0:5:30];
        targets = {};

        xMax = 150;
        yMax = 1.0;
        topMargin = 0.050;
        leftMargin = 0.050;
        xSpace = 0.004;
        ySpace = 0.030;

        fszYLabel = 12;
        fszTicks = 11;
        fszXTitle = 14;

        binSize = 10;
        runSmoothing = false;
        smoothWindowSize = 3;
        timeUnit = 'min';

        data = {}; %2D mtx of data point structs

    end

    %%
    methods

        %%
        function obj = reallocateDataMtx(obj)
            if ~isempty(obj.timePoints)
                tpCount = size(obj.timePoints, 2);
            else
                obj.data = {};
            end

            if ~isempty(obj.targets)
                trgCount = size(obj.targets, 2);
            else
                obj.data = {};
            end

            obj.data = cell(trgCount, tpCount);
        end

        %%
        function [obj, figHandle] = render(obj, figHandle)
            figure(figHandle);
            hold on;

            if ~isempty(obj.timePoints)
                tpCount = size(obj.timePoints, 2);
            else
                return;
            end

            if ~isempty(obj.targets)
                trgCount = size(obj.targets, 2);
            else
                return;
            end

            ww = ((1.0 - obj.leftMargin) / tpCount) - (obj.xSpace * 1.1);
            hh = ((1.0 - obj.topMargin) / trgCount) - (obj.ySpace * 1.1);
            yy = 1.0 - hh - obj.topMargin;

            xtick_dist = round(obj.xMax ./ 4);
            xtick_vals = [0:xtick_dist:obj.xMax];
            xtickCount = size(xtick_vals, 2);
            xtick_vals(xtickCount) = obj.xMax;
            xtick_str_all = cell(1, xtickCount);
            for i = 1:xtickCount
                xtick_str_all{i} = num2str(xtick_vals(i));
            end
            xtick_str_trimmed = xtick_str_all;
            xtick_str_trimmed{1} = '';
            clear i

            %n = 1;
            subplot(trgCount, tpCount, 1);
            for y = 1:trgCount
                trgInfo = obj.targets{y};
                xx = obj.leftMargin;
                for x = 1:tpCount
                    %subplot(trgCount, tpCount, n);
                    subplot('Position', [xx yy ww hh]);
                    hold on;

                    dataGroup = obj.data{y, x};
                    if ~isempty(dataGroup)
                        repCount = size(dataGroup.reps, 2);
                        for r = 1:repCount
                            rcolor = trgInfo.colors(r,:);
                            lwidth = trgInfo.lineWidth(r);
                            lstyle = trgInfo.lineStyle{r};

                            xData = dataGroup.reps(r).x;
                            yData = dataGroup.reps(r).y;
                            
                            if ~isempty(xData)
                                if obj.runSmoothing
                                    yData = smoothdata(yData, 'movmean', obj.smoothWindowSize);
                                end

                                plot(xData, yData,...
                                    'Color', rcolor, 'LineWidth', lwidth, 'LineStyle', lstyle);
                            else
                                plot(0, 0,...
                                    'Color', rcolor, 'LineWidth', lwidth, 'LineStyle', lstyle);
                            end
                            clear rcolor lwidth lstyle
                        end
                    else
                        if x == tpCount
                            repCount = size(trgInfo.lineWidth, 2);
                            for r = 1:repCount
                                rcolor = trgInfo.colors(r,:);
                                lwidth = trgInfo.lineWidth(r);
                                lstyle = trgInfo.lineStyle{r};
                                plot(0, 0,...
                                    'Color', rcolor, 'LineWidth', lwidth, 'LineStyle', lstyle);
                                clear rcolor lwidth lstyle
                            end
                        end
                    end
                    
                    xlim([0 obj.xMax]);
                    ylim([0 obj.yMax]);

                    %Apply/remove titles or labels as appropriate
                    ax = gca;
                    ax.FontSize = obj.fszTicks;
                    if x > 1
                        set(gca,'YTickLabel',[]);
                    end

                    if y < trgCount
                        set(gca,'XTickLabel',[]);
                    else
                        %Reduce ticks shown to avoid crowding
                        set(gca,'XTick',xtick_vals);
                        if x == 1
                            set(gca,'XTickLabel',xtick_str_all);
                        else
                            set(gca,'XTickLabel',xtick_str_trimmed);
                        end
                    end

                    if y == 1
                        %Top titles
                        if x == 1
                            title([num2str(obj.timePoints(x)) ' ' obj.timeUnit], 'FontSize', obj.fszXTitle);
                        else
                            title(num2str(obj.timePoints(x)), 'FontSize', obj.fszXTitle);
                        end
                        
                    end

                    if x == 1
                        %Side titles
                        if ~isempty(trgInfo.subtitle)
                            ylabel({trgInfo.name; trgInfo.subtitle}, 'FontWeight', 'bold', 'FontSize', obj.fszYLabel);
                        else
                            ylabel(trgInfo.name, 'FontWeight', 'bold', 'FontSize', obj.fszYLabel);
                        end
                    end

                    if x == tpCount
                        %Legend
                        legend(trgInfo.repNames);
                    end

                    %n = n+1;
                    xx = xx + ww + obj.xSpace;
                end
                yy = yy - hh - obj.ySpace;
            end
        end

        %%
        function [obj, figHandle] = renderJointProbHeatmap(obj, figHandle, targetPairs)
            %Uses stored raw counts
            %Target pairs is a cell vector of joint pair structs
            %If at a target/timepoint/rep the number of cells does not
            %match between the two targets, will skip!

            %Combine replicates for each heatmap subplot?
            figure(figHandle);
            clf;

            if ~isempty(obj.timePoints)
                tpCount = size(obj.timePoints, 2);
            else
                return;
            end

            if ~isempty(targetPairs)
                trgPairCount = size(targetPairs, 2);
            else
                return;
            end

            hmYSpace = obj.ySpace * 1.5;
            hmXSpace = obj.xSpace * 0.3;
            hmLeft = obj.leftMargin * 0.6;
            ww = ((1.0 - (hmLeft * 2)) / tpCount) - (hmXSpace * 1.1);
            hh = ((1.0 - obj.topMargin) / trgPairCount) - (hmYSpace * 1.1);
            yy = 1.0 - hh - obj.topMargin;

            %Bin labels
            bin1 = [0:obj.binSize:(obj.xMax - obj.binSize)];
            binCount = size(bin1, 2);
            emptyLabels = cell(1, binCount);
            for i = 1:binCount; emptyLabels{i} = ' '; end
            edgeLabels_all = emptyLabels;
            edgeLabels_all{1} = '0';
            edgeLabels_all{binCount} = num2str(obj.xMax - obj.binSize);
            mididx = round(binCount./2);
            midval = bin1(mididx);
            edgeLabels_all{mididx} = num2str(midval);
            edgeLabels_trimmed = edgeLabels_all;
            edgeLabels_trimmed{1} = '';

            %subplot(trgPairCount, tpCount, 1);
            for y = 1:trgPairCount
                trgIndexX = targetPairs{y}.xTargetIndex;
                trgIndexY = targetPairs{y}.yTargetIndex;
                heatLog = targetPairs{y}.heatLogScale;
                trgInfoX = obj.targets{trgIndexX};
                trgInfoY = obj.targets{trgIndexY};
                xx = hmLeft;

                hmHandles = cell(1, tpCount);
                cMax = 0;
                for x = 1:tpCount
                    %subplot(trgCount, tpCount, n);
                    subplot('Position', [xx yy ww hh]);
                    %hold on;

                    dataGroupX = obj.data{trgIndexX, x};
                    dataGroupY = obj.data{trgIndexY, x};
                    [hmBins, ~] = ProbDistroPlots.heatmapBin(dataGroupX, dataGroupY, obj.xMax, obj.binSize, ~heatLog);
                    myMax = max(hmBins, [], 'all', 'omitnan');
                    if myMax > cMax; cMax = myMax; end
                    clear myMax

                    hold off;
                    hm = heatmap(bin1, bin1, hmBins);
                    hm.Colormap = turbo;
                    hm.CellLabelColor = 'none';
                    hm.GridVisible = 'off';
                    hm.FontSize = obj.fszTicks;
                    if heatLog
                        hm.ColorScaling = 'log';
                    else
                        hm.ColorLimits = [0.0 1.0];
                    end

                    if x < tpCount
                        hm.ColorbarVisible = 'off';
                    end

                    %https://www.mathworks.com/matlabcentral/answers/513971-how-can-i-modify-the-x-and-y-axis-tick-labels-for-a-heatmap-chart
                    if x == 1
                        hm.YDisplayLabels = edgeLabels_all;
%                         if y == 1
%                             hm.YDisplayLabels = edgeLabels_all;
%                         else
%                             hm.YDisplayLabels = edgeLabels_trimmed;
%                         end
                    else 
                        hm.YDisplayLabels = emptyLabels;
                    end
                    if y == trgPairCount
                        if x == 1
                            hm.XDisplayLabels = edgeLabels_all;
                        else
                            hm.XDisplayLabels = edgeLabels_trimmed;
                        end
                    else
                        hm.XDisplayLabels = emptyLabels;
                    end

                    if y == 1
                        %Top titles
                        if x == 1
                            title([num2str(obj.timePoints(x)) ' ' obj.timeUnit]);
                        else
                            title(num2str(obj.timePoints(x)));
                        end

                    end

                    if x == 1
                        %Side titles
                        if ~isempty(trgInfoY.subtitle)
                            ylabel([trgInfoY.name ' [' trgInfoY.subtitle ']']);
                        else
                            ylabel(trgInfoY.name);
                        end

                        if ~isempty(trgInfoX.subtitle)
                            xlabel([trgInfoX.name ' [' trgInfoX.subtitle ']']);
                        else
                            xlabel(trgInfoX.name);
                        end
                    end
                    
                    %n = n+1;
                    hmHandles{x} = hm;
                    xx = xx + ww + hmXSpace;
                end
                if heatLog
                    cLim = log(cMax);
                    for x = 1:tpCount
                        hmHandles{x}.ColorLimits = [0.0 cLim];
                    end
                end
                clear hmHandles

                yy = yy - hh - hmYSpace;
            end
        end

        %%
        %Removes all existing data!
        function obj = setTimePoints(obj, val)
            obj.timePoints = val;
            obj = obj.reallocateDataMtx();
        end

        %%
        function obj = loadRawCountSet(obj, rawCounts, targetIndex, timepointIndex, replicate)
            dataGroup = obj.data{targetIndex, timepointIndex};
            if isempty(dataGroup)
                dataGroup = ProbDistroPlots.genDataPointStruct(replicate);
            end

            myhisto = dataGroup.reps(replicate);
            cellCount = size(rawCounts, 2);
            binEdges = [0:obj.binSize:obj.xMax];
            [myhisto.y, myhisto.x] = histcounts(rawCounts, binEdges);
            myhisto.y = myhisto.y ./ cellCount;
            myhisto.x = myhisto.x(1:size(myhisto.y, 2));
            myhisto.rawCounts = rawCounts;

            dataGroup.reps(replicate) = myhisto;
            obj.data{targetIndex, timepointIndex} = dataGroup;
        end

    end

    %%
    methods (Static)
        
        %%
        function jointPairStruct = genJointPairStruct(xIdx, yIdx, useLog)
            if nargin < 1; xIdx = 0; end
            if nargin < 2; yIdx = 0; end
            if nargin < 3; useLog = false; end

            jointPairStruct = struct();
            jointPairStruct.xTargetIndex = xIdx;
            jointPairStruct.yTargetIndex = yIdx;
            jointPairStruct.heatLogScale = useLog;
        end

        %%
        function tstruct = genTargetInfoStruct(replicateCount)
            tstruct = struct();
            tstruct.name = 'Target';
            tstruct.subtitle = [];
            tstruct.colors = zeros(replicateCount, 3);
            tstruct.lineWidth = repmat(2, 1, replicateCount);

            tstruct.lineStyle = cell(1, replicateCount);
            for i = 1:replicateCount; tstruct.lineStyle{i} = '-'; end

            tstruct.repNames = cell(1, replicateCount);
            for i = 1:replicateCount; tstruct.repNames{i} = ['Rep. ' num2str(i)]; end
        end

        %%
        function dstruct = genDataPointStruct(replicateCount)
            dstruct = struct();
            dstruct.reps(1, replicateCount) = ProbDistroPlots.genHistStruct();
            for i = 1:replicateCount-1
                dstruct.reps(i) = ProbDistroPlots.genHistStruct();
            end
        end

        %%
        function dstruct = genHistStruct()
            dstruct = struct();
            dstruct.x = [];
            dstruct.y = [];
            dstruct.rawCounts = [];
        end

        %%
        function xMax = suggestXMax(rawCounts, roundTo)
            p95 = prctile(rawCounts, 95, 'all');
            xMax = ((p95 + (roundTo - 1)) ./ roundTo) .* roundTo;
        end

        %%
        function [hmBins, okay] = heatmapBin(dataGroupX, dataGroupY, xMax, binSize, normalize)
            bins1 = [0:binSize:xMax];
            bin1Count = size(bins1, 2) - 1;
            hmBins = zeros(bin1Count, bin1Count);
            okay = false;

            if isempty(dataGroupX); return; end
            if isempty(dataGroupY); return; end
            if ~isfield(dataGroupX, 'reps'); return; end
            if ~isfield(dataGroupY, 'reps'); return; end

            cellCount = 0;
            repCountX = size(dataGroupX.reps, 2);
            repCountY = size(dataGroupY.reps, 2);
            repCount = min(repCountX, repCountY);
            if repCount < 1; return; end
            for r = 1:repCount
                dataX = dataGroupX.reps(r);
                dataY = dataGroupY.reps(r);

                if isempty(dataX.rawCounts); continue; end
                if isempty(dataY.rawCounts); continue; end
                 
                ctX = size(dataX.rawCounts, 2);
                ctY = size(dataY.rawCounts, 2);

                if ctX ~= ctY; continue; end

                cellCount = cellCount + ctX;

                gtX = dataX.rawCounts >= bins1';
                binX = sum(gtX, 1); %Trim out of range. Don't bin.
                binX(binX > bin1Count) = NaN;
                gtY = dataY.rawCounts >= bins1';
                binY = sum(gtY, 1);
                binY(binY > bin1Count) = NaN;

                for i = 1:ctX
                    xx = binX(i);
                    yy = binY(i);
                    if isfinite(xx) & isfinite(yy)
                        hmBins(yy, xx) = hmBins(yy, xx) + 1;
                    end
                end
            end

            if cellCount < 1; return; end
            if normalize
                hmBins = hmBins ./ cellCount;
            end
            okay = true;
        end

    end

end