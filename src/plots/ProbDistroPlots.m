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

        inclReps = [];

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
            clf;
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

                xMaxUse = obj.xMax;
                yMaxUse = obj.yMax;
                if ~isnan(trgInfo.xMaxOverride)
                    xMaxUse = trgInfo.xMaxOverride;

                    if ~isnan(trgInfo.xTickDistOverride)
                        xtick_dist_use = trgInfo.xTickDistOverride;
                    else
                        xtick_dist_use = round(xMaxUse ./ 4);
                    end
                    xtick_vals_use = [0:xtick_dist_use:xMaxUse];
                    xtickCount_use = size(xtick_vals_use, 2);
                    xtick_vals_use(xtickCount_use) = xMaxUse;
                    xtick_str_all_use = cell(1, xtickCount_use);
                    for i = 1:xtickCount_use
                        xtick_str_all_use{i} = num2str(xtick_vals_use(i));
                    end
                    xtick_str_trimmed_use = xtick_str_all_use;
                    xtick_str_trimmed_use{1} = '';
                    clear i
                else
                    xtick_vals_use = xtick_vals;
                    xtick_str_all_use = xtick_str_all;
                    xtick_str_trimmed_use = xtick_str_trimmed;
                end

                xx = obj.leftMargin;
                for x = 1:tpCount
                    %subplot(trgCount, tpCount, n);
                    subplot('Position', [xx yy ww hh]);
                    hold on;

                    dataGroup = obj.data{y, x};
                    if ~isempty(dataGroup)
                        repCount = size(dataGroup.reps, 2);
                        for r = 1:repCount
                            if ~isempty(obj.inclReps)
                                if ~ismember(r, obj.inclReps)
                                    continue;
                                end
                            end

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
                                if ~isempty(obj.inclReps)
                                    if ~ismember(r, obj.inclReps)
                                        continue;
                                    end
                                end

                                rcolor = trgInfo.colors(r,:);
                                lwidth = trgInfo.lineWidth(r);
                                lstyle = trgInfo.lineStyle{r};
                                plot(0, 0,...
                                    'Color', rcolor, 'LineWidth', lwidth, 'LineStyle', lstyle);
                                clear rcolor lwidth lstyle
                            end
                        end
                    end
                    
                    xlim([0 xMaxUse]);
                    ylim([0 yMaxUse]);

                    %Apply/remove titles or labels as appropriate
                    ax = gca;
                    ax.FontSize = obj.fszTicks;
                    if x > 1
                        set(gca,'YTickLabel',[]);
                    end

                    if (y < trgCount) & isnan(trgInfo.xMaxOverride)
                        set(gca,'XTickLabel',[]);
                    else
                        %Reduce ticks shown to avoid crowding
                        set(gca,'XTick',xtick_vals_use);
                        if x == 1
                            set(gca,'XTickLabel',xtick_str_all_use);
                        else
                            set(gca,'XTickLabel',xtick_str_trimmed_use);
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

                    if ~isnan(trgInfo.binSizeOverride) & (x == 1)
                        yShift = yMaxUse ./ 10;
                        xShift = xMaxUse ./ 2;
                        text(xMaxUse - xShift, yMaxUse - yShift, ['Bin Size: ' num2str(trgInfo.binSizeOverride)]);
                        clear xShift yShift
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
        function [obj, figHandle] = renderJointProbHeatmapOnePair(obj, figHandle, targetPair)
            figure(figHandle);
            clf;

            if ~isempty(obj.timePoints)
                tpCount = size(obj.timePoints, 2);
            else
                return;
            end

            if isempty(targetPair); return; end

            %Determine replicate count
            xtrg = obj.targets{targetPair.xTargetIndex};
            ytrg = obj.targets{targetPair.yTargetIndex};
            replCount = min(size(xtrg.repNames, 2), size(ytrg.repNames, 2));

            hmYSpace = obj.ySpace * 1.5;
            hmXSpace = obj.xSpace * 0.3;
            hmLeft = obj.leftMargin;
            ww = ((1.0 - (hmLeft * 2)) / tpCount) - (hmXSpace * 1.1);
            hh = ((1.0 - obj.topMargin) / replCount) - (hmYSpace * 1.1);
            yy = 1.0 - hh - obj.topMargin;

            for y = 1:replCount
                xx = hmLeft;
                hmHandles = cell(1, tpCount);
                cMax = 0;
                for x = 1:tpCount
                    subplot('Position', [xx yy ww hh]);

                    dataGroupX = obj.data{targetPair.xTargetIndex, x};
                    dataGroupY = obj.data{targetPair.yTargetIndex, x};
                    [hmBins, xBins, yBins, ~] = ProbDistroPlots.heatmapBin(dataGroupX, dataGroupY, y, ...
                        targetPair.xMax, targetPair.yMax,...
                        targetPair.xBinSize, targetPair.yBinSize, ~targetPair.heatLogScale);
                    myMax = max(hmBins, [], 'all', 'omitnan');
                    if myMax > cMax; cMax = myMax; end
                    clear myMax

                    xBinCount = size(xBins, 2) - 1;
                    yBinCount = size(yBins, 2) - 1;
                    xBins = xBins(1:xBinCount);
                    yBins = yBins(1:yBinCount);

                    hold off;
                    hm = heatmap(xBins, yBins, hmBins);
                    hm.Colormap = turbo;
                    hm.CellLabelColor = 'none';
                    hm.GridVisible = 'off';
                    hm.FontSize = obj.fszTicks;
                    if targetPair.heatLogScale
                        hm.ColorScaling = 'log';
                    else
                        if isnan(targetPair.colorScaleCap)
                            hm.ColorLimits = [0.0 1.0];
                        else
                            hm.ColorLimits = [0.0 targetPair.colorScaleCap];
                        end
                    end

                    if x < tpCount
                        hm.ColorbarVisible = 'off';
                    end

                    %Generate labels
                    xlbl = cell(1, xBinCount);
                    for n = 1:xBinCount; xlbl{n} = ""; end
                    ylbl = cell(1, yBinCount);
                    for n = 1:yBinCount; ylbl{n} = ""; end

                    lblVal = 0;
                    if x > 1
                        lblVal = targetPair.xLblInterval;
                    end
                    while lblVal < targetPair.xMax
                        idx = find((xBins == lblVal), 1);
                        if ~isempty(idx)
                            xlbl{idx} = num2str(lblVal);
                        end
                        lblVal = lblVal + targetPair.xLblInterval;
                    end

                    if x == 1
                        lblVal = 0;
                        if y > 1
                            lblVal = targetPair.yLblInterval;
                        end
                        while lblVal < targetPair.yMax
                            idx = find((yBins == lblVal), 1);
                            if ~isempty(idx)
                                ylbl{idx} = num2str(lblVal);
                            end
                            lblVal = lblVal + targetPair.yLblInterval;
                        end
                    end

                    hm.XDisplayLabels = xlbl;
                    hm.YDisplayLabels = ylbl;

                    %https://www.mathworks.com/matlabcentral/answers/554917-heatmap-axis-labels-printing-vertically
                    hAx = hm.NodeChildren(3);
                    hAx.XAxis.TickLabelRotation=0;

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
                        if ~isempty(ytrg.repNames{y})
                            if ~isempty(ytrg.subtitle)
                                ylabel([ytrg.repNames{y}; ytrg.name ' [' ytrg.subtitle ']']);
                            else
                                ylabel([ytrg.repNames{y}; ytrg.name]);
                            end
                        else
                            if ~isempty(ytrg.subtitle)
                                ylabel([ytrg.name ' [' ytrg.subtitle ']']);
                            else
                                ylabel(ytrg.name);
                            end
                        end
                    end

                    if ~isempty(xtrg.subtitle)
                        xlabel([xtrg.name ' [' xtrg.subtitle ']']);
                    else
                        xlabel(xtrg.name);
                    end

                    hmHandles{x} = hm;
                    xx = xx + ww + hmXSpace;
                end
                if targetPair.heatLogScale
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

            trgInfo = obj.targets{targetIndex};
            xMaxUse = obj.xMax;
            binSizeUse = obj.binSize;
            if ~isnan(trgInfo.xMaxOverride)
                xMaxUse = trgInfo.xMaxOverride;
            end
            if ~isnan(trgInfo.binSizeOverride)
                binSizeUse = trgInfo.binSizeOverride;
            end

            myhisto = dataGroup.reps(replicate);
            cellCount = size(rawCounts, 2);
            binEdges = [0:binSizeUse:xMaxUse];
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
            jointPairStruct.colorScaleCap = NaN;

            jointPairStruct.xMax = 150;
            jointPairStruct.yMax = 150;
            jointPairStruct.xBinSize = 10;
            jointPairStruct.yBinSize = 10;
            jointPairStruct.xLblInterval = 50;
            jointPairStruct.yLblInterval = 50;
        end

        %%
        function tstruct = genTargetInfoStruct(replicateCount)
            tstruct = struct();
            tstruct.name = 'Target';
            tstruct.subtitle = [];
            tstruct.colors = zeros(replicateCount, 3);
            tstruct.lineWidth = repmat(2, 1, replicateCount);
            tstruct.xMaxOverride = NaN;
            tstruct.yMaxOverride = NaN;
            tstruct.binSizeOverride = NaN;
            tstruct.xTickDistOverride = NaN;

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
        function [hmBins, xBins, yBins, okay] = heatmapBin(dataGroupX, dataGroupY, repl, xMax, yMax, xBinSize, yBinSize, normalize)
            xBins = [0:xBinSize:xMax];
            yBins = [0:yBinSize:yMax];
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
            
            %Glue reps together and take double hist
            if isempty(repl) | (repl < 1) | (repl > repCount)
                for r = 1:repCount
                    cellCount = cellCount + size(dataGroupX.reps(r).rawCounts, 2);
                end
                X = NaN(1, cellCount);
                Y = NaN(1, cellCount);
                pos = 1;
                for r = 1:repCount
                    rCellCount = size(dataGroupX.reps(r).rawCounts, 2);
                    endPos = pos + rCellCount;
                    X(pos:endPos) = dataGroupX.reps(r).rawCounts(:);
                    Y(pos:endPos) = dataGroupY.reps(r).rawCounts(:);
                    pos = endPos + 1;
                end
            else
                X = dataGroupX.reps(repl).rawCounts(:);
                Y = dataGroupY.reps(repl).rawCounts(:);
            end

            cellCount = size(X ,1);
            
            [hmBins,xBins,yBins] = histcounts2(X,Y,xBins,yBins);
            hmBins = hmBins';

            if cellCount < 1; return; end
            if normalize
                hmBins = hmBins ./ cellCount;
            end
            okay = true;
        end

    end

end