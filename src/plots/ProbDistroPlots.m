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
        xSpace = 0.008;
        ySpace = 0.030;

        fszYLabel = 12;
        fszTicks = 11;
        fszXTitle = 14;

        binSize = 10;
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
                            lstyle = trgInfo.lineStyle(r);
                            plot(dataGroup.reps(r).x, dataGroup.reps(r).y,...
                                'Color', rcolor, 'LineWidth', lwidth, 'LineStyle', lstyle);
                            clear rcolor lwidth lstyle
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

                    %n = n+1;
                    xx = xx + ww + obj.xSpace;
                end
                yy = yy - hh - obj.ySpace;
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

            dataGroup.reps(replicate) = myhisto;
            obj.data{targetIndex, timepointIndex} = dataGroup;
        end

    end

    %%
    methods (Static)
        
        %%
        function tstruct = genTargetInfoStruct(replicateCount)
            tstruct = struct();
            tstruct.name = 'Target';
            tstruct.subtitle = [];
            tstruct.colors = zeros(replicateCount, 3);
            tstruct.lineWidth = ones(1, replicateCount);
            tstruct.lineStyle = repmat('-', 1, replicateCount);
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
        end

    end

end