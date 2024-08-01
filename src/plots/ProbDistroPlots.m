%%
classdef ProbDistroPlots

    %%
    properties
        timePoints = [0:5:30];
        targets = [];

        xMax = 150;
        yMax = 1.0;
        binSize = 10;
        timeUnit = 'min';

        data = []; %2D mtx of data point structs

    end

    %%
    methods

        %%
        function obj = reallocateDataMtx(obj)
            if ~isempty(obj.timePoints)
                tpCount = size(obj.timePoints, 2);
            else
                obj.data = [];
            end

            if ~isempty(obj.targets)
                trgCount = size(obj.targets, 2);
            else
                obj.data = [];
            end

            obj.data(trgCount, tpCount) = struct();
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

            n = 1;
            subplot(trgCount, tpCount, n);
            for yy = 1:trgCount
                trgInfo = obj.targets(yy);
                for xx = 1:tpCount
                    %subplot(yy, xx, n);
                    subplot('Position', n);
                    hold on;

                    dataGroup = obj.data(yy, xx);
                    if ~isempty(dataGroup)
                        repCount = size(dataGroup.reps, 2);
                        for r = 1:repCount
                            rcolor = trgInfo.colors(r);
                            lwidth = trgInfo.lineWidth(r);
                            lstyle = trgInfo.lineStyle(r);
                            plot(dataGroup.reps(r).x, dataGroup.reps(r).y,...
                                'Color', rcolor, 'LineWidth', lwidth, 'LineStyle', lstyle);
                            clear rcolor lwidth lstyle
                        end
                    end
                    
                    xlim(0, obj.xMax);
                    ylim(0, obj.yMax);

                    %Apply/remove titles or labels as appropriate
                    if xx > 1
                        set(gca,'YTickLabel',[]);
                    end
                    if yy < trgCount
                        set(gca,'XTickLabel',[]);
                    end
                    if yy == 1
                        %Top titles
                        if xx == 1
                            title([num2str(obj.timePoints(xx)) ' ' obj.timeUnit]);
                        else
                            title(num2str(obj.timePoints(xx)));
                        end
                        
                    end
                    if xx == 1
                        %Side titles
                        ylabel(trgInfo.name);
                    end
                end
            end

            
        end

        %%
        %Removes all existing data!
        function obj = setTimePoints(obj, val)
            obj.timePoints = val;
            obj = obj.reallocateDataMtx();
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
            dstruct.reps(replicateCount, 1) = genHistStruct();
            for i = 1:replicateCount-1
                dstruct.reps(i) = genHistStruct();
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