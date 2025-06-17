%%
%
%Much of this code ported from RNA_Threshold_SpotSelector
classdef SpotCallVisualization
    
%%
%Modes:
%   Results View
%       Draws spot circles for all spots in threshold & z range from
%       results calltable. If in single slice mode, spots on current z
%       plane are drawn in red. Spots below are in yellow, above in green.
%       
%       Basically, this is for assessing a call set.
%
%   Reference Select
%       Similar to results view mode, but there is no threshold value and
%       the table is read/write. 
%       Will take a new approach for z-snapping, maybe run a local
%       imregionalmax around a manual selection? That way don't need
%       previous computer call set.
%
%   Results Compare (Callset vs. Callset)
%       Compare two tool callsets. Calls for one will have spots drawn in
%       red, the other blue with magenta being agreements and cyan being
%       xy within z rad agreements.
%   
%   Results Compare (Callset vs. Reference)
%       Compare a tool callset to a reference callset. True positives are
%       shown in green, false positives in red, and false negatives in
%       yellow. (Note this is a different color scheme from
%       RNA_Threshold_SpotSelector, mostly to avoid confusion).
%
%Results from other tools can be visualized as long as they have been
%imported into a callTable formatted like TrueSpot's.

    %%
    properties
        callTable; %TrueSpot output call table
        callTableCompare;
        referenceTable; %Just an xyz matrix in a struct like the queues. Don't need a table.

        zMode; %Max proj(1) or single slice(2)?
        zNearRad; %In single slice mode, how many slices up and down are calls considered "near" and drawn in brighter colors?
        zShowRad; %In single slice mode, how many slices up and down to see calls from 
        zMIPColor; %Bool. If set, spot circle color is slightly lighter/darker dep. on z position.
        workRegion; %Outside this region, pretend spots don't exist.

        matchRad_xy;
        matchRad_z;
        matchMinTh;

        visCommon; %Instance of VisCommon to get color tables from

        %%Internal scratch
        refAddQueue;
        refRemoveQueue;
    end

    %%
    methods

        function obj = initializeMe(obj)
            obj.zMode = 1;
            obj.zNearRad = 2;
            obj.zShowRad = 7;
            obj.zMIPColor = false;
            obj.matchRad_xy = 4;
            obj.matchRad_z = 2;
            obj.matchMinTh = 20;

            obj.workRegion = struct();
            obj.workRegion.x_min = 0;
            obj.workRegion.x_max = 0;
            obj.workRegion.y_min = 0;
            obj.workRegion.y_max = 0;
            obj.workRegion.z_min = 0;
            obj.workRegion.z_max = 0;

            obj.refRemoveQueue = []; %Bool array for each row in ref data

            obj.refAddQueue = struct();
            obj.refAddQueue.queue = []; %Matrix w/ cols x,y,z,flag
            obj.refAddQueue.capacity = 0;
            obj.refAddQueue.used = 0;

            obj = obj.setRefSetEmpty();

            obj.callTableCompare = struct();
            obj.callTableCompare.callTable = [];
            obj.callTableCompare.boolAOnly = [];
            obj.callTableCompare.boolBOnly = [];
            obj.callTableCompare.boolFullMatchA = [];
            obj.callTableCompare.boolXYMatchA = [];
            obj.callTableCompare.checkedThA = 0;
            obj.callTableCompare.checkedThB = 0;

            obj.visCommon = VisCommon; %Don't initialize.
        end

        function [figHandle, okay] = drawResultsBasic(obj, figHandle, thVal, z)
            if nargin < 4; z = 0; end
            okay = false;
            if isempty(obj); return; end
            %if isempty(figHandle); return; end

            if obj.zMode == 1
                %MaxProj
                rowBool = SpotCallVisualization.filterCallsToValidRange(obj.callTable, obj.workRegion, 0, 0, thVal);
                if nnz(rowBool) > 0
                    xx = obj.callTable{rowBool, 'isnap_x'};
                    yy = obj.callTable{rowBool, 'isnap_y'};

                    if ~isempty(figHandle); figure(figHandle); end
                    hold on;
                    
                    if obj.zMIPColor & ~isempty(obj.visCommon) & obj.visCommon.isInitialized
                        zz = obj.callTable{rowBool, 'isnap_z'};
                        zMax = max(zz, [], 'all', 'omitnan');
                        obj.visCommon.idims.z = max(obj.visCommon.idims.z, zMax);
                        if size(obj.visCommon.colortbl_red, 1) < zMax
                            obj = obj.generateColorTables(obj.visCommon.idims.z + 1);
                        end
                        clear zMax
                        zc = obj.visCommon.idims.z/2;

                        for z = 1:obj.visCommon.idims.z
                            rowBoolz = (zz == z);
                            if nnz(rowBoolz) > 0
                                zdist = min(size(obj.visCommon.colortbl_red, 1), (abs(zc - z) * 2) + 1);
                                plot(xx(rowBoolz), yy(rowBoolz),...
                                'LineStyle','none','Marker','o',...
                                'MarkerEdgeColor',obj.visCommon.colortbl_red(zdist, :),'markersize',10);
                            end
                        end

                    else
                        plot(xx, yy,'LineStyle','none','Marker','o','MarkerEdgeColor','r','markersize',10);
                    end
                end
                okay = true;
            elseif obj.zMode == 2
                rowBool = SpotCallVisualization.filterCallsToValidRange(obj.callTable,...
                    obj.workRegion, max((z - obj.zShowRad), 1), z + obj.zShowRad, thVal);

                if nnz(rowBool) > 0
                    darkYellow = [0.60, 0.50, 0.0];
                    brightYellow = [1.0, 1.0, 0.0];
                    darkGreen = [0.0, 0.50, 0.20];
                    brightGreen = [0.0, 1.0, 0.0];

                    if ~isempty(obj.visCommon)
                        darkYellow = obj.visCommon.colorBelowFar;
                        brightYellow = obj.visCommon.colorBelowNear;
                        darkGreen = obj.visCommon.colorAboveFar;
                        brightGreen = obj.visCommon.colorAboveNear;
                    end

                    if ~isempty(figHandle); figure(figHandle); end
                    hold on;

                    groupBool = and(rowBool, obj.callTable{:, 'isnap_z'} < (z - obj.zNearRad));
                    if nnz(groupBool) > 0
                        plot(obj.callTable{groupBool, 'isnap_x'}, obj.callTable{groupBool, 'isnap_y'},...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',darkYellow,'markersize',10);
                    end

                    groupBool = and(rowBool, obj.callTable{:, 'isnap_z'} > (z + obj.zNearRad));
                    if nnz(groupBool) > 0
                        plot(obj.callTable{groupBool, 'isnap_x'}, obj.callTable{groupBool, 'isnap_y'},...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',darkGreen,'markersize',10);
                    end

                    groupBool = and(rowBool, obj.callTable{:, 'isnap_z'} >= (z - obj.zNearRad));
                    groupBool = and(groupBool, obj.callTable{:, 'isnap_z'} < z);
                    if nnz(groupBool) > 0
                        plot(obj.callTable{groupBool, 'isnap_x'}, obj.callTable{groupBool, 'isnap_y'},...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',brightYellow,'markersize',10);
                    end

                    groupBool = and(rowBool, obj.callTable{:, 'isnap_z'} <= (z + obj.zNearRad));
                    groupBool = and(groupBool, obj.callTable{:, 'isnap_z'} > z);
                    if nnz(groupBool) > 0
                        plot(obj.callTable{groupBool, 'isnap_x'}, obj.callTable{groupBool, 'isnap_y'},...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',brightGreen,'markersize',10);
                    end

                    groupBool = and(rowBool, obj.callTable{:, 'isnap_z'} == z);
                    if nnz(groupBool) > 0
                        plot(obj.callTable{groupBool, 'isnap_x'}, obj.callTable{groupBool, 'isnap_y'},...
                            'LineStyle','none','Marker','o','MarkerEdgeColor','r','markersize',10);
                    end
                end
                okay = true;
            end
        end

        function [figHandle, okay] = drawResultsTTCompare(obj, figHandle, thValA, thValB, z)
            okay = false;
            if isempty(obj); return; end
            %if isempty(figHandle); return; end
            if isempty(obj.callTable); return; end
            if isempty(obj.callTableCompare); return; end
            if isempty(obj.callTableCompare.callTable); return; end

            if (obj.callTableCompare.checkedThA ~= thValA) | (obj.callTableCompare.checkedThB ~= thValB)
                obj = obj.matchCompareSpots(thValA, thValB);
            end

            tableA = obj.callTable;
            tableB = obj.callTableCompare.callTable;

            if obj.zMode == 1
                %Draw all.
                rowBool1 = SpotCallVisualization.filterCallsToValidRange(tableA, obj.workRegion, 0, 0, thValA);
                rowBool2 = SpotCallVisualization.filterCallsToValidRange(tableB, obj.workRegion, 0, 0, thValB);
            elseif obj.zMode == 2
                %Only draw for this slice.
                rowBool1 = SpotCallVisualization.filterCallsToValidRange(tableA, obj.workRegion, z, z, thValA);
                rowBool2 = SpotCallVisualization.filterCallsToValidRange(tableB, obj.workRegion, z, z, thValB);
            else
                return;
            end

            if ~isempty(figHandle); figure(figHandle); end
            hold on;
            groupBool = and(rowBool2, obj.callTableCompare.boolBOnly);
            if nnz(groupBool) > 0
                plot(tableB{groupBool, 'isnap_x'}, tableB{groupBool, 'isnap_y'},...
                    'LineStyle','none','Marker','o','MarkerEdgeColor','blue','markersize',10);
            end

            groupBool = and(rowBool1, obj.callTableCompare.boolAOnly);
            if nnz(groupBool) > 0
                plot(tableA{groupBool, 'isnap_x'}, tableA{groupBool, 'isnap_y'},...
                    'LineStyle','none','Marker','o','MarkerEdgeColor','red','markersize',10);
            end

            groupBool = and(rowBool1, obj.callTableCompare.boolXYMatchA);
            if nnz(groupBool) > 0
                plot(tableA{groupBool, 'isnap_x'}, tableA{groupBool, 'isnap_y'},...
                    'LineStyle','none','Marker','o','MarkerEdgeColor','cyan','markersize',10);
            end

            groupBool = and(rowBool1, obj.callTableCompare.boolFullMatchA);
            if nnz(groupBool) > 0
                plot(tableA{groupBool, 'isnap_x'}, tableA{groupBool, 'isnap_y'},...
                    'LineStyle','none','Marker','o','MarkerEdgeColor','magenta','markersize',10);
            end
            okay = true;
        end

        function [figHandle, okay] = drawResultsTRCompare(obj, figHandle, thVal, z)
            okay = false;
            if isempty(obj); return; end
            %if isempty(figHandle); return; end
            if isempty(obj.callTable); return; end
            if isempty(obj.referenceTable); return; end
            if isempty(obj.referenceTable.data); return; end

            if ~obj.referenceTable.checkedVs
                obj = obj.matchReferenceSpots();
            end

            tableC = obj.callTable;
            tableR = obj.referenceTable.data(1:obj.referenceTable.used,:);

            if obj.zMode == 1
                rowBoolCall = SpotCallVisualization.filterCallsToValidRange(tableC, obj.workRegion, 0, 0, thVal);
                rowBoolRef = SpotCallVisualization.filterCallsToValidRangeMtx(tableR, obj.workRegion, 0, 0);
            elseif obj.zMode == 2
                rowBoolCall = SpotCallVisualization.filterCallsToValidRange(obj.callTable, obj.workRegion, z, z, thVal);
                rowBoolRef = SpotCallVisualization.filterCallsToValidRangeMtx(tableR, obj.workRegion, z, z);
            else
                return;
            end

            if ~isempty(figHandle); figure(figHandle); end
            hold on;

            %True positives
            groupBool = and(rowBoolCall, obj.referenceTable.truePos);
            if nnz(groupBool) > 0
                plot(tableC{groupBool, 'isnap_x'}, tableC{groupBool, 'isnap_y'},...
                    'LineStyle','none','Marker','o','MarkerEdgeColor','green','markersize',10);
            end

            %False positives
            groupBool = and(rowBoolCall, ~obj.referenceTable.truePos);
            if nnz(groupBool) > 0
                plot(tableC{groupBool, 'isnap_x'}, tableC{groupBool, 'isnap_y'},...
                    'LineStyle','none','Marker','o','MarkerEdgeColor','red','markersize',10);
            end

            %False negatives
            groupBool = and(rowBoolRef, obj.referenceTable.falseNeg);
            if nnz(groupBool) > 0
                plot(tableR(groupBool, 1), tableR(groupBool, 2),...
                    'LineStyle','none','Marker','o','MarkerEdgeColor','yellow','markersize',10);
            end

            okay = true;
        end

        function [figHandle, okay] = drawReferenceSet(obj, figHandle, z)
            okay = false;
            if isempty(obj); return; end
            %if isempty(figHandle); return; end
            if isempty(obj.referenceTable); return; end
            if isempty(obj.referenceTable.data); return; end

            tableR = obj.referenceTable.data(1:obj.referenceTable.used,:);

            if obj.zMode == 1
                rowBoolRef = SpotCallVisualization.filterCallsToValidRangeMtx(tableR, obj.workRegion, 0, 0);

                if nnz(rowBoolRef) > 0
                    if ~isempty(figHandle); figure(figHandle); end
                    hold on;
                    plot(tableR(rowBoolRef, 1), tableR(rowBoolRef, 2),'LineStyle','none','Marker','o','MarkerEdgeColor','r','markersize',10);
                end
            elseif obj.zMode == 2
                %Draw neighboring slice calls as with call only
                rowBoolRef = SpotCallVisualization.filterCallsToValidRangeMtx(tableR,...
                    obj.workRegion, max((z - obj.zShowRad), 1), z + obj.zShowRad);

                if nnz(rowBoolRef) > 0
                    darkYellow = [0.60, 0.50, 0.0];
                    brightYellow = [1.0, 1.0, 0.0];
                    darkGreen = [0.0, 0.50, 0.20];
                    brightGreen = [0.0, 1.0, 0.0];

                    if ~isempty(obj.visCommon)
                        darkYellow = obj.visCommon.colorBelowFar;
                        brightYellow = obj.visCommon.colorBelowNear;
                        darkGreen = obj.visCommon.colorAboveFar;
                        brightGreen = obj.visCommon.colorAboveNear;
                    end

                    if ~isempty(figHandle); figure(figHandle); end
                    hold on;

                    groupBool = and(rowBoolRef, tableR(:, 3) < (z - obj.zNearRad));
                    if nnz(groupBool) > 0
                        plot(tableR(groupBool, 1), tableR(groupBool, 2),...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',darkYellow,'markersize',10);
                    end

                    groupBool = and(rowBoolRef, tableR(:, 3) > (z + obj.zNearRad));
                    if nnz(groupBool) > 0
                        plot(tableR(groupBool, 1), tableR(groupBool, 2),...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',darkGreen,'markersize',10);
                    end

                    groupBool = and(rowBoolRef, tableR(:, 3) >= (z - obj.zNearRad));
                    groupBool = and(groupBool, tableR(:, 3) < z);
                    if nnz(groupBool) > 0
                        plot(tableR(groupBool, 1), tableR(groupBool, 2),...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',brightYellow,'markersize',10);
                    end

                    groupBool = and(rowBoolRef, tableR(:, 3) <= (z + obj.zNearRad));
                    groupBool = and(groupBool, tableR(:, 3) > z);
                    if nnz(groupBool) > 0
                        plot(tableR(groupBool, 1), tableR(groupBool, 2),...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',brightGreen,'markersize',10);
                    end

                    groupBool = and(rowBoolRef, tableR(:, 3) == z);
                    if nnz(groupBool) > 0
                        plot(tableR(groupBool, 1), tableR(groupBool, 2),...
                            'LineStyle','none','Marker','o','MarkerEdgeColor','red','markersize',10);
                    end
                end
                okay = true;
            end
        end

        function [figHandle, okay] = drawResultsBasic3D(obj, figHandle, thVal)
            %TODO
            okay = false;
        end

        function [figHandle, okay] = drawResultsTTCompare3D(obj, figHandle, thValA, thValB)
            %TODO
            okay = false;
        end

        function [figHandle, okay] = drawResultsTRCompare3D(obj, figHandle, thVal)
            %TODO
            okay = false;
        end

        function [figHandle, okay] = drawReferenceSet3D(obj, figHandle)
            %TODO
            okay = false;
        end

        function obj = setRefSetEmpty(obj)
            obj.referenceTable = struct();
            obj.referenceTable.data = []; %Matrix w/ cols x,y,z
            obj.referenceTable.capacity = 0;
            obj.referenceTable.used = 0;
            obj.referenceTable.truePos = []; %Bool vec for call table rows
            obj.referenceTable.falseNeg = []; %Bool vec for ref table rows
            obj.referenceTable.checkedVs = false;
        end

        function obj = expandAddQueue(obj)
            newCap = obj.refAddQueue.capacity + 256;
            newTable = zeros(newCap, 4);
            newTable(1:obj.refAddQueue.used, :) = obj.refAddQueue.queue(1:obj.refAddQueue.used, :);

            obj.refAddQueue.capacity = newCap;
            obj.refAddQueue.queue = newTable;
        end

        function [obj, figHandle] = onRefModeClick(obj, figHandle, x, y, z, allowAdd, allowRemove)
            %Draws new circle and updates queue
            if nargin < 5; z = 0; end
            if nargin < 6; allowAdd = true; end
            if nargin < 7; allowRemove = true; end
            if isempty(figHandle); return; end
            if isempty(obj.referenceTable); return; end

            %Get idims
            X = 1024;
            Y = 1024;
            Z = 1;
            if ~isempty(obj.visCommon) & obj.visCommon.isInitialized
                X = obj.visCommon.idims.x;
                Y = obj.visCommon.idims.y;
                Z = obj.visCommon.idims.z;
            else
                if ~isempty(obj.callTable)
                    X = max(obj.callTable{:, 'isnap_x'}, [], 'all', 'omitnan');
                    Y = max(obj.callTable{:, 'isnap_y'}, [], 'all', 'omitnan');
                    Z = max(obj.callTable{:, 'isnap_z'}, [], 'all', 'omitnan');
                elseif ~isempty(obj.referenceTable.data)
                    X = max(obj.referenceTable.data(:,1), [], 'all', 'omitnan');
                    Y = max(obj.referenceTable.data(:,2), [], 'all', 'omitnan');
                    Z = max(obj.referenceTable.data(:,3), [], 'all', 'omitnan');
                end
            end

            %Click radius, based on figure scale
            cRad = SpotCallVisualization.calculateClickRadius(figHandle, X, Y);
            
            %Find existing spots within that radius - both in queues and in
            %ref table
            %Also, ignore anything outside mask, if active
            searchRegion = struct();
            searchRegion.x_min = max(1, x - cRad);
            searchRegion.x_max = min(X, x + cRad);
            searchRegion.y_min = max(1, y - cRad);
            searchRegion.y_max = min(Y, y + cRad);

            if z == 0 | (obj.zMode == 1)
                searchRegion.z_min = 1;
                searchRegion.z_max = Z;
            else
                searchRegion.z_min = z;
                searchRegion.z_max = z;
            end

            if ~isempty(obj.workRegion)
                if obj.workRegion.x_min > 0
                    searchRegion.x_min = max(searchRegion.x_min, obj.workRegion.x_min);
                end
                if obj.workRegion.x_max > 0
                    searchRegion.x_max = min(searchRegion.x_max, obj.workRegion.x_max);
                end
                if obj.workRegion.y_min > 0
                    searchRegion.y_min = max(searchRegion.y_min, obj.workRegion.y_min);
                end
                if obj.workRegion.y_max > 0
                    searchRegion.y_max = min(searchRegion.y_max, obj.workRegion.y_max);
                end
                if obj.workRegion.z_min > 0
                    searchRegion.z_min = max(searchRegion.z_min, obj.workRegion.z_min);
                end
                if obj.workRegion.z_max > 0
                    searchRegion.z_max = min(searchRegion.z_max, obj.workRegion.z_max);
                end
            end

            figure(figHandle);
            hold on;

            %Check existing table
            noHits = true;
            if ~isempty(obj.referenceTable.data)
                if isempty(obj.refRemoveQueue)
                    obj.refRemoveQueue = false(size(obj.referenceTable.data, 1), 1);
                end

                rowBool = SpotCallVisualization.filterCallsToValidRangeMtx(...
                    obj.referenceTable.data, searchRegion);
                if nnz(rowBool) > 0
                    noHits = false;
                    if allowRemove
                        obj.refRemoveQueue(rowBool) = ~obj.refRemoveQueue(rowBool);
                        %Draw magenta
                        drawBool = and(~obj.refRemoveQueue, rowBool);
                        if nnz(drawBool) > 0
                            xx = obj.referenceTable.data(drawBool, 1);
                            yy = obj.referenceTable.data(drawBool, 2);
                            plot(xx, yy,'LineStyle','none','Marker','o','MarkerEdgeColor','magenta','markersize',10);
                        end

                        %Draw white
                        drawBool = and(obj.refRemoveQueue, rowBool);
                        if nnz(drawBool) > 0
                            xx = obj.referenceTable.data(drawBool, 1);
                            yy = obj.referenceTable.data(drawBool, 2);
                            plot(xx, yy,'LineStyle','none','Marker','o','MarkerEdgeColor','white','markersize',10);
                        end
                    end
                end
            end

            %Check add queue (Always allow)
            if ~isempty(obj.refAddQueue.queue) & (obj.refAddQueue.used > 0)
                rowBool = SpotCallVisualization.filterCallsToValidRangeMtx(...
                    obj.refAddQueue.queue, searchRegion);
                rowBool((obj.refAddQueue.used + 1):size(obj.refAddQueue, 1)) = false;

                if nnz(rowBool) > 0
                    noHits = false;
                    obj.refAddQueue.queue(rowBool) = ~obj.refAddQueue.queue(rowBool);

                    %Draw cyan
                    drawBool = and(obj.refAddQueue.queue(:,4), rowBool);
                    if nnz(drawBool) > 0
                        xx = obj.referenceTable.data(drawBool, 1);
                        yy = obj.referenceTable.data(drawBool, 2);
                        plot(xx, yy,'LineStyle','none','Marker','o','MarkerEdgeColor','cyan','markersize',10);
                    end

                    %Draw white
                    drawBool = and(~obj.refAddQueue.queue(:,4), rowBool);
                    if nnz(drawBool) > 0
                        xx = obj.referenceTable.data(drawBool, 1);
                        yy = obj.referenceTable.data(drawBool, 2);
                        plot(xx, yy,'LineStyle','none','Marker','o','MarkerEdgeColor','white','markersize',10);
                    end
                end
            end

            %Add click target if there are no matches
            if noHits & allowAdd
                if isempty(obj.refAddQueue.queue)
                    obj.refAddQueue.queue = zeros(256, 4);
                    obj.refAddQueue.capacity = 256;
                    obj.refAddQueue.used = 0;
                end

                if obj.refAddQueue.used >= obj.refAddQueue.capacity
                    obj = obj.expandAddQueue();
                end

                if z == 0
                    z = round(Z ./ 2);
                end
                newIndex = obj.refAddQueue.used + 1;
                obj.refAddQueue.queue(newIndex, 1) = round(x);
                obj.refAddQueue.queue(newIndex, 2) = round(y);
                obj.refAddQueue.queue(newIndex, 3) = z;
                obj.refAddQueue.queue(newIndex, 4) = true;
                obj.refAddQueue.used = newIndex;

                plot(x, y,'LineStyle','none','Marker','o','MarkerEdgeColor','cyan','markersize',10);
            end

        end

        function obj = onRefModeAccept(obj)
            %updates ref table to reflect click queue
            %Does not redraw anything though. The draw function must be
            %   called again!

            %First, run removes
            if ~isempty(obj.refRemoveQueue) & ~isempty(obj.referenceTable.data)
                if nnz(obj.refRemoveQueue) > 0
                    rCount = nnz(obj.refRemoveQueue);
                    nowCount = size(obj.referenceTable.data,1);
                    obj.referenceTable.data = obj.referenceTable.data(~obj.refRemoveQueue, :);
                    obj.referenceTable.capacity = size(obj.referenceTable.data,1);
                    obj.referenceTable.used = nowCount - rCount;
                end
            end
            obj.refRemoveQueue = [];

            %Then adds
            if ~isempty(obj.refAddQueue.queue) & (obj.refAddQueue.used > 0)
                %Whittle down to those to actually add.
                addActual = obj.refAddQueue.queue(1:obj.refAddQueue.used, :);
                rowbool = logical(addActual(:,4));
                if nnz(rowbool) > 0
                    addActual = addActual(rowbool, 1:3);
                    obj.referenceTable.data = [obj.referenceTable.data; addActual];
                    obj.referenceTable.capacity = size(obj.referenceTable.data, 1);
                    obj.referenceTable.used = obj.referenceTable.capacity;
                end
            end

            obj.refAddQueue.queue = [];
            obj.refAddQueue.capacity = 0;
            obj.refAddQueue.used = 0;

        end

        function [obj, snappedCount, removedCount] = refModeSnap(obj, filteredImg, zOnly, anyZ, removeUnsnapped, verbose)
            if nargin < 3; zOnly = true; end
            if nargin < 4; anyZ = false; end
            if nargin < 5; removeUnsnapped = false; end
            if nargin < 6; verbose = false; end
            snappedCount = 0;
            removedCount = 0;

            %Runs spot snap (looking for local maxima)
            if isempty(obj); return; end
            if isempty(filteredImg); return; end
            if isempty(obj.referenceTable); return; end
            if isempty(obj.referenceTable.data); return; end

            refCount = obj.referenceTable.used;
            keepRow = true(refCount, 1);

            X = size(filteredImg, 2);
            Y = size(filteredImg, 1);
            Z = size(filteredImg, 3);

            %I don't know any way around looping offhand. Maybe fix this
            %later.
            for i = 1:refCount
                x = obj.referenceTable.data(i,1);
                y = obj.referenceTable.data(i,2);
                z = obj.referenceTable.data(i,3);

                %Calculate search box...
                if zOnly
                    x_min = x; x_max = x;
                    y_min = y; y_max = y;
                else
                    x_min = max(x - obj.matchRad_xy, 1);
                    x_max = min(x + obj.matchRad_xy, X);
                    y_min = max(y - obj.matchRad_xy, 1);
                    y_max = min(y + obj.matchRad_xy, Y);
                end

                if anyZ
                    z_min = 1; z_max = Z;
                else
                    z_min = max(z - obj.matchRad_z, 1);
                    z_max = min(z + obj.matchRad_z, Z);
                end

                ibox = filteredImg(y_min:y_max, x_min:x_max, z_min:z_max);
                [mVal, mIdx] = max(ibox, [], 'all', 'omitnan');
                if removeUnsnapped & (mVal < obj.matchMinTh)
                    keepRow(i, 1) = false;
                else
                    if verbose; fprintf('[%d, %d, %d] snapped to ', x, y, z); end
                    [yy, xx, zz] = ind2sub(size(ibox), mIdx);
                    x = xx + x_min - 1;
                    y = yy + y_min - 1;
                    z = zz + z_min - 1;
                    filteredImg(y,x,z) = 0; %Zero out in local copy so don't find again.
                    obj.referenceTable.data(i,1) = x;
                    obj.referenceTable.data(i,2) = y;
                    obj.referenceTable.data(i,3) = z;
                    snappedCount = snappedCount + 1;
                    if verbose; fprintf('[%d, %d, %d]\n', x, y, z); end
                end
            end

            removedCount = nnz(~keepRow);
            if removedCount > 0
                obj.referenceTable.data = obj.referenceTable.data(keepRow, :);
                obj.referenceTable.used = size(obj.referenceTable.data, 1);
                obj.referenceTable.capacity = obj.referenceTable.used;
                obj.referenceTable.truePos = []; %Bool vec for call table rows
                obj.referenceTable.falseNeg = [];
                obj.referenceTable.checkedVs = false;
            end
        end

        function obj = matchCompareSpots(obj, thValA, thValB)
            %This needs to be rerun every time a new threshold for EITHER
            %tool is tried...
            if isempty(obj); return; end
            if isempty(obj.callTable); return; end
            if isempty(obj.callTableCompare); return; end
            if isempty(obj.callTableCompare.callTable); return; end
            
            tableA = obj.callTable;
            tableB = obj.callTableCompare.callTable;
            obj.callTableCompare.checkedThA = thValA;
            obj.callTableCompare.checkedThB = thValB;

            thOkayA = tableA{:, 'dropout_thresh'} >= thValA;
            thOkayB = tableB{:, 'dropout_thresh'} >= thValB;
            
            obj.callTableCompare.boolAOnly = false(sizeof(tableA, 1), 1);
            obj.callTableCompare.boolBOnly = false(sizeof(tableB, 1), 1);
            obj.callTableCompare.boolFullMatchA = false(sizeof(tableA, 1), 1);
            obj.callTableCompare.boolXYMatchA = false(sizeof(tableA, 1), 1);

            if nnz(thOkayA) < 1
                if nnz(thOkayB) > 0
                    %Everything is B only
                    obj.callTableCompare.boolBOnly(:,1) = true;
                end
                return;
            elseif nnz(thOkayB) < 1
                %Everything is A only
                obj.callTableCompare.boolAOnly(:,1) = true;
                return;
            end

            tableATh = tableA{thOkayA, :};
            tableBTh = tableB{thOkayB, :};
            mapA = find(thOkayA);
            mapB = find(thOkayB);

            %Hold one as reference against the other
            tableBR = NaN(size(tableBTh, 1), 3);
            tableBR(:,1) = tableBTh{:, 'isnap_x'};
            tableBR(:,2) = tableBTh{:, 'isnap_y'};
            tableBR(:,3) = tableBTh{:, 'isnap_z'};

            refAssign3 = RNACoords.matchSpots(tableATh, tableBR, obj.matchRad_xy, obj.matchRad_z, thValA);
            refAssign2 = RNACoords.matchSpots2D(tableATh, tableBR, obj.matchRad_xy, thValA);

            %I have no idea if this index wizardry will work. Remember to
            %test it.
            match3 = refAssign3(and(isfinite(refAssign3), (refAssign3 ~= 0)));
            match3 = mapA(match3);
            obj.callTableCompare.boolFullMatchA(match3) = true;

            match2 = refAssign2(and(isfinite(refAssign2), (refAssign2 ~= 0)));
            match2 = mapA(match2);
            obj.callTableCompare.boolXYMatchA(match2) = true;
            obj.callTableCompare.boolXYMatchA = and(obj.callTableCompare.boolXYMatchA, ~obj.callTableCompare.boolFullMatchA);

            obj.callTableCompare.boolAOnly = and(~obj.callTableCompare.boolXYMatchA, ~obj.callTableCompare.boolFullMatchA);

            noMatch = or(~isfinite(refAssign3), (refAssign3 == 0));
            noMatch = and(noMatch, or(~isfinite(refAssign2), (refAssign2 == 0)));
            noMatch = mapB(noMatch);
            obj.callTableCompare.boolBOnly(noMatch) = true;
        end

        function obj = matchReferenceSpots(obj)
            if isempty(obj); return; end
            if isempty(obj.callTable); return; end
            if isempty(obj.referenceTable); return; end
            if isempty(obj.referenceTable.data); return; end

            tableC = obj.callTable;
            tableR = obj.referenceTable.data(1:obj.referenceTable.used,:);

            ref_assign =  RNACoords.matchSpots(tableC, tableR, obj.matchRad_xy, obj.matchRad_z, obj.matchMinTh);
            obj.referenceTable.falseNeg = or(ref_assign < 1, ~isfinite(ref_assign));
            hitBool = ~obj.referenceTable.falseNeg;
            if nnz(hitBool) > 0
                obj.referenceTable.truePos = false(size(tableC, 1), 1);
                obj.referenceTable.truePos(ref_assign(hitBool), 1) = true;
            end

            obj.referenceTable.checkedVs = true;
        end

        function obj = reset(obj)
            %TODO
        end

    end

     %%
    methods(Static)

        function rowBool = filterCallsToValidRangeMtx(callMtx, validRegion, zMin, zMax)
            if nargin < 3; zMin = 0; end
            if nargin < 4; zMax = 0; end

            rowBool = [];
            if isempty(callMtx); return; end

            rowBool = true(size(callMtx, 1), 1);
            if ~isempty(validRegion)
                if validRegion.x_min > 1
                    rowBool = and(rowBool, callMtx(:,1) >= validRegion.x_min);
                end
                if validRegion.x_max > 0
                    rowBool = and(rowBool, callMtx(:,1) <= validRegion.x_max);
                end
                if validRegion.y_min > 1
                    rowBool = and(rowBool, callMtx(:,2) >= validRegion.y_min);
                end
                if validRegion.y_max > 0
                    rowBool = and(rowBool, callMtx(:,2) <= validRegion.y_max);
                end
                if validRegion.z_min > 1
                    rowBool = and(rowBool, callMtx(:,3) >= validRegion.z_min);
                end
                if validRegion.z_max > 0
                    rowBool = and(rowBool, callMtx(:,3) <= validRegion.z_max);
                end
            end

            if zMin > 1
                rowBool = and(rowBool, callMtx(:,3) >= zMin);
            end

            if zMax > 0
                rowBool = and(rowBool, callMtx(:,3) <= zMax);
            end
        end

        function rowBool = filterCallsToValidRange(callTable, validRegion, zMin, zMax, thVal)
            rowBool = [];
            if isempty(callTable); return; end

            rowBool = true(size(callTable, 1), 1);
            if ~isempty(validRegion)
                if validRegion.x_min > 1
                    rowBool = and(rowBool, callTable{:, 'isnap_x'} >= validRegion.x_min);
                end
                if validRegion.x_max > 0
                    rowBool = and(rowBool, callTable{:, 'isnap_x'} <= validRegion.x_max);
                end
                if validRegion.y_min > 1
                    rowBool = and(rowBool, callTable{:, 'isnap_y'} >= validRegion.y_min);
                end
                if validRegion.y_max > 0
                    rowBool = and(rowBool, callTable{:, 'isnap_y'} <= validRegion.y_max);
                end
                if validRegion.z_min > 1
                    rowBool = and(rowBool, callTable{:, 'isnap_z'} >= validRegion.z_min);
                end
                if validRegion.z_max > 0
                    rowBool = and(rowBool, callTable{:, 'isnap_z'} <= validRegion.z_max);
                end
            end

            if zMin > 1
                rowBool = and(rowBool, callTable{:, 'isnap_z'} >= zMin);
            end

            if zMax > 0
                rowBool = and(rowBool, callTable{:, 'isnap_z'} <= zMax);
            end

            if thVal > 0
                rowBool = and(rowBool, callTable{:, 'dropout_thresh'} >= thVal);
            end
        end

        function [figHandle, okay] = draw3DPlot(figHandle, mtx, base_color, idims)  
            okay = false;
            if isempty(figHandle); return; end
            if isempty(mtx); return; end
            if isempty(base_color)
                base_color = [1 0 0];
            end
            if isempty(idims)
                idims = struct();
                idims.x = round(max(mtx(:,1), [], 'all', 'omitnan'));
                idims.y = round(max(mtx(:,2), [], 'all', 'omitnan'));
                idims.z = round(max(mtx(:,3), [], 'all', 'omitnan'));
            end

            figure(figHandle);
            hold on;
            
            mtx = double(mtx);
            X = double(idims.x);
            Y = double(idims.y);
            Z = double(idims.z);
            max_rad = sqrt(X^2 + Y^2 + Z^2);
            
            count = size(mtx, 1);
            for i = 1:count
                %No wonder this thing is so slow...
                %Fix this eventually...
                x = mtx(i,1);
                y = mtx(i,2);
                z = mtx(i,3);
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
        function ax = get3DPlotAxes(idims)
            if isempty(idims)
                idims = struct();
                idims.x = 1;
                idims.y = 1;
                idims.z = 1;
            end

            ax = axes();
            ax.XLim = [1, idims.x];
            ax.YLim = [1, idims.y];
            ax.ZLim = [1, idims.z];
            ax.XAxisLocation = 'top';
            ax.YDir = 'reverse';
            ax.XLabel.String = 'X Axis';
            ax.YLabel.String = 'Y Axis';
            ax.ZLabel.String = 'Z Axis';
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

    end


end