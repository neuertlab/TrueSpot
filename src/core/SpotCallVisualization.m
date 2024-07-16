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

            obj.refRemoveQueue = struct();
            obj.refRemoveQueue.queue = []; %1D table of row indices in reftable
            obj.refRemoveQueue.capacity = 0;
            obj.refRemoveQueue.used = 0;

            obj.refAddQueue = struct();
            obj.refAddQueue.queue = []; %Matrix w/ cols x,y,z
            obj.refAddQueue.capacity = 0;
            obj.refAddQueue.used = 0;

            obj.referenceTable = struct();
            obj.referenceTable.data = []; %Matrix w/ cols x,y,z
            obj.referenceTable.capacity = 0;
            obj.referenceTable.used = 0;
            obj.referenceTable.truePos = []; %Bool vec for call table rows
            obj.referenceTable.falseNeg = []; %Bool vec for ref table rows
            obj.referenceTable.checkedVs = false;

            obj.callTableCompare = struct();
            obj.callTableCompare.callTable = [];
            obj.callTableCompare.boolAOnly = [];
            obj.callTableCompare.boolBOnly = [];
            obj.callTableCompare.boolFullMatchA = [];
            obj.callTableCompare.boolXYMatchA = [];
            obj.callTableCompare.checkedThA = 0;
            obj.callTableCompare.checkedThB = 0;

            obj.visCommon = VisCommon.empty();
        end

        function [figHandle, okay] = drawResultsBasic(obj, figHandle, thVal, z)
            okay = false;
            if isempty(obj); return; end
            if isempty(figHandle); return; end

            if obj.zMode == 1
                %MaxProj
                rowBool = SpotCallVisualization.filterCallsToValidRange(obj.callTable, obj.workRegion, 0, 0, thVal);
                if nnz(rowBool) > 0
                    xx = obj.callTable{rowBool, 'isnap_x'};
                    yy = obj.callTable{rowBool, 'isnap_y'};

                    figure(figHandle);
                    hold on;
                    plot(xx, yy,'LineStyle','none','Marker','o','MarkerEdgeColor','r','markersize',10);
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

                    figure(figHandle);
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
            if isempty(figHandle); return; end
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

            figure(figHandle);
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
            if isempty(figHandle); return; end
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

            figure(figHandle);
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
            if isempty(figHandle); return; end
            if isempty(obj.referenceTable); return; end
            if isempty(obj.referenceTable.data); return; end

            tableR = obj.referenceTable.data(1:obj.referenceTable.used,:);

            if obj.zMode == 1
                rowBoolRef = SpotCallVisualization.filterCallsToValidRangeMtx(tableR, obj.workRegion, 0, 0);

                if nnz(rowBoolRef) > 0
                    figure(figHandle);
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

                    figure(figHandle);
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

                    groupBool = and(rowBool, tableR(:, 3) >= (z - obj.zNearRad));
                    groupBool = and(groupBool, tableR(:, 3) < z);
                    if nnz(groupBool) > 0
                        plot(tableR(groupBool, 1), tableR(groupBool, 2),...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',brightYellow,'markersize',10);
                    end

                    groupBool = and(rowBool, tableR(:, 3) <= (z + obj.zNearRad));
                    groupBool = and(groupBool, tableR(:, 3) > z);
                    if nnz(groupBool) > 0
                        plot(tableR(groupBool, 1), tableR(groupBool, 2),...
                            'LineStyle','none','Marker','o','MarkerEdgeColor',brightGreen,'markersize',10);
                    end

                    groupBool = and(rowBool, tableR(:, 3) == z);
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

        function [obj, figHandle] = onRefModeClick(obj, figHandle, x, y)
            %TODO Draws new circle and updates queue
        end

        function obj = onRefModeAccept(obj)
            %TODO updates ref table to reflect click queue
        end

        function obj = refModeSnap(obj, zOnly, filteredImg)
            %TODO Runs spot snap (looking for local maxima)
        end

        function obj = matchCompareSpots(obj, thValA, thValB)
            %TODO
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

            %Hold one as reference against the other
            tableBR = NaN(size(tableBTh, 1), 3);
            tableBR(:,1) = tableBTh{:, 'isnap_x'};
            tableBR(:,2) = tableBTh{:, 'isnap_y'};
            tableBR(:,3) = tableBTh{:, 'isnap_z'};

            refAssign3 = RNACoords.matchSpots(tableATh, tableBR, obj.matchRad_xy, obj.matchRad_z, thValA);
            refAssign2 = RNACoords.matchSpots2D(tableATh, tableBR, obj.matchRad_xy, thValA);

            %TODO
            %These indices also need to be mapped back to original tables
            %at some point


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

    end


end