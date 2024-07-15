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
        referenceTable;

        zMode; %Max proj(1) or single slice(2)?
        zNearRad; %In single slice mode, how many slices up and down are calls considered "near" and drawn in brighter colors?
        zShowRad; %In single slice mode, how many slices up and down to see calls from 
        workRegion; %Outside this region, pretend spots don't exist.

        visCommon; %Instance of VisCommon to get color tables from

        %%Internal scratch

    end

    %%
    methods

        function obj = initializeMe(obj)
            obj.zMode = 1;
            obj.zNearRad = 2;
            obj.zShowRad = 7;

            obj.workRegion = struct();
            obj.workRegion.x_min = 0;
            obj.workRegion.x_max = 0;
            obj.workRegion.y_min = 0;
            obj.workRegion.y_max = 0;
            obj.workRegion.z_min = 0;
            obj.workRegion.z_max = 0;

            obj.visCommon = VisCommon.empty();
        end

        function [figHandle, okay] = drawResultsBasic(obj, figHandle, thVal, z)
            %TODO
            okay = false;
        end

        function [figHandle, okay] = drawResultsTTCompare(obj, figHandle, thVal1, thVal2, z)
            %TODO
            okay = false;
        end

        function [figHandle, okay] = drawResultsTRCompare(obj, figHandle, thVal, z)
            %TODO
            okay = false;
        end

        function [figHandle, okay] = drawReferenceSet(obj, figHandle, z)
            %TODO
            okay = false;
        end

        function [figHandle, okay] = drawResultsBasic3D(obj, figHandle, thVal)
            %TODO
            okay = false;
        end

        function [figHandle, okay] = drawResultsTTCompare3D(obj, figHandle, thVal1, thVal2)
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

        function obj = refModeSnap(obj, zOnly)
            %TODO Runs spot snap (looking for local maxima)
        end

    end

     %%
    methods(Static)
    end


end