%
%%
classdef VisMultiCh

    %Ability to cache image channels to disk instead of holding in mem.

    %%
    properties

        %Disk & Caching
        imagePath = [];
        channelCount = 0;
        cacheDir = []; %Defaults to directory of input image
        cacheNameStem = [];
        preloadedMemTarget = 0x100000000;

        %View Region
        visibleRegion = RNAUtils.genCoordRangeStruct(true);

        %View Options (Toggles)
        useFilt = false; %Toggle between filtered and raw image channel
        maxProj = true; %Toggle between single slice view and max projection
        globalContrast = true; %(Slice view only) Toggle between contrast scale from full image or just that slice

        imageLayerOn = true;
        cellLayerOn = false;
        nucLayerOn = false;
        spotCircleLayerOn = false;
        quantFitLayerOn = false;
        cloudLayerOn = false;

        %Analysis Data Viewing
        cellsegMaskPath = [];
        myChannels = []; %Array of channel structs

        visCommon = VisCommon;

    end

    %%
    methods

        %% ================== Init ===============================

        %% ================== Destroy ===============================

        %% ================== Caching ===============================

        function value = getCachedMemEstimate(obj)
            value = 0;
            for c = 1:obj.channelCount
                channelInfo = obj.myChannels{c};
                if ~isempty(channelInfo)
                    value = value + channelInfo.cacheUsed;
                end
            end
        end

        function obj = clearMem(obj, amt)
            %TODO
            targetMem = obj.preloadedMemTarget - amt;
            currentMem = obj.getCachedMemEstimate();

            if currentMem <= targetMem; return; end

            %Clear any currently invisible channels
            for c = 1:obj.channelCount
                channelInfo = obj.myChannels{c};
                if ~isempty(channelInfo)
                    if ~channelInfo.isVisible
                        currentChCache = channelInfo.cacheUsed;
                        if currentChCache > 0
                            channelInfo = VisMultiCh.clearLoadedImageData(channelInfo);
                            obj.myChannels{c} = channelInfo;
                            clear channelInfo;
                            currentMem = currentMem - currentChCache;
                            if currentMem <= targetMem; return; end
                        end
                    end
                end
            end

            %Clear any filtered if in raw mode or raw if in filtered mode

            %Clear any stale max proj renders or out of range z slices

            

        end

        function [obj, maxProj] = getMaxProj(obj, channelNumber, z_min, z_max, filt)
            %TODO
            channelInfo = obj.myChannels{channelNumber};
        end

        function [obj, slice] = getZSlice(obj, channelNumber, z, filt)
            %TODO
        end

        function [obj, stack] = getZSlices(obj, channelNumber, z_min, z_max, filt)
            %TODO
        end

        %% ================== Rendering ===============================
        function [obj, channelBW2D] = renderChannelImage(obj, channelNumber)
            %TODO
        end

    end

    %%
    methods (Static)
        %% ================== Init ===============================
        function visualizer = newVisMultiCh(tifPath, channelCount, cacheDir, loadListener)
            if nargin < 3; cacheDir = []; end
            if nargin < 4; loadListener = []; end
            %TODO

            visualizer = VisMultiCh;
            if isempty(cacheDir)
                [tifDir, tifName, ~] = fileparts(tifPath);
                visualizer.cacheDir = tifDir;
                visualizer.cacheNameStem = tifName;
                clear tifDir tifName
            else
                [~, tifName, ~] = fileparts(tifPath);
                visualizer.cacheDir = cacheDir;
                visualizer.cacheNameStem = tifName;
                clear tifName
            end

            visualizer.channelCount = channelCount;
            visualizer.myChannels = cell(1, channelCount);
            for c = 1:channelCount
                chStruct = VisMultiCh.genEmptyChannelStruct();
                chStruct.channelNumber = c;
                chStruct.viewName = ['CH' num2str(c)];
                chStruct.rawCacheFilePath = [visualizer.cacheDir filesep visualizer.cacheNameStem '_ch' num2str(c) '_viscacheraw.mat'];
                chStruct.filtCacheFilePath = [visualizer.cacheDir filesep visualizer.cacheNameStem '_ch' num2str(c) '_viscachefilt.mat'];
                visualizer.myChannels{c} = chStruct;
            end

            %If image is small enough, load and pre-cache raw data
            %TODO

        end


        %% ================== Structs ===============================

        function channelStruct = genEmptyChannelStruct()
            channelStruct = struct();
            channelStruct.viewName = [];
            channelStruct.isVisible = false;
            channelStruct.channelNumber = 0;
            channelStruct.rawCacheFilePath = [];
            channelStruct.filtCacheFilePath = [];
            channelStruct.spotsrunFilePath = [];
            channelStruct.quantFilePath = [];

            %Cache
            channelStruct.cacheUsed = 0; %Approximate memory used by loaded image data.
            channelStruct.cacheRaw = []; %Cells of cachedSliceStructs - one per Z slice
            channelStruct.cacheFilt = [];

            channelStruct.maxProj_zMin = 0;
            channelStruct.maxProj_zMax = 0;
            channelStruct.maxProjRaw = [];
            channelStruct.maxProjFilt = [];
        end

        function cachedSliceStruct = genEmptyCachedSliceStruct()
            cachedSliceStruct = struct();
            cachedSliceStruct.z = 0;
            cachedSliceStruct.imageData = [];
        end

        %% ================== Caching ===============================

        function channelStruct = clearLoadedImageData(channelStruct)
            if isempty(channelStruct); return; end
            channelStruct.maxProj_zMin = 0;
            channelStruct.maxProj_zMax = 0;
            channelStruct.maxProjRaw = [];
            channelStruct.maxProjFilt = [];

            Z = size(channelStruct.cacheRaw, 2);
            for z = 1:Z
                sliceDat = channelStruct.cacheRaw{z};
                if ~isempty(sliceDat)
                    if ~isempty(sliceDat.imageData)
                        sliceDat.imageData = [];
                    end
                end
                channelStruct.cacheRaw{z} = sliceDat;
            end

            Z = size(channelStruct.cacheFilt, 2);
            for z = 1:Z
                sliceDat = channelStruct.cacheFilt{z};
                if ~isempty(sliceDat)
                    if ~isempty(sliceDat.imageData)
                        sliceDat.imageData = [];
                    end
                end
                channelStruct.cacheFilt{z} = sliceDat;
            end
            clear sliceDat;

            channelStruct.cacheUsed = 0;
        end

    end

end