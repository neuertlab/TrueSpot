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
        currentSlice = 1;

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

        %% ================== Getters ===============================

        %% ================== Setters ===============================

        function obj = setChannelFilterGaussRad(obj, channelNumber, gaussRad)
            if channelNumber < 1; return; end
            if channelNumber > obj.channelCount; return; end

            channelInfo = obj.myChannels{channelNumber};
            if channelInfo.appliedGaussRad == gaussRad; return; end %Don't bother with other processing.

            %Clear enough cache space for expected filtered image size...
            if obj.useFilt
                slice_size_est = 2 * obj.visCommon.idims.x * obj.visCommon.idims.y;
                stack_size_est = slice_size_est * (obj.visCommon.idims.z + 2); %For two max projections.
                obj = obj.clearMem(stack_size_est);
            end
            channelInfo = VisMultiCh.updateChannelFilter(channelInfo, obj.useFilt);
            
            obj.myChannels{channelNumber} = channelInfo;
        end

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
            clear currentChCache channelInfo

            %Clear any filtered if in raw mode or raw if in filtered mode
            for c = 1:obj.channelCount
                channelInfo = obj.myChannels{c};
                if ~isempty(channelInfo)
                    currentChCache = channelInfo.cacheUsed;
                    if currentChCache > 0
                        if obj.useFilt
                            channelInfo = VisMultiCh.clearLoadedImageData(channelInfo, false, true);
                        else
                            channelInfo = VisMultiCh.clearLoadedImageData(channelInfo, true, false);
                        end
                        obj.myChannels{c} = channelInfo;
                        currentMem = currentMem - (currentChCache - channelInfo.cacheUsed);
                        clear channelInfo;
                        if currentMem <= targetMem; return; end
                    end
                end
            end
            clear currentChCache channelInfo

            %Clear any stale max proj renders or out of range z slices
            z_min = obj.visibleRegion.z_min;
            z_max = obj.visibleRegion.z_max;
            if isnan(z_min); z_min = 1; end
            if isnan(z_max); z_max = obj.visCommon.idims.z; end
            for c = 1:obj.channelCount
                channelInfo = obj.myChannels{c};
                if ~isempty(channelInfo)
                    currentChCache = channelInfo.cacheUsed;
                    if(currentChCache > 0)
                        channelInfo = VisMultiCh.clearCachedZOutsideRange(channelInfo, z_min, z_max);
                        obj.myChannels{c} = channelInfo;
                        currentMem = currentMem - (currentChCache - channelInfo.cacheUsed);
                        clear channelInfo;
                        if currentMem <= targetMem; return; end
                    end
                end
            end

            %Clear any z slices or max proj not literally currently being
            %used for render
            if ~obj.maxProj
                for c = 1:obj.channelCount
                    channelInfo = obj.myChannels{c};
                    if ~isempty(channelInfo)
                        currentChCache = channelInfo.cacheUsed;
                        if(currentChCache > 0)
                            channelInfo = VisMultiCh.clearCachedZOutsideRange(channelInfo, obj.currentSlice, obj.currentSlice);
                            obj.myChannels{c} = channelInfo;
                            currentMem = currentMem - (currentChCache - channelInfo.cacheUsed);
                            clear channelInfo;
                            if currentMem <= targetMem; return; end
                        end
                    end
                end
            else
                for c = 1:obj.channelCount
                    channelInfo = obj.myChannels{c};
                    if ~isempty(channelInfo)
                        currentChCache = channelInfo.cacheUsed;
                        if(currentChCache > 0)
                            channelInfo = VisMultiCh.clearCachedZSlices(channelInfo);
                            obj.myChannels{c} = channelInfo;
                            currentMem = currentMem - (currentChCache - channelInfo.cacheUsed);
                            clear channelInfo;
                            if currentMem <= targetMem; return; end
                        end
                    end
                end
            end

        end

        function [obj, maxProjData] = getMaxProj(obj, channelNumber, z_min, z_max, filt)
            maxProjData = [];
            channelInfo = obj.myChannels{channelNumber};
            if (channelInfo.maxProj_zMin == z_min) & (channelInfo.maxProj_zMax == z_max)
                %May already be loaded
                if filt
                    if ~isempty(channelStruct.maxProjFilt)
                        maxProjData = channelStruct.maxProjFilt;
                        return;
                    end
                else
                    if ~isempty(channelStruct.maxProjRaw)
                        maxProjData = channelStruct.maxProjRaw;
                        return;
                    end
                end
            else
                %Cached max proj does not match current range.
                prevMem = 0;
                if ~isempty(channelInfo.maxProjFilt)
                    prevMem = prevMem + VisMultiCh.estimateMemSize(channelInfo.maxProjFilt);
                    channelInfo.maxProjFilt = [];
                end
                if ~isempty(channelInfo.maxProjRaw)
                    prevMem = prevMem + VisMultiCh.estimateMemSize(channelInfo.maxProjRaw);
                    channelInfo.maxProjRaw = [];
                end
                channelInfo.cacheUsed = channelInfo.cacheUsed - prevMem;
                channelInfo.maxProj_zMin = z_min;
                channelInfo.maxProj_zMax = z_max;
            end

            %If it reaches here, need to load.
            mpFieldName = ['mp_' num2str(z_min) '_' num2str(z_max)];
            if filt
                if isfile(channelInfo.filtCacheFilePath)
                    load(channelInfo.filtCacheFilePath, 'maxproj');
                    if isfield(maxproj, mpFieldName)
                        maxProjData = maxproj.(mpFieldName);
                        channelInfo.maxProjFilt = maxProjData;
                        channelInfo.cacheUsed = channelInfo.cacheUsed + VisMultiCh.estimateMemSize(maxProjData);
                    else
                        %Not previously generated.
                        load(channelInfo.filtCacheFilePath, 'imgFilt');
                        maxProjData = max(imgFilt(:,:,z_min:z_max), [], 3, 'omitnan');
                        clear imgFilt
                        channelInfo.maxProjFilt = maxProjData;
                        channelInfo.cacheUsed = channelInfo.cacheUsed + VisMultiCh.estimateMemSize(maxProjData);
                        maxproj.(mpFieldName) = maxProjData;
                        save(channelInfo.filtCacheFilePath, 'maxproj', '-append', '-v7.3');
                    end
                    clear maxproj
                else
                    slice_size_est = 2 * obj.visCommon.idims.x * obj.visCommon.idims.y;
                    stack_size_est = slice_size_est * (obj.visCommon.idims.z + 2); %For two max projections.
                    obj = obj.clearMem(stack_size_est);
                    channelInfo = VisMultiCh.updateChannelFilter(channelInfo, true);
                    maxProjData = channelInfo.maxProjFilt;
                end
            else
                if isfile(channelInfo.rawCacheFilePath)
                    load(channelInfo.rawCacheFilePath, 'maxproj');
                    if isfield(maxproj, mpFieldName)
                        maxProjData = maxproj.(mpFieldName);
                        channelInfo.maxProjRaw = maxProjData;
                        channelInfo.cacheUsed = channelInfo.cacheUsed + VisMultiCh.estimateMemSize(maxProjData);
                    else
                        %Not previously generated.
                        load(channelInfo.rawCacheFilePath, 'imageData');
                        maxProjData = max(imageData(:,:,z_min:z_max), [], 3, 'omitnan');
                        clear imageData
                        channelInfo.maxProjRaw = maxProjData;
                        channelInfo.cacheUsed = channelInfo.cacheUsed + VisMultiCh.estimateMemSize(maxProjData);
                        maxproj.(mpFieldName) = maxProjData;
                        save(channelInfo.rawCacheFilePath, 'maxproj', '-append', '-v7.3');
                    end
                    clear maxproj
                end %If file doesn't exist, just return empty matrix.
            end

            obj.myChannels{channelNumber} = channelInfo; %Cache has been modified.
        end

        function [obj, slice] = getZSlice(obj, channelNumber, z, filt)
            %slice = [];
            channelInfo = obj.myChannels{channelNumber};
            if filt
                %Gauss rad should be correct if updated correctly since
                %cache is purged when refiltering is done.
                sliceInfo = channelInfo.cacheFilt{z};
                if ~isempty(sliceInfo.imageData)
                    slice = sliceInfo.imageData;
                else
                    if ~isfile(channelInfo.filtCacheFilePath)
                        channelInfo = VisMultiCh.updateChannelFilter(channelInfo, true);
                        sliceInfo = channelInfo.cacheFilt{z};
                        slice = sliceInfo.imageData;
                    else
                        estSliceSize = 2 * obj.visCommon.idims.x * obj.visCommon.idims.y;
                        obj = obj.clearMem(estSliceSize);
                        load(channelInfo.filtCacheFilePath, 'imgFilt');
                        slice = imgFilt(:,:,z);
                        sliceInfo.imageData = slice;
                        clear imgFilt
                        channelInfo.cacheFilt{z} = sliceInfo;
                        channelInfo.cacheUsed = channelInfo.cacheUsed + VisMultiCh.estimateMemSize(slice);
                    end

                end
            else
                sliceInfo = channelInfo.cacheRaw{z};
                if ~isempty(sliceInfo.imageData)
                    slice = sliceInfo.imageData;
                else
                    %Load :(
                    estSliceSize = channelInfo.rawVoxelBytes * obj.visCommon.idims.x * obj.visCommon.idims.y;
                    obj = obj.clearMem(estSliceSize);
                    load(channelInfo.rawCacheFilePath, 'imageData');
                    slice = imageData(:,:,z);
                    sliceInfo.imageData = slice;
                    clear imageData
                    channelInfo.cacheRaw{z} = sliceInfo;
                    channelInfo.cacheUsed = channelInfo.cacheUsed + VisMultiCh.estimateMemSize(slice);
                end
            end

            obj.myChannels{channelNumber} = channelInfo;
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
            if nargin < 4; cacheDir = []; end
            if nargin < 5; loadListener = []; end
            %LoadListener must implement functions:
            %void onTifReadStart(void)
            %void onTifReadEnd(void)

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
            %Right now we can only load the whole TIF at once anyway soooooo
            if ~isempty(loadListener); loadListener.onTifReadStart(); end
            [stack, img_read] = tiffread2(tifpath);
            Z = img_read/visualizer.channelCount;
            Y = size(stack(1,1).data,1);
            X = size(stack(1,1).data,2);

            idims = struct('x', X, 'y', Y, 'z', Z);
            visualizer.visCommon = visualizer.visCommon.initializeMe(idims);

            chBuffer = zeros(Y,X,Z);
            if isa(stack(1,1).data, 'uint16')
                chBuffer = uint16(chBuffer);
            end
            for c = 1:channelCount
                zz = c;
                chStruct = visualizer.myChannels{c};
                for z = 1:Z
                    chBuffer(:,:,z) = stack(1,zz).data;
                    zz = zz + visualizer.channelCount;
                end
                chStruct = VisMultiCh.genDiskCacheRawImage(chStruct, chBuffer);
                visualizer.myChannels{c} = chStruct;
            end

            clear stack
            if ~isempty(loadListener); loadListener.onTifReadEnd(); end

        end


        %% ================== Structs ===============================

        function channelStruct = genEmptyChannelStruct()
            channelStruct = struct();
            channelStruct.viewName = [];
            channelStruct.isVisible = false;
            channelStruct.channelNumber = 0;
            channelStruct.gaussRad = 7;
            channelStruct.rawCacheFilePath = [];
            channelStruct.filtCacheFilePath = [];
            channelStruct.spotsrunFilePath = [];
            channelStruct.quantFilePath = [];

            %Cache
            channelStruct.cacheUsed = 0; %Approximate memory used by loaded image data.
            channelStruct.cacheRaw = []; %Cells of cachedSliceStructs - one per Z slice
            channelStruct.cacheFilt = [];
            channelStruct.appliedGaussRad = 0;
            channelStruct.rawVoxelBytes = 0;

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

        function value = estimateMemSize(data)
            value = 0;
            if isempty(data); return; end
            X = size(data, 2);
            Y = size(data, 1);
            Z = 1;
            if ndims(data) > 2
                Z = size(data, 3);
            end

            elementSize = 1;
            if isdouble(data)
                elementSize = 8; 
            elseif issingle(data)
                elementSize = 4;
            elseif isa(data,'uint16')
                elementSize = 2;
            elseif isa(data,'int16')
                elementSize = 2;
            elseif isa(data,'uint32')
                elementSize = 4;
            elseif isa(data,'int32')
                elementSize = 4;
            end

            value = X .* Y .* Z .* elementSize;
        end

        function channelStruct = clearCachedZSlices(channelStruct)
            if isempty(channelStruct); return; end

            memCleared = 0;
            Z = size(channelStruct.cacheRaw, 2);
            for z = 1:Z
                sliceDat = channelStruct.cacheFilt{z};
                if ~isempty(sliceDat)
                    if ~isempty(sliceDat.imageData)
                        clearAmt = VisMultiCh.estimateMemSize(sliceDat.imageData);
                        sliceDat.imageData = [];
                        memCleared = memCleared + clearAmt;
                    end
                end
                channelStruct.cacheFilt{z} = sliceDat;

                sliceDat = channelStruct.cacheRaw{z};
                if ~isempty(sliceDat)
                    if ~isempty(sliceDat.imageData)
                        clearAmt = VisMultiCh.estimateMemSize(sliceDat.imageData);
                        sliceDat.imageData = [];
                        memCleared = memCleared + clearAmt;
                    end
                end
                channelStruct.cacheRaw{z} = sliceDat;
            end
            
            channelStruct.cacheUsed = channelStruct.cacheUsed - memCleared;
        end

        function channelStruct = clearCachedZOutsideRange(channelStruct, z_min_valid, z_max_valid)
            if isempty(channelStruct); return; end
            
            %Max projections
            memCleared = 0;
            if ~isnan(channelStruct.maxProj_zMin) & (channelStruct.maxProj_zMin > 0)
                if channelStruct.maxProj_zMin ~= z_min_valid
                    clearAmt = VisMultiCh.estimateMemSize(channelStruct.maxProjFilt);
                    channelStruct.maxProjFilt = [];
                    clearAmt = VisMultiCh.estimateMemSize(channelStruct.maxProjRaw) + clearAmt;
                    channelStruct.maxProjRaw = [];
                    channelStruct.maxProj_zMin = 0;
                    channelStruct.maxProj_zMax = 0;
                    memCleared = memCleared + clearAmt;
                end
            end

            if ~isnan(channelStruct.maxProj_zMax) & (channelStruct.maxProj_zMax > 0)
                if channelStruct.maxProj_zMax ~= z_max_valid
                    clearAmt = VisMultiCh.estimateMemSize(channelStruct.maxProjFilt);
                    channelStruct.maxProjFilt = [];
                    clearAmt = VisMultiCh.estimateMemSize(channelStruct.maxProjRaw) + clearAmt;
                    channelStruct.maxProjRaw = [];
                    channelStruct.maxProj_zMin = 0;
                    channelStruct.maxProj_zMax = 0;
                    memCleared = memCleared + clearAmt;
                end
            end

            %Z Slices
            for z = 1:(z_min_valid - 1)
                sliceDat = channelStruct.cacheFilt{z};
                if ~isempty(sliceDat)
                    if ~isempty(sliceDat.imageData)
                        clearAmt = VisMultiCh.estimateMemSize(sliceDat.imageData);
                        sliceDat.imageData = [];
                        memCleared = memCleared + clearAmt;
                    end
                end
                channelStruct.cacheFilt{z} = sliceDat;

                sliceDat = channelStruct.cacheRaw{z};
                if ~isempty(sliceDat)
                    if ~isempty(sliceDat.imageData)
                        clearAmt = VisMultiCh.estimateMemSize(sliceDat.imageData);
                        sliceDat.imageData = [];
                        memCleared = memCleared + clearAmt;
                    end
                end
                channelStruct.cacheRaw{z} = sliceDat;
            end

            Z = size(channelStruct.cacheRaw, 2);
            for z = (z_max_valid + 1):Z
                sliceDat = channelStruct.cacheFilt{z};
                if ~isempty(sliceDat)
                    if ~isempty(sliceDat.imageData)
                        clearAmt = VisMultiCh.estimateMemSize(sliceDat.imageData);
                        sliceDat.imageData = [];
                        memCleared = memCleared + clearAmt;
                    end
                end
                channelStruct.cacheFilt{z} = sliceDat;

                sliceDat = channelStruct.cacheRaw{z};
                if ~isempty(sliceDat)
                    if ~isempty(sliceDat.imageData)
                        clearAmt = VisMultiCh.estimateMemSize(sliceDat.imageData);
                        sliceDat.imageData = [];
                        memCleared = memCleared + clearAmt;
                    end
                end
                channelStruct.cacheRaw{z} = sliceDat;
            end
            
            channelStruct.cacheUsed = channelStruct.cacheUsed - memCleared;
        end
      
        function channelStruct = clearLoadedImageData(channelStruct, clearFilt, clearRaw)
            if nargin < 2; clearFilt = true; end
            if nargin < 3; clearRaw = true; end
            if isempty(channelStruct); return; end
            
            if clearFilt
                memCleared = 0;
                if ~isempty(channelStruct.maxProjFilt)
                    clearAmt = VisMultiCh.estimateMemSize(channelStruct.maxProjFilt);
                    channelStruct.maxProjFilt = [];
                    memCleared = memCleared + clearAmt;
                end
                Z = size(channelStruct.cacheFilt, 2);
                for z = 1:Z
                    sliceDat = channelStruct.cacheFilt{z};
                    if ~isempty(sliceDat)
                        if ~isempty(sliceDat.imageData)
                            clearAmt = VisMultiCh.estimateMemSize(sliceDat.imageData);
                            sliceDat.imageData = [];
                            memCleared = memCleared + clearAmt;
                        end
                    end
                    channelStruct.cacheFilt{z} = sliceDat;
                end
                clear sliceDat;
                channelStruct.cacheUsed = channelStruct.cacheUsed - memCleared;
            end

            if clearRaw
                memCleared = 0;
                if ~isempty(channelStruct.maxProjRaw)
                    clearAmt = VisMultiCh.estimateMemSize(channelStruct.maxProjRaw);
                    channelStruct.maxProjRaw = [];
                    memCleared = memCleared + clearAmt;
                end
                Z = size(channelStruct.cacheRaw, 2);
                for z = 1:Z
                    sliceDat = channelStruct.cacheRaw{z};
                    if ~isempty(sliceDat)
                        if ~isempty(sliceDat.imageData)
                            clearAmt = VisMultiCh.estimateMemSize(sliceDat.imageData);
                            sliceDat.imageData = [];
                            memCleared = memCleared + clearAmt;
                        end
                    end
                    channelStruct.cacheRaw{z} = sliceDat;
                end
                clear sliceDat;
                channelStruct.cacheUsed = channelStruct.cacheUsed - memCleared;
            end

            if isempty(channelStruct.maxProjFilt) & isempty(channelStruct.maxProjRaw)
                channelStruct.maxProj_zMin = 0;
                channelStruct.maxProj_zMax = 0;
            end
        end

        function channelStruct = genDiskCacheRawImage(channelStruct, imgCh)
            if isempty(channelStruct); return; end
            if isempty(imgCh); return; end

            Z = size(imgCh, 3);
            channelStruct.cacheRaw = cell(1, Z);
            channelStruct.cacheFilt = cell(1, Z);
            for z = 1:Z
                sliceData = VisMultiCh.genEmptyCachedSliceStruct();
                sliceData.z = z;
                channelStruct.cacheRaw{z} = sliceData;

                sliceData = VisMultiCh.genEmptyCachedSliceStruct();
                sliceData.z = z;
                channelStruct.cacheFilt{z} = sliceData;
            end

            %Raw
            imageData = imgCh; %Changing the name for save
            clear imgCh;

            if isdouble(imageData)
                channelStruct.rawVoxelBytes = 8;
            elseif isa(imageData, 'uint16')
                channelStruct.rawVoxelBytes = 2;
            elseif isa(imageData, 'int16')
                channelStruct.rawVoxelBytes = 2;
            elseif isa(imageData, 'uint32')
                channelStruct.rawVoxelBytes = 4;
            elseif isa(imageData, 'int32')
                channelStruct.rawVoxelBytes = 4;
            elseif issingle(imageData)
                channelStruct.rawVoxelBytes = 4;
            elseif isa(imageData, 'uint8')
                channelStruct.rawVoxelBytes = 1;
            elseif isa(imageData, 'int8')
                channelStruct.rawVoxelBytes = 1;
            end

            %Raw max proj
            maxproj = struct();
            mpName = ['mp_1_' num2str(Z)];
            maxproj.(mpName) = max(imageData, [], 3, 'omitnan');
            save(channelStruct.rawCacheFilePath, 'imageData', 'maxproj', '-v7.3');
            clear maxproj imageData
        end

        function channelStruct = updateChannelFilter(channelStruct, holdInMem)
            if isempty(channelStruct); return; end
            if channelStruct.gaussRad == channelStruct.appliedGaussRad
                return;
            end

            [cacheDir, cacheName, ~] = fileparts(channelStruct.filtCacheFilePath);
            dpcStem = [cacheDir filesep cacheName];
            channelStruct = VisMultiCh.clearLoadedImageData(channelStruct, true, false);
            load(channelStruct.rawCacheFilePath, 'imageData');
            [imgFilt] = RNA_Threshold_SpotDetector.run_spot_detection_pre(...
                imageData, dpcStem, true, channelStruct.gaussRad, false);
            clear imageData

            Z = size(imgFilt, 3);

            maxproj = struct();
            mpName = ['mp_1_' num2str(Z)];
            maxproj.(mpName) = max(imgFilt, [], 3, 'omitnan');

            if (channelStruct.maxProj_zMin > 0)
                mp_min = channelStruct.maxProj_zMin;
                mp_max = channelStruct.maxProj_zMax;
                if (mp_min) > 1 | (mp_max < Z)
                    mpName = ['mp_' num2str(mp_min) '_' num2str(mp_max)];
                    maxproj.(mpName) = max(imgFilt(:,:,mp_min:mp_max), [], 3, 'omitnan');
                end

                if holdInMem
                    channelStruct.maxProjFilt = maxproj.(mpName); %Cleared before, so imageData should be empty
                    nowSize = VisMultiCh.estimateMemSize(channelStruct.maxProjFilt);
                    channelStruct.cacheUsed = channelStruct.cacheUsed + nowSize;
                    clear nowSize
                end
            end

            %Cache new slices
            if holdInMem
                for z = 1:Z
                    sliceInfo = channelStruct.cacheFilt{z}; %Cleared before, so imageData should be empty
                    sliceInfo.imageData = imgFilt(:,:,z);
                    nowSize = VisMultiCh.estimateMemSize(sliceInfo.imageData);
                    channelStruct.cacheFilt{z} = sliceInfo;
                    channelStruct.cacheUsed = channelStruct.cacheUsed + nowSize;
                end
            end

            save(channelStruct.filtCacheFilePath, 'imgFilt', 'maxproj', '-v7.3');
            clear maxproj imgFilt

            channelStruct.appliedGaussRad = channelStruct.gaussRad;
        end


    end

end