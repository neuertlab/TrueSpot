%
%%
classdef BaseImageRenderer

    %%
    properties
        stackData8; %Raw stack (x,y,z,c) as uint8

        channelRenderInfo; %ChannelRenderStruct vector
        channelFiltCache; %FilterCacheStruct vector

        clutStore;
    end

    %%
    methods

        %---------------------------- Init ---------------------------

        function obj = initDefaultCLUTStore(obj)
            defoClutPath = './util/clut.mat';
            if isfile(defoClutPath)
                load(defoClutPath, 'clutSet');
                obj.clutStore = clutSet;
            else
                %Could not be found. Default to just "grey"
                gg = double(0:255);
                gg = gg ./ 255.0;
                greyClut = zeros(256, 3);
                greyClut(:, 1) = gg;
                greyClut(:, 2) = gg;
                greyClut(:, 3) = gg;
                obj.clutStore.grey = greyClut;
            end
        end
        
        %---------------------------- Load Source Images ---------------------------

        function obj = loadTif(obj, tifPath, channelCount, seriesIndex, subImageIndex)
            if nargin < 3; channelCount = 1; end
            if nargin < 4; seriesIndex = 1; end
            if nargin < 5; subImageIndex = 1; end

            loadTifOpStruct = LoadTifPlusStructGen();
            loadTifOpStruct.chCount = channelCount;
            loadResults = LoadTifPlusSingle(tifPath, loadTifOpStruct, seriesIndex, subImageIndex);

            channelCount = size(loadResults.imageData, 2);
            X = loadResults.imageDims.x;
            Y = loadResults.imageDims.y;
            Z = loadResults.imageDims.z;
            obj.stackData8 = uint8(zeros(Y, X, Z, channelCount));

            obj.channelRenderInfo(channelCount) = genChannelRenderStruct();
            obj.channelFiltCache(channelCount) = genFilterCacheStuct();

            %Rescale and import based on type
            for c = 1:channelCount
                chDat = loadResults.imageData{c};
                dead_pix_info = RNAUtils.detectDeadPixels(chDat, false);

                if isa(chDat, 'uint8')
                    %As-is
                    obj.stackData8(:,:,:,c) = chDat;
                elseif isa(chDat, 'uint16')
                    %Rescale based on 0-65535
                    chDat = double(chDat) ./ double(0xffff);
                    chDat = uint8(round(chDat .* 255.0));
                    obj.stackData8(:,:,:,c) = chDat;
                else
                    %Rescale based on minimum and maximum
                    chDat = BaseImageRenderer.linearRescaleToU8FullRange(chDat);
                end
                
                %Structs
                crStruct = genChannelRenderStruct();
                crStruct.alpha = 255;

                %Auto detect raw contrast
                [crStruct.contrastMin, crStruct.contrastMax] = ...
                    BaseImageRenderer.autoPickRawContrastRange(chDat);

                obj.channelRenderInfo(c) = crStruct;

                cfStruct = genFilterCacheStuct();
                cfStruct.deadPixInfo = dead_pix_info;
                obj.channelFiltCache(c) = cfStruct;
            end
        end

        %---------------------------- Easy Get/Set Channel States ---------------------------

        %---------------------------- Filtered Image Cache ---------------------------

        function obj = clearCachedFilteredChannels(obj)
            if ~isempty(obj.channelFiltCache)
                cCount = size(obj.channelFiltCache, 2);
                for c = 1:cCount
                    cfStruct = obj.channelFiltCache(c);
                    cfStruct.filtData8 = [];
                    obj.channelFiltCache(c) = cfStruct;
                end
            end
        end

        function obj = setChannelFilterGaussRad(obj, channel, gaussRad)
            if isempty(obj.channelFiltCache); return; end
            if channel < 1; return; end
            cCount = size(obj.channelFiltCache, 2);
            if channel > cCount; return; end

            cfStruct = obj.channelFiltCache(channel);
            oldRad = cfStruct.gaussRad;
            cfStruct.gaussRad = gaussRad;
            if oldRad ~= gaussRad
                cfStruct.filtData8 = [];
            end
            obj.channelFiltCache(channel) = cfStruct;
        end

        function obj = setChannelFilterDPC(obj, channel, dpcBool)
            if isempty(obj.channelFiltCache); return; end
            if channel < 1; return; end
            cCount = size(obj.channelFiltCache, 2);
            if channel > cCount; return; end

            cfStruct = obj.channelFiltCache(channel);
            oldBool = cfStruct.useDPC;
            cfStruct.useDPC = dpcBool;
            if oldBool ~= dpcBool
                cfStruct.filtData8 = [];
            end
            obj.channelFiltCache(channel) = cfStruct;
        end

        function obj = applyFilterToChannelAndCache(obj, channel, autoUpdateContrastBool)
            if isempty(obj.channelFiltCache); return; end
            if channel < 1; return; end
            cCount = size(obj.channelFiltCache, 2);
            if channel > cCount; return; end

            cfStruct = obj.channelFiltCache(channel);
            rawDat = obj.stackData8(:,:,:,channel);

            crStruct = obj.channelRenderInfo(channel);
            if crStruct.useDPC & ~isempty(cfStruct.deadPixInfo)
                rawDat = RNAUtils.cleanDeadPixels(rawDat, cfStruct.deadPixInfo, false);
            end

            filtDat = RNA_Threshold_Common.applyGaussianFilter(rawDat, cfStruct.gaussRad, cfStruct.gaussAmt);
            filtDat = RNA_Threshold_Common.applyEdgeDetectFilter(filtDat);
            filtDat = RNA_Threshold_Common.blackoutBorders(filtDat, cfStruct.gaussRad + 1, 0);
            cfStruct.filtData8 = BaseImageRenderer.linearRescaleToU8FullRange(filtDat);

            if autoUpdateContrastBool
                [crStruct.fContrastMin, crStruct.fContrastMax] = ...
                    BaseImageRenderer.autoPickFiltContrastRange(cfStruct.filtData8);
                obj.channelRenderInfo(channel) = crStruct;
            end
            
            obj.channelFiltCache(channel) = cfStruct;
        end

        %---------------------------- Render ---------------------------

        function bw8Image = renderSingleChannel(obj, channel, zSlice)
            if nargin < 3; zSlice = 0; end

            bw8Image = [];

            if isempty(obj.channelRenderInfo); return; end
            if channel < 1; return; end
            cCount = size(obj.channelRenderInfo, 2);
            if channel > cCount; return; end

            crStruct = obj.channelRenderInfo(channel);

            %zSlice 0 means max projection
            if crStruct.filtMode
                cfStruct = obj.channelFiltCache(channel);
                if isempty(cfStruct.filtData8)
                    obj = obj.applyFilterToChannelAndCache(channel, false);
                    cfStruct = obj.channelFiltCache(channel);
                end
                chDat = cfStruct.filtData8;
            else
                chDat = obj.stackData8(:,:,:,channel);
            end

            if zSlice > 0
                sliceData = chDat(:,:,zSlice);
            else
                sliceData = max(chDat, [], 3, 'omitnan');
            end

            %Linear rescale
            if crStruct.filtMode
                rsMax = (crStruct.fContrastMax - crStruct.fContrastMin);
                sliceData = sliceData - crStruct.fContrastMin;
            else
                rsMax = (crStruct.contrastMax - crStruct.contrastMin);
                sliceData = sliceData - crStruct.contrastMin;
            end

            sliceData = sliceData ./ rsMax;
            sliceData(sliceData < 0.0) = 0.0;
            sliceData(sliceData > 1.0) = 1.0;

            sliceData = sliceData .* 255.0;
            bw8Image = uint8(round(sliceData));
        end

        function rgb8Image = renderRawComposite(obj, zSlice)
            if nargin < 2; zSlice = 0; end

            cCount = size(obj.channelRenderInfo, 2);
            % alphaTotal = 0.0;
            % 
            % for c = 1:cCount
            %     crStruct = obj.channelRenderInfo(c);
            %     alphaTotal = alphaTotal + crStruct.alpha;
            % end

            for c = 1:cCount
                crStruct = obj.channelRenderInfo(c);
                if crStruct.alpha <= 0
                    continue;
                end

                chDat = obj.stackData8(:,:,:,c);
                if zSlice > 0
                    sliceData = chDat(:,:,zSlice);
                else
                    sliceData = max(chDat, [], 3, 'omitnan');
                end

                %Rescale to contrast (0.0 - 1.0)
                rsMax = (crStruct.contrastMax - crStruct.contrastMin);
                sliceData = sliceData - crStruct.contrastMin;
                sliceData = sliceData ./ rsMax;
                sliceData(sliceData < 0.0) = 0.0;
                sliceData(sliceData > 1.0) = 1.0;

                %Weight of each pixel applied to alpha for each channel
                %TODO
            end
        end

    end

    %%
    methods(Static)

        %---------------------------- Structs ---------------------------

        function fcStruct = genFilterCacheStuct()
            fcStruct.gaussRad = 5;
            fcStruct.gaussAmt = 2;
            fcStruct.useDPC = true;
            fcStruct.deadPixInfo = []; %Applies to raw channel, so don't need to recalculate
            fcStruct.filtData8 = [];
        end

        function crStruct = genChannelRenderStruct()
            crStruct.alpha = 0;
            crStruct.clutName = 'grey';
            crStruct.contrastMin = 0;
            crStruct.contrastMax = 255;

            %For LoG view - only applicable in single channel view
            crStruct.filtMode = false; 
            crStruct.fContrastMin = 0;
            crStruct.fContrastMin = 255;

            %crStruct.filtCache = BaseImageRenderer.genFilterCacheStuct();
        end

    %---------------------------- Image Rendering ---------------------------

        %%
        function [cmin, cmax] = autoPickRawContrastRange(channelStackData8)
            maxproj = double(max(channelStackData8, [], 3, 'omitnan'));
            cmin = min(maxproj(:));
            cmax = median(maxproj(:)) + round(10 * std(maxproj(:)));
        end

        %%
        function [cmin, cmax] = autoPickFiltContrastRange(channelStackData8)
            maxproj = double(max(channelStackData8, [], 3, 'omitnan'));
            cmin = median(maxproj(:)) - round(0 * std(maxproj(:)));
            cmax = median(maxproj(:)) + round(10 * std(maxproj(:)));
        end

        %---------------------------- uint8 ---------------------------

        function image8 = linearRescaleToU8FullRange(imageIn)
            %Basically, minimum (whatever it is) is set to 0, and maximum
            %to 255
            valMin = min(imageIn, [], 'all', 'omitnan');
            imageIn = imageIn - valMin; %Normalize min to 0
            valMax = max(imageIn, [], 'all', 'omitnan');
            imageIn = imageIn ./ valMax; %Normalize max to 1

            image8 = uint8(round(imageIn .* 255.0));
        end

    end
end