%
%%
classdef TrueSpotProjSettings

    %%
    properties (Constant)
        SAVE_VERSION = 1;
    end

    %%
    properties
        paths;
        channelInfo;
        cellsegSettings;
        metadata;
        previewState;
    end

    methods

        %% ========================== I/O ==========================
        function saveMe(obj, filePath)
            saveVer = obj.SAVE_VERSION;
            timeStamp = datetime;

            paths = obj.paths;
            channelInfo = obj.channelInfo;
            cellsegSettings = obj.cellsegSettings;
            metadata = obj.metadata;
            previewState = obj.previewState;
            save(filePath, 'saveVer', 'timeStamp', 'paths', 'channelInfo', ...
                'cellsegSettings', 'metadata', 'previewState');
        end

    end

    %%
    methods (Static)

        %% ========================== Inner Structs ==========================

        function pathsInfo = genPathsStruct()
            pathsInfo = struct();
            pathsInfo.inputPath = []; %Can be image, directory, or table
            pathsInfo.controlPath = []; %Can be image, directory, or table
            pathsInfo.outputPath = []; %Directory

            %External cellseg
            %No since these are per file, they'd have to be tabled...
%             pathsInfo.extCellMask = []; %2D or 3D image file
%             pathsInfo.extNucMask = []; %2D or 3D image file
        end

        function metaInfo = genMetadataStruct()
            metaInfo = struct();
            metaInfo.species = [];
            metaInfo.cellType = [];
            metaInfo.voxelSize = struct('x', 0, 'y', 0, 'z', 0);
            metaInfo.batchName = [];
        end

        function channelInfo = genChannelInfoStruct()
            channelInfo = struct();
            channelInfo.channelCount = 0;
            channelInfo.nucMarkerChannel = 0;
            channelInfo.lightChannel = 0;
            channelInfo.controlChannelCount = 0;
            channelInfo.controlNucMarkerChannel = 0;
            channelInfo.controlLightChannel = 0;
        end

        function cellsegSettings = genCellsegSettingsStruct()
            cellsegSettings = struct();
            cellsegSettings.presetName = [];
            cellsegSettings.nucZMin = 0;
            cellsegSettings.nucZMax = 0;
            cellsegSettings.lightZMin = 0;
            cellsegSettings.lightZMax = 0;
            cellsegSettings.cszmin = 0;
            cellsegSettings.cszmax = 0;
            cellsegSettings.nszmin = 0;
            cellsegSettings.nszmax = 0;
            cellsegSettings.xtrim = NaN;
            cellsegSettings.ytrim = NaN;
            cellsegSettings.nzrange = NaN;
            cellsegSettings.nthsmpl = NaN;
            cellsegSettings.ncutoff = NaN;
            cellsegSettings.ndxy = NaN;
            cellsegSettings.outputCellMaskPNG = false;
            cellsegSettings.outputCellMaskTIF = false;
            cellsegSettings.outputNucMaskPNG = false;
            cellsegSettings.outputNucMaskTIF = false;
            cellsegSettings.overwrite = true;
            cellsegSettings.dumpSettings = false;
        end

        %% ========================== I/O ==========================
        function batchSettings = newBatchSettings()
            batchSettings = TrueSpotProjSettings;

            %Tables are used to match samples, controls, and light/nuc
            %channels stored in different files
            batchSettings.paths = TrueSpotProjSettings.genPathsStruct();

            batchSettings.metadata = TrueSpotProjSettings.genMetadataStruct();
            batchSettings.channelInfo = TrueSpotProjSettings.genChannelInfoStruct();
            batchSettings.cellsegSettings = TrueSpotProjSettings.genCellsegSettingsStruct();

            batchSettings.previewState = struct();
            batchSettings.previewState.currentImageFile = [];
            batchSettings.previewState.viewRegion = RNAUtils.genCoordRangeStruct(true);
            batchSettings.previewState.lightChannelVisible = false;
            batchSettings.previewState.dapiChannelVisible = false;
            batchSettings.previewState.toggleMaxProjection = true;
            batchSettings.previewState.toggleGlobalContrast = true;
            batchSettings.previewState.toggleFilterOn = false;
            batchSettings.previewState.toggleCellMaskOutline = false; %If true, cell masks are transparent overlay. If false, outline.
            batchSettings.previewState.toggleLayerCellMask = false;
            batchSettings.previewState.toggleLayerNucMask = false;
            batchSettings.previewState.toggleLayerCellNum = false;
            batchSettings.previewState.toggleLayerSpotCalls = true;
            batchSettings.previewState.toggleLayerQuantFits = false;
            batchSettings.previewState.toggleLayerClouds = false;
        end

        function batchSettings = loadBatchSettings(filePath)
            batchSettings = [];
            if ~isfile(filePath); return; end

            load(filePath, 'paths', 'channelInfo', ...
                'cellsegSettings', 'metadata', 'previewState');

            batchSettings.paths = paths;
            batchSettings.channelInfo = channelInfo;
            batchSettings.cellsegSettings = cellsegSettings;
            batchSettings.metadata = metadata;
            batchSettings.previewState = previewState;
        end

        

    end


end