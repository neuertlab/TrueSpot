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

        %% ========================== I/O ==========================
        function batchSettings = newBatchSettings()
            batchSettings = TrueSpotProjSettings;

            %Tables are used to match samples, controls, and light/nuc
            %channels stored in different files
            batchSettings.paths = struct();
            batchSettings.paths.inputPath = []; %Can be image, directory, or table
            batchSettings.paths.controlPath = []; %Can be image, directory, or table
            batchSettings.paths.outputPath = []; %Directory

            batchSettings.metadata = struct();
            batchSettings.metadata.species = [];
            batchSettings.metadata.cellType = [];
            batchSettings.metadata.voxelSize = struct('x', 0, 'y', 0, 'z', 0);
            batchSettings.metadata.batchName = [];

            batchSettings.channelInfo = struct();
            batchSettings.channelInfo.channelCount = 0;
            batchSettings.channelInfo.nucMarkerChannel = 0;
            batchSettings.channelInfo.lightChannel = 0;
            batchSettings.channelInfo.controlChannelCount = 0;
            batchSettings.channelInfo.controlNucMarkerChannel = 0;
            batchSettings.channelInfo.controlLightChannel = 0;

            batchSettings.cellsegSettings = struct();
            batchSettings.cellsegSettings.presetName = [];
            batchSettings.cellsegSettings.nucZMin = 0;
            batchSettings.cellsegSettings.nucZMax = 0;
            batchSettings.cellsegSettings.lightZMin = 0;
            batchSettings.cellsegSettings.lightZMax = 0;
            batchSettings.cellsegSettings.cszmin = 0;
            batchSettings.cellsegSettings.cszmax = 0;
            batchSettings.cellsegSettings.nszmin = 0;
            batchSettings.cellsegSettings.nszmax = 0;
            batchSettings.cellsegSettings.xtrim = NaN;
            batchSettings.cellsegSettings.ytrim = NaN;
            batchSettings.cellsegSettings.nzrange = NaN;
            batchSettings.cellsegSettings.nthsmpl = NaN;
            batchSettings.cellsegSettings.ncutoff = NaN;
            batchSettings.cellsegSettings.ndxy = NaN;
            batchSettings.cellsegSettings.outputCellMaskPNG = false;
            batchSettings.cellsegSettings.outputCellMaskTIF = false;
            batchSettings.cellsegSettings.outputNucMaskPNG = false;
            batchSettings.cellsegSettings.outputNucMaskTIF = false;
            batchSettings.cellsegSettings.overwrite = true;
            batchSettings.cellsegSettings.dumpSettings = false;

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