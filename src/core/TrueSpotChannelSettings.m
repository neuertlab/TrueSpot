%
%%
classdef TrueSpotChannelSettings

    %%
    properties (Constant)
        SAVE_VERSION = 1;
    end

    %%
    properties
        channelIndex = 0;
        metadata;
        spotCountSettings; %Includes quant
        thresholdSettings;
        options;
        previewState;
    end

    %%
    methods

        %% ========================== I/O ==========================
        function saveMe(obj, filePath)
            saveVer = obj.SAVE_VERSION;
            timeStamp = datetime;
            channelIndex = obj.channelIndex;
            metadata = obj.metadata;
            spotCountSettings = obj.spotCountSettings;
            thresholdSettings = obj.thresholdSettings;
            options = obj.options;
            previewState = obj.previewState;

            save(filePath, 'saveVer', 'timeStamp', 'channelIndex', 'metadata', ...
                'spotCountSettings', 'thresholdSettings', 'options', 'previewState');
        end

    end

    %%
    methods (Static)

        %% ========================== Inner Structs ==========================

        function metaInfo = genMetadataStruct()
            metaInfo = struct();
            metaInfo.targetName = [];
            metaInfo.probeName = [];
            metaInfo.targetMolType = [];
            metaInfo.pointSize = struct('x', 0, 'y', 0, 'z', 0);
        end

        function spotCountSettings = genCountSettingsStruct()
            spotCountSettings = struct();
            spotCountSettings.gaussRad = 7;
            spotCountSettings.zMin = 0;
            spotCountSettings.zMax = 0;
            spotCountSettings.useDPC = true;
            spotCountSettings.spotDetectThreads = 1;
            spotCountSettings.quantThreads = 1;
            spotCountSettings.controlImageChannel = 0;
            spotCountSettings.zAdj = NaN;
            spotCountSettings.fitRadXY = 4;
            spotCountSettings.fitRadZ = 2;
            spotCountSettings.quantNoRefilter = true;
            spotCountSettings.quantNoClouds = true;
            spotCountSettings.quantNoCells = false;
            spotCountSettings.quantCellZero = true;
            spotCountSettings.nucMaskLevel = 2; %2 is medium. 0 is 2D.
        end

        function thresholdSettings = genThresholdSettingsStruct()
            thresholdSettings = struct();
            thresholdSettings.thMin = 0; %Auto
            thresholdSettings.thMax = 0; %Auto
            thresholdSettings.preset = NaN;
            thresholdSettings.thParams = RNAThreshold.genEmptyThresholdParamStruct();
            thresholdSettings.manualTh = 0; %Picked by user in GUI
            thresholdSettings.manualThRangeMin = 0; %Picked by user in GUI
            thresholdSettings.manualThRangeMax = 0; %Picked by user in GUI
            thresholdSettings.autoBatchThMid = 0; %Derived from batch statistics
            thresholdSettings.autoBatchThLo = 0;
            thresholdSettings.autoBatchThHi = 0;
        end

        function optionsStruct = genOptionsStruct()
            optionsStruct = struct();
            optionsStruct.csvzero = false;
            optionsStruct.csvrange = false;
            optionsStruct.csvthonly = false;
            optionsStruct.csvfull = false;
            optionsStruct.runparamtxt = false;
            optionsStruct.overwrite = true;
        end

        %% ========================== I/O ==========================
        function chSettings = newChannelSettings()
            chSettings = TrueSpotChannelSettings;

            chSettings.metadata = TrueSpotChannelSettings.genMetadataStruct();
            chSettings.spotCountSettings = TrueSpotChannelSettings.genCountSettingsStruct();
            chSettings.thresholdSettings = TrueSpotChannelSettings.genThresholdSettingsStruct();
            chSettings.options = TrueSpotChannelSettings.genOptionsStruct();

            chSettings.previewState = struct();
            chSettings.previewState.workRegion = RNAUtils.genCoordRangeStruct(true);
            chSettings.previewState.viewRegion = RNAUtils.genCoordRangeStruct(true);
            chSettings.previewState.viewThreshold = 0;
            chSettings.previewState.channelVisible = false;
            chSettings.previewState.contrastMin = 0; %Value mapped to 0
            chSettings.previewState.contrastMax = 0; %Value mapped to 255
            chSettings.previewState.contrastMinFilt = 0; %Repeat, for filtered
            chSettings.previewState.contrastMaxFilt = 0;

        end

        function chSettings = loadChannelSettings(filePath)
            chSettings = [];
            if ~isfile(filePath); return; end

            chSettings = TrueSpotChannelSettings;
            load(filePath, 'channelIndex', 'metadata', ...
                'spotCountSettings', 'thresholdSettings', 'options', 'previewState');

            chSettings.channelIndex = channelIndex;
            chSettings.metadata = metadata;
            chSettings.spotCountSettings = spotCountSettings;
            chSettings.thresholdSettings = thresholdSettings;
            chSettings.options = options;
            chSettings.previewState = previewState;
        end


    end


end