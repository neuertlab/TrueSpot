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

        %% ========================== I/O ==========================
        function chSettings = newChannelSettings()
            chSettings = TrueSpotChannelSettings;

            chSettings.metadata = struct();
            chSettings.metadata.targetName = [];
            chSettings.metadata.probeName = [];
            chSettings.metadata.targetMolType = [];
            chSettings.metadata.pointSize = struct('x', 0, 'y', 0, 'z', 0);

            chSettings.spotCountSettings = struct();
            chSettings.spotCountSettings.gaussRad = 7;
            chSettings.spotCountSettings.zMin = 0;
            chSettings.spotCountSettings.zMax = 0;
            chSettings.spotCountSettings.useDPC = true;
            chSettings.spotCountSettings.spotDetectThreads = 1;
            chSettings.spotCountSettings.quantThreads = 1;
            chSettings.spotCountSettings.controlImageChannel = 0;
            chSettings.spotCountSettings.zAdj = NaN;
            chSettings.spotCountSettings.fitRadXY = 4;
            chSettings.spotCountSettings.fitRadZ = 2;
            chSettings.spotCountSettings.quantNoRefilter = true;
            chSettings.spotCountSettings.quantNoClouds = true;
            chSettings.spotCountSettings.quantNoCells = false;
            chSettings.spotCountSettings.quantCellZero = true;
            chSettings.spotCountSettings.nucMaskLevel = 2; %2 is medium. 0 is 2D.

            chSettings.thresholdSettings = struct();
            chSettings.thresholdSettings.thMin = 0; %Auto
            chSettings.thresholdSettings.thMax = 0; %Auto
            chSettings.thresholdSettings.preset = 0;
            chSettings.thresholdSettings.thParams = RNAThreshold.genEmptyThresholdParamStruct();
            chSettings.thresholdSettings.manualTh = 0; %Picked by user in GUI
            chSettings.thresholdSettings.manualThRangeMin = 0; %Picked by user in GUI
            chSettings.thresholdSettings.manualThRangeMax = 0; %Picked by user in GUI
            chSettings.thresholdSettings.autoBatchThMid = 0; %Derived from batch statistics
            chSettings.thresholdSettings.autoBatchThLo = 0;
            chSettings.thresholdSettings.autoBatchThHi = 0;

            chSettings.options = struct();
            chSettings.options.csvzero = false;
            chSettings.options.csvrange = false;
            chSettings.options.csvthonly = false;
            chSettings.options.csvfull = false;
            chSettings.options.runparamtxt = false;
            chSettings.options.overwrite = true;

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