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

        %% ========================== Run ==========================

        function obj = runSpotDetect(obj, parentProject, imgName, imgPath, ctrlPath, listener)
            if nargin < 6; listener = []; end

            funcTag = 'TrueSpotChannelSettings.runSpotDetect';
            scratchList = cell(1, 256);
            argPos = 1;

            scratchList{argPos} = '-input'; argPos = argPos + 1;
            scratchList{argPos} = imgPath; argPos = argPos + 1;
            scratchList{argPos} = '-imgname'; argPos = argPos + 1;
            scratchList{argPos} = imgName; argPos = argPos + 1;

            scratchList{argPos} = '-chtotal'; argPos = argPos + 1;
            scratchList{argPos} = num2str(parentProject.channelInfo.channelCount); argPos = argPos + 1;
            scratchList{argPos} = '-chsamp'; argPos = argPos + 1;
            scratchList{argPos} = num2str(obj.channelIndex); argPos = argPos + 1;
            scratchList{argPos} = '-chtrans'; argPos = argPos + 1;
            scratchList{argPos} = num2str(parentProject.channelInfo.lightChannel); argPos = argPos + 1;

            stackDir = [parentProject.paths.outputPath filesep imgName];
            chDir = [stackDir filesep 'CH' num2str(obj.channelIndex)];

            %scratchList{argPos} = '-outdir'; argPos = argPos + 1;
            %scratchList{argPos} = chDir; argPos = argPos + 1;
            scratchList{argPos} = '-outstem'; argPos = argPos + 1;
            scratchList{argPos} = [chDir filesep 'TS_' imgName '_CH' num2str(obj.channelIndex)]; argPos = argPos + 1;

            cellMaskPath = [stackDir filesep 'CellSeg_' imgName '.mat'];
%             if ~isempty(extCellMaskPath)
%                 cellMaskPath = extCellMaskPath;
%             end
            if isfile(cellMaskPath)
                scratchList{argPos} = '-cellseg'; argPos = argPos + 1;
                scratchList{argPos} = cellMaskPath; argPos = argPos + 1;
            end
            if ~isempty(ctrlPath)
                scratchList{argPos} = '-ctrlimg'; argPos = argPos + 1;
                scratchList{argPos} = ctrlPath; argPos = argPos + 1;
                scratchList{argPos} = '-chctrtotal'; argPos = argPos + 1;
                scratchList{argPos} = num2str(parentProject.channelInfo.controlChannelCount); argPos = argPos + 1;
                scratchList{argPos} = '-chctrsamp'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.controlImageChannel); argPos = argPos + 1;
            end

            %Tuning params
            if (obj.spotCountSettings.gaussRad > 0)
                scratchList{argPos} = '-gaussrad'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.gaussRad); argPos = argPos + 1;
            end
            if (obj.spotCountSettings.zMin > 0)
                scratchList{argPos} = '-zmin'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.zMin); argPos = argPos + 1;
            end
            if (obj.spotCountSettings.zMax > 0)
                scratchList{argPos} = '-zmax'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.zMax); argPos = argPos + 1;
            end
            if ~obj.spotCountSettings.useDPC
                scratchList{argPos} = '-nodpc'; argPos = argPos + 1;
            end
            if (obj.thresholdSettings.thMin > 0)
                scratchList{argPos} = '-thmin'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.thresholdSettings.thMin); argPos = argPos + 1;
            else
                scratchList{argPos} = '-autominth'; argPos = argPos + 1;
            end
            if (obj.thresholdSettings.thMax > 0)
                scratchList{argPos} = '-thmax'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.thresholdSettings.thMax); argPos = argPos + 1;
            else
                scratchList{argPos} = '-automaxth'; argPos = argPos + 1;
            end

            if ~isnan(obj.thresholdSettings.preset)
                if obj.thresholdSettings.preset >= 0
                    scratchList{argPos} = '-sensitivity'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.thresholdSettings.preset); argPos = argPos + 1;
                else
                    scratchList{argPos} = '-precision'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.thresholdSettings.preset * -1); argPos = argPos + 1;
                end
            else
                %Gather all the threshold tuning params :)
                scratchList{argPos} = '-mfmin'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.thresholdSettings.thParams.mad_factor_min); argPos = argPos + 1;
                scratchList{argPos} = '-mfmax'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.thresholdSettings.thParams.mad_factor_max); argPos = argPos + 1;

                winmin = min(obj.thresholdSettings.thParams.window_sizes, [], 'all');
                winmax = max(obj.thresholdSettings.thParams.window_sizes, [], 'all');
                winincr = min(diff(obj.thresholdSettings.thParams.window_sizes), [], 'all');
                scratchList{argPos} = '-wszmin'; argPos = argPos + 1;
                scratchList{argPos} = num2str(winmin); argPos = argPos + 1;
                scratchList{argPos} = '-wszmax'; argPos = argPos + 1;
                scratchList{argPos} = num2str(winmax); argPos = argPos + 1;
                scratchList{argPos} = '-wszincr'; argPos = argPos + 1;
                scratchList{argPos} = num2str(winincr); argPos = argPos + 1;

                scratchList{argPos} = '-thmw'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.thresholdSettings.thParams.madth_weight); argPos = argPos + 1;
                scratchList{argPos} = '-thfw'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.thresholdSettings.thParams.fit_weight); argPos = argPos + 1;
                scratchList{argPos} = '-thiw'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.thresholdSettings.thParams.fit_ri_weight); argPos = argPos + 1;

                if obj.thresholdSettings.thParams.test_data
                    scratchList{argPos} = '-usespc'; argPos = argPos + 1;
                end
                if obj.thresholdSettings.thParams.test_diff
                    scratchList{argPos} = '-usedfc'; argPos = argPos + 1;
                end
                if (obj.thresholdSettings.thParams.std_factor ~= 0.0)
                    scratchList{argPos} = '-stdfac'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.thresholdSettings.thParams.std_factor); argPos = argPos + 1;
                end

                scratchList{argPos} = '-logproj'; argPos = argPos + 1;
                if (obj.thresholdSettings.thParams.log_proj_mode == 0)
                    scratchList{argPos} = 'None'; argPos = argPos + 1;
                elseif (obj.thresholdSettings.thParams.log_proj_mode == 1)
                    scratchList{argPos} = 'All'; argPos = argPos + 1;
                elseif (obj.thresholdSettings.thParams.log_proj_mode == 2)
                    scratchList{argPos} = 'FitOnly'; argPos = argPos + 1;
                end
            end

            if (obj.spotCountSettings.spotDetectThreads > 0)
                scratchList{argPos} = '-threads'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.spotDetectThreads); argPos = argPos + 1;
            end

            %Output options
            if obj.options.csvzero | obj.options.csvrange | obj.options.csvthonly | obj.options.csvfull
                scratchList{argPos} = '-csvout'; argPos = argPos + 1;
                scratchList{argPos} = [chDir filesep 'LocalMaxCalls_' imgName '.csv']; argPos = argPos + 1;
            end

            if obj.options.runparamtxt
                scratchList{argPos} = '-runparamout'; argPos = argPos + 1;
                scratchList{argPos} = [chDir filesep 'SpotDetectRunParams_' imgName '.txt']; argPos = argPos + 1;
            end
            if obj.options.overwrite
                scratchList{argPos} = '-ovrw'; argPos = argPos + 1;
            end
            if obj.options.csvzero
                scratchList{argPos} = '-csvzero'; argPos = argPos + 1;
            end
            if obj.options.csvrange
                scratchList{argPos} = '-csvrange'; argPos = argPos + 1;
            end
            if obj.options.csvthonly
                scratchList{argPos} = '-csvthonly'; argPos = argPos + 1;
            end
            if obj.options.csvfull
                scratchList{argPos} = '-csvfull'; argPos = argPos + 1;
            end

            %Metadata
            voxSize = parentProject.metadata.voxelSize;
            scratchList{argPos} = '-voxelsize'; argPos = argPos + 1;
            scratchList{argPos} = ['(' num2str(voxSize.x) ',' num2str(voxSize.y) ',' num2str(voxSize.z) ')']; argPos = argPos + 1;
            spotSize = obj.metadata.pointSize;
            scratchList{argPos} = '-expspotsize'; argPos = argPos + 1;
            scratchList{argPos} = ['(' num2str(spotSize.x) ',' num2str(spotSize.y) ',' num2str(spotSize.z) ')']; argPos = argPos + 1;

            if ~isempty(obj.metadata.targetName)
                scratchList{argPos} = '-target'; argPos = argPos + 1;
                scratchList{argPos} = obj.metadata.targetName; argPos = argPos + 1;
            end
            if ~isempty(obj.metadata.probeName)
                scratchList{argPos} = '-probetype'; argPos = argPos + 1;
                scratchList{argPos} = obj.metadata.probeName; argPos = argPos + 1;
            end
            if ~isempty(obj.metadata.targetMolType)
                scratchList{argPos} = '-targettype'; argPos = argPos + 1;
                scratchList{argPos} = obj.metadata.targetMolType; argPos = argPos + 1;
            end
            if ~isempty(parentProject.metadata.species)
                scratchList{argPos} = '-species'; argPos = argPos + 1;
                scratchList{argPos} = parentProject.metadata.species; argPos = argPos + 1;
            end
            if ~isempty(parentProject.metadata.cellType)
                scratchList{argPos} = '-celltype'; argPos = argPos + 1;
                scratchList{argPos} = parentProject.metadata.cellType; argPos = argPos + 1;
            end
            scratchList{argPos} = '-verbose'; argPos = argPos + 1;

            myArgs = scratchList(1:(argPos-1));
            if ~isempty(listener)
                %Log generated command
                cmdString = strjoin(myArgs, ' ');
                listener.logMessage(['[' funcTag ']' 'Generated argument list: ' cmdString]);
            end

            Main_RNASpots(myArgs{:});
        end

        function obj = runQuant(obj, parentProject, imgName, explicitImagePath, explicitCellSegPath, listener)
            if nargin < 6; listener = []; end

            funcTag = 'TrueSpotChannelSettings.runQuant';
            scratchList = cell(1, 256);
            argPos = 1;

            stackDir = [parentProject.paths.outputPath filesep imgName];
            chDir = [stackDir filesep 'CH' num2str(obj.channelIndex)];
            
            outstem = [chDir filesep 'TS_' imgName '_CH' num2str(obj.channelIndex)];
            runPath = [outstem '_rnaspotsrun.mat'];

            if isfile(runPath)
                scratchList{argPos} = '-runinfo'; argPos = argPos + 1;
                scratchList{argPos} = runPath; argPos = argPos + 1;
            else
                %Probably no data if no spotsrun, but try anyway.
                scratchList{argPos} = '-coordtable'; argPos = argPos + 1;
                scratchList{argPos} = [chDir filesep 'TS_' imgName '_callTable.mat']; argPos = argPos + 1;
                scratchList{argPos} = '-chcount'; argPos = argPos + 1;
                scratchList{argPos} = num2str(parentProject.channelInfo.channelCount); argPos = argPos + 1;
                scratchList{argPos} = '-ch'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.channelIndex); argPos = argPos + 1;
            end

            if ~isempty(explicitImagePath)
                scratchList{argPos} = '-tif'; argPos = argPos + 1;
                scratchList{argPos} = explicitImagePath; argPos = argPos + 1;
            end
            if ~isempty(explicitCellSegPath)
                scratchList{argPos} = '-cellsegpath'; argPos = argPos + 1;
                scratchList{argPos} = explicitCellSegPath; argPos = argPos + 1;
            end

            if (obj.thresholdSettings.manualThRangeMin > 0) | (obj.thresholdSettings.manualThRangeMax > 0)
                if (obj.thresholdSettings.manualThRangeMin > 0)
                    scratchList{argPos} = '-thmin'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.thresholdSettings.manualThRangeMin); argPos = argPos + 1;
                end
                if (obj.thresholdSettings.manualThRangeMax > 0)
                    scratchList{argPos} = '-thmax'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.thresholdSettings.manualThRangeMax); argPos = argPos + 1;
                end
            else
                if obj.thresholdSettings.manualTh > 0
                    scratchList{argPos} = '-thmin'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.thresholdSettings.manualTh - 100); argPos = argPos + 1;
                    scratchList{argPos} = '-thmax'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.thresholdSettings.manualTh + 100); argPos = argPos + 1;
                end
            end

            %Tuning
            if(obj.spotCountSettings.quantThreads > 1)
                scratchList{argPos} = '-workers'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.quantThreads); argPos = argPos + 1;
            end
            if ~isnan(obj.spotCountSettings.zAdj)
                scratchList{argPos} = '-zadj'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.zAdj); argPos = argPos + 1;
            end
            if(obj.spotCountSettings.fitRadXY > 0)
                scratchList{argPos} = '-radxy'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.fitRadXY); argPos = argPos + 1;
            end
            if(obj.spotCountSettings.fitRadZ > 0)
                scratchList{argPos} = '-radz'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.spotCountSettings.fitRadZ); argPos = argPos + 1;
            end
            if obj.spotCountSettings.quantNoRefilter
                scratchList{argPos} = '-norefilter'; argPos = argPos + 1;
            end
            if obj.spotCountSettings.quantNoClouds
                scratchList{argPos} = '-noclouds'; argPos = argPos + 1;
            end
            if obj.spotCountSettings.quantCellZero
                scratchList{argPos} = '-cellzero'; argPos = argPos + 1;
            end

            myArgs = scratchList(1:(argPos-1));
            if ~isempty(listener)
                %Log generated command
                cmdString = strjoin(myArgs, ' ');
                listener.logMessage(['[' funcTag ']' 'Generated argument list: ' cmdString]);
            end

            Main_RNAQuant(myArgs{:});
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
            spotCountSettings.fitRadXY = 0;
            spotCountSettings.fitRadZ = 0;
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