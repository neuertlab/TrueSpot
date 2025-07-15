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

        inputTable; %Held in memory only, not saved.
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

        function obj = loadInputTable(obj, filePath)
            %A csv or tsv table

            delim = '\t';
            if endsWith(filePath, '.csv'); delim = ','; end

            [varNames, ~] = TrueSpotProjSettings.getInputTableFields();
            temp_table = readtable(filePath,'Delimiter',delim,'ReadVariableNames',true);

            %Fill in missing variables with empty columns
            isVar = RNAUtils.isTableVariable(temp_table, varNames);
            varCount = size(varNames, 2);
            for i = 1:varCount
                if ~isVar(i)
                    temp_table{:, varNames{i}} = "";
                end
            end

            obj.inputTable = temp_table;

        end

        function obj = scanInputDir(obj, dirPath)
            [varNames, varTypes] = TrueSpotProjSettings.getInputTableFields();

            dirList = dir(dirPath);
            dirTable = struct2table(dirList);
            dirTable = convertvars(dirTable, 'name', 'string');

            fileBool = ~dirTable{:,'isdir'};
            fNames = dirTable{:, 'name'};
            tifBool = endsWith(fNames, '.tif');
            passBool = and(fileBool, tifBool);
            iCount = nnz(passBool);
            clear fileBool tifBool fNames dirList

            table_size = [iCount size(varNames,2)];
            obj.inputTable = table('Size', table_size, 'VariableTypes', varTypes, 'VariableNames', varNames);
            obj.inputTable{:, 'ImageFilePath'} = dirTable{passBool, 'name'};
            obj.inputTable{:, 'ImageName'} = replace(dirTable{passBool, 'name'}, '.tif', '');
            obj.inputTable{:, 'ControlFilePath'} = "";
            obj.inputTable{:, 'NucDataFilePath'} = "";
            obj.inputTable{:, 'LightDataFilePath'} = "";
        end

        function exportInputTable(obj, filePath)
            if ~isempty(obj.inputTable); return; end

            delim = '\t';
            if endsWith(filePath, '.csv'); delim = ','; end
            writetable(obj.inputTable, filePath, 'Delimiter',delim)
        end

        %% ========================== Run ==========================

        function obj = initializeForRun(obj)
            %Loads and fills in input table (ie. if common control)
            if isfolder(obj.paths.inputPath)
                obj = obj.scanInputDir(obj.paths.inputPath);
                if ~isempty(obj.paths.controlPath)
                    obj.inputTable{:, 'ControlFilePath'} = obj.paths.controlPath;
                end
            else
                %Tif image or table
                if endsWith(obj.paths.inputPath, '.tif')
                    [varNames, varTypes] = TrueSpotProjSettings.getInputTableFields();
                    table_size = [1 size(varNames,2)];
                    obj.inputTable = table('Size', table_size, 'VariableTypes', varTypes, 'VariableNames', varNames);
                    if ~isempty(obj.paths.controlPath)
                        obj.inputTable{:, 'ControlFilePath'} = obj.paths.controlPath;
                    end
                else
                    %Try to read as table
                    obj = obj.loadInputTable(obj.paths.inputPath);
                end
            end

            %Prepare output. Export the inputtable as well for user.
            if ~isfolder(obj.paths.outputPath)
                mkdir(obj.paths.outputPath);
            end
            obj.exportInputTable([obj.paths.outputPath filesep 'tsInput.csv']);
        end

        function obj = runCellSeg(obj, listener)
            %TODO
            if ~isempty(obj.inputTable)
                if ~isempty(listener)
                    listener.sendMessage('[TrueSpotProjSettings.runCellSeg] ERROR: Batch has not been initialized.');
                end
                return;
            end

            if ~isempty(listener)
                listener.sendMessage('[TrueSpotProjSettings.runCellSeg] Starting cell segmentation...');
            end

            %TODO

        end

        function obj = runSpotDetect(obj, channels, listener)
            %TODO
        end

        function obj = runQuant(obj, channels, listener)
            %TODO
            %Remember to dump a count dump xml
        end

        function obj = runQC(obj, listener)
            %TODO
        end

        function obj = runCountDump(obj, listener)
            %TODO
        end

    end

    %%
    methods (Static)

        %% ========================== Utility ==========================

        function [varNames, varTypes] = getInputTableFields()
            varNames = {'ImageName' 'ImageFilePath' 'ControlFilePath' 'NucDataFilePath' 'LightDataFilePath'};
            varTypes = {'string' 'string' 'string' 'string' 'string'};
        end
        
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