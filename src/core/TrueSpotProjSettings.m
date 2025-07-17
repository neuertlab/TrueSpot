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
            %obj.inputTable{:, 'ImageFilePath'} = dirTable{passBool, 'name'};
            obj.inputTable{:, 'ImageName'} = replace(dirTable{passBool, 'name'}, '.tif', '');
            obj.inputTable{:, 'ImageFilePath'} = arrayfun(@(tifFileName) strjoin([dirPath filesep tifFileName], ''), dirTable{passBool, 'name'});
            obj.inputTable{:, 'ControlFilePath'} = "";
            obj.inputTable{:, 'NucDataFilePath'} = "";
            obj.inputTable{:, 'ExtCellSegPath'} = "";
            obj.inputTable{:, 'ExtNucSegPath'} = "";
        end

        function exportInputTable(obj, filePath)
            if isempty(obj.inputTable); return; end

            delim = '\t';
            if endsWith(filePath, '.csv'); delim = ','; end
            writetable(obj.inputTable, filePath, 'Delimiter',delim)
        end

        %% ========================== Utils ==========================

        function value = getInputTableValue(obj, row, fieldName)
            value = [];
            rawVal = obj.inputTable{row, fieldName};
            if isempty(rawVal); return; end
            if ismissing(rawVal); return; end
            if (rawVal == ""); return; end
            value = char(rawVal);
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
            if nargin < 2; listener = []; end

            funcTag = 'TrueSpotProjSettings.runCellSeg';

            if isempty(obj.inputTable)
                GUIUtils.sendMessage(listener, funcTag, 'ERROR: Batch has not been initialized.');
                return;
            end

            GUIUtils.sendMessage(listener, funcTag, 'Starting cell segmentation...');

            TUNING_PARAMS = {'cszmin' 'cszmax' 'nszmin' 'nszmax' 'xtrim' 'ytrim' ...
                'nzrange' 'nthsmpl' 'ncutoff' 'ndxy'};
            tuningParamCount = size(TUNING_PARAMS, 2);

            stackCount = size(obj.inputTable, 1);
            scratchList = cell(1, 128);
            for i = 1:stackCount
                GUIUtils.sendMessage(listener, funcTag, sprintf('Working on image %d of %d ...', i, stackCount));

                iname = obj.getInputTableValue(i, 'ImageName');
                altLightPath = obj.getInputTableValue(i, 'LightDataFilePath');
                
                argPos = 1;
                scratchList{argPos} = '-input'; argPos = argPos + 1;
                if ~isempty(altLightPath)
                    scratchList{argPos} = altLightPath; argPos = argPos + 1;
                    scratchList{argPos} = '-chtotal'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.channelInfo.lightStackChannelCount); argPos = argPos + 1;
                else
                    scratchList{argPos} = obj.getInputTableValue(i, 'ImageFilePath'); argPos = argPos + 1;
                    scratchList{argPos} = '-chtotal'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.channelInfo.channelCount); argPos = argPos + 1;
                end

                altNucPath = obj.getInputTableValue(i, 'NucDataFilePath');
                if ~isempty(altNucPath)
                    scratchList{argPos} = '-innuc'; argPos = argPos + 1;
                    scratchList{argPos} = altNucPath; argPos = argPos + 1;
                    scratchList{argPos} = '-chtotnuc'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.channelInfo.nucStackChannelCount); argPos = argPos + 1;
                end

                scratchList{argPos} = '-chlight'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.channelInfo.lightChannel); argPos = argPos + 1;
                scratchList{argPos} = '-chnuc'; argPos = argPos + 1;
                scratchList{argPos} = num2str(obj.channelInfo.nucMarkerChannel); argPos = argPos + 1;

                outdirPath = [obj.paths.outputPath filesep iname];
                scratchList{argPos} = '-outpath'; argPos = argPos + 1;
                scratchList{argPos} = outdirPath; argPos = argPos + 1;

                scratchList{argPos} = '-imgname'; argPos = argPos + 1;
                scratchList{argPos} = iname; argPos = argPos + 1;

                if (obj.cellsegSettings.lightZMin > 0)
                    scratchList{argPos} = '-lightzmin'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.cellsegSettings.lightZMin); argPos = argPos + 1;
                end
                if (obj.cellsegSettings.lightZMax > 0)
                    scratchList{argPos} = '-lightzmax'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.cellsegSettings.lightZMax); argPos = argPos + 1;
                end
                if (obj.cellsegSettings.nucZMin > 0)
                    scratchList{argPos} = '-nuczmin'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.cellsegSettings.nucZMin); argPos = argPos + 1;
                end
                if (obj.cellsegSettings.nucZMax > 0)
                    scratchList{argPos} = '-nuczmin'; argPos = argPos + 1;
                    scratchList{argPos} = num2str(obj.cellsegSettings.nucZMax); argPos = argPos + 1;
                end
                if (obj.cellsegSettings.outputCellMaskPNG)
                    scratchList{argPos} = '-ocellmask'; argPos = argPos + 1;
                    scratchList{argPos} = [outdirPath filesep 'CellMask_' iname '.png']; argPos = argPos + 1;
                end
                if (obj.cellsegSettings.outputCellMaskTIF)
                    scratchList{argPos} = '-ocellmask'; argPos = argPos + 1;
                    scratchList{argPos} = [outdirPath filesep 'CellMask_' iname '.tif']; argPos = argPos + 1;
                end
                if (obj.cellsegSettings.outputNucMaskPNG)
                    scratchList{argPos} = '-onucmask'; argPos = argPos + 1;
                    scratchList{argPos} = [outdirPath filesep 'NucMask_' iname '.png']; argPos = argPos + 1;
                end
                if (obj.cellsegSettings.outputNucMaskTIF)
                    scratchList{argPos} = '-oonucmask'; argPos = argPos + 1;
                    scratchList{argPos} = [outdirPath filesep 'NucMask_' iname '.tif']; argPos = argPos + 1;
                end
                if (obj.cellsegSettings.dumpSettings)
                    scratchList{argPos} = '-dumpsummary'; argPos = argPos + 1;
                end
                if (obj.cellsegSettings.overwrite)
                    scratchList{argPos} = '-ovrw'; argPos = argPos + 1;
                end

                checkTuningParams = true;
                if ~isempty(obj.cellsegSettings.presetName)
                    %Check if preset exists
                    if isfile(['.' filesep 'cellsegTemplates' filesep obj.cellsegSettings.presetName '.mat'])
                        checkTuningParams = false;
                        scratchList{argPos} = '-template'; argPos = argPos + 1;
                        scratchList{argPos} = obj.cellsegSettings.presetName; argPos = argPos + 1;
                    end
                    %Otherwise save the template
                end

                if(checkTuningParams)
                    for j = 1:tuningParamCount
                        paramName = TUNING_PARAMS{j};
                        if (obj.cellsegSettings.(paramName) > 0)
                            scratchList{argPos} = ['-' paramName]; argPos = argPos + 1;
                            scratchList{argPos} = num2str(obj.cellsegSettings.(paramName)); argPos = argPos + 1;
                        end
                    end

                    if ~isempty(obj.cellsegSettings.presetName)
                        scratchList{argPos} = '-savetmpl'; argPos = argPos + 1;
                        scratchList{argPos} = obj.cellsegSettings.presetName; argPos = argPos + 1;
                    end
                end

                %Call the cellseg main
                myArgs = scratchList(1:(argPos-1));
                Main_CellSegConsole(myArgs{:});
            end
        end

        function obj = runSpotDetect(obj, channels, listener)
            if nargin < 3; listener = []; end

            funcTag = 'TrueSpotProjSettings.runSpotDetect';

            if isempty(obj.inputTable)
                GUIUtils.sendMessage(listener, funcTag, 'ERROR: Batch has not been initialized.');
                return;
            end

            GUIUtils.sendMessage(listener, funcTag, 'Starting spot detection...');
            stackCount = size(obj.inputTable, 1);
            channelCount = size(channels, 2);
            for i = 1:stackCount
                GUIUtils.sendMessage(listener, funcTag, sprintf('Working on image %d of %d...', i, stackCount));
                iname = obj.getInputTableValue(i, 'ImageName');
                ipath = obj.getInputTableValue(i, 'ImageFilePath');
                ctrlpath = obj.getInputTableValue(i, 'ControlFilePath');

                iOutDir = [obj.paths.outputPath filesep iname];

                %Import cellseg, if needed
                cellSegPath = [iOutDir filesep 'CellSeg_' iname '.mat'];
                cellMaskPath = obj.getInputTableValue(i, 'ExtCellSegPath');
                nucMaskPath = obj.getInputTableValue(i, 'ExtNucSegPath');
                if ~isempty(cellMaskPath) | ~isempty(nucMaskPath)
                    if ~isempty(listener)
                        listener.logMessage(['[' funcTag '] External cell or nuc mask provided! Converting...']);
                    end

                    CellSeg.importCellSegData(cellMaskPath, nucMaskPath, cellSegPath);
                end

                for c = 1:channelCount
                    GUIUtils.sendMessage(listener, funcTag, sprintf('Working on image %d of %d (Channel %d of %d)...', i, stackCount, c, channelCount));
                    if iscell(channels)
                        channelSpecs = channels{c};
                    else
                        channelSpecs = channels(c);
                    end

                    channelSpecs.runSpotDetect(obj, iname, ipath, ctrlpath, listener);
                end
            end

            %Generate an xml for quant dumps
            xmlpath = [obj.paths.outputPath filesep 'countInfo.xml'];
            fhandle = fopen(xmlpath, 'w');
            fprintf(fhandle, '<?xml version="1.0" encoding="UTF-8"?>\n');
            fprintf(fhandle, '<BatchSamples>\n');
            for c = 1:channelCount
                if iscell(channels)
                    channelSpecs = channels{c};
                else
                    channelSpecs = channels(c);
                end
                if ~isempty(channelSpecs)
                    fprintf(fhandle, '\t<SampleChannel Name="CH%d" ChannelIndex="%d">\n', ...
                        channelSpecs.channelIndex, channelSpecs.channelIndex);
                    fprintf(fhandle, '\t<QueryParams>\n');
                    fprintf(fhandle, '\t\t<QueryParam Key="ChannelNo" Value="%d"/>\n', channelSpecs.channelIndex);
                    fprintf(fhandle, '\t\t<QueryParam Key="INameContains" Value="_CH%d"/>\n', channelSpecs.channelIndex);
                    if ~isempty(channelSpecs.metadata.targetName)
                        fprintf(fhandle, '\t\t<QueryParam Key="TargetName" Value="%s"/>\n', channelSpecs.metadata.targetName);
                    end
                    if ~isempty(channelSpecs.metadata.probeName)
                        fprintf(fhandle, '\t\t<QueryParam Key="ProbeName" Value="%s"/>\n', channelSpecs.metadata.probeName);
                    end
                    fprintf(fhandle, '\t</QueryParams>\n');
                    fprintf(fhandle, '\t\t<CountThresholds ThMid="0" ThLow="0" ThHigh="0"/>\n');
                    fprintf(fhandle, '\t</SampleChannel>\n');
                end
            end
            fprintf(fhandle, '</BatchSamples>\n');
            fclose(fhandle);

            %Call Main_AnalyzeBatchThresholds after all run as well!
            Main_AnalyzeBatchThresholds('-input', obj.paths.outputPath);
        end

        function obj = runQuant(obj, channels, listener)
            if nargin < 3; listener = []; end

            funcTag = 'TrueSpotProjSettings.runQuant';

            if isempty(obj.inputTable)
                GUIUtils.sendMessage(listener, funcTag, 'ERROR: Batch has not been initialized.');
                return;
            end

            GUIUtils.sendMessage(listener, funcTag, 'Starting fit/quantification...');
            stackCount = size(obj.inputTable, 1);
            channelCount = size(channels, 2);
            for i = 1:stackCount
                GUIUtils.sendMessage(listener, funcTag, sprintf('Working on image %d of %d...', i, stackCount));
                iname = obj.getInputTableValue(i, 'ImageName');
                for c = 1:channelCount
                    GUIUtils.sendMessage(listener, funcTag, sprintf('Working on image %d of %d (Channel %d of %d)...', i, stackCount, c, channelCount));
                    if iscell(channels)
                        channelSpecs = channels{c};
                    else
                        channelSpecs = channels(c);
                    end

                    channelSpecs.runQuant(obj, iname, [], [], listener);
                end
            end

            %TODO Quant counts?
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
            varNames = {'ImageName' 'ImageFilePath' 'ControlFilePath' ...
                'NucDataFilePath' 'LightDataFilePath' 'ExtCellSegPath' 'ExtNucSegPath'};
            varTypes = {'string' 'string' 'string' 'string' 'string' 'string' 'string'};
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
            channelInfo.nucStackChannelCount = 0; %If separate file
            channelInfo.lightStackChannelCount = 0; %If separate file
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