%
%%
classdef TrueSpotXML

     %%
    methods (Static)

        function boolres = isAbsolutePath(path)
            if ispc
                boolres = contains(path, ':');
            else
                boolres = startsWith(path, '/');
            end
        end

        function abspath = rel2absPath(relpath, wd)
            abspath = [wd filesep relpath]; %Quick way
            abspath = replace(abspath, '/', filesep);
        end

        %% ========================== Read ==========================
        function textData = getElementText(xmlElement)
            sChild = getFirstChild(xmlElement);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.TEXT_NODE
                    textData = char(sChild.getData());
                    textData = replace(textData, '"', '');
                    return;
                end
                sChild = getNextSibling(sChild);
            end
        end

        function numberVal = getElementTextAsNumber(xmlElement, defoVal)
            if nargin < 2; defoVal = NaN; end
            numberVal = defoVal;
            valStr = TrueSpotXML.getElementText(xmlElement);
            if ~isempty(valStr)
                numberVal = Force2Num(valStr);
            end
        end

        function numberVal = getBoolAttribute(xmlElement, attrKey, defoVal)
            if nargin < 3; defoVal = false; end
            numberVal = defoVal;
            attrStr = char(getAttribute(xmlElement, attrKey));
            if ~isempty(attrStr)
                numberVal = Force2Bool(attrStr);
            end
        end

        function numberVal = getNumberAttribute(xmlElement, attrKey, defoVal)
            if nargin < 3; defoVal = NaN; end
            numberVal = defoVal;
            attrStr = char(getAttribute(xmlElement, attrKey));
            if ~isempty(attrStr)
                numberVal = Force2Num(attrStr);
            end
        end

        function metaStruct = readMetaNode(xmlNode, metaStruct)
            if nargin < 2; metaStruct = []; end
            if isempty(metaStruct)
                metaStruct = TrueSpotChannelSettings.genMetadataStruct();
                metaStruct.species = [];
                metaStruct.cellType = [];
                metaStruct.voxelSize = struct('x', 0, 'y', 0, 'z', 0);
                metaStruct.batchName = [];
            end

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'Species')
                        metaStruct.species = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'CellType')
                        metaStruct.cellType = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'VoxelDimsNano')
                        metaStruct.voxelSize.x = TrueSpotXML.getNumberAttribute(sChild, 'X', 0);
                        metaStruct.voxelSize.y = TrueSpotXML.getNumberAttribute(sChild, 'Y', 0);
                        metaStruct.voxelSize.z = TrueSpotXML.getNumberAttribute(sChild, 'Z', 0);
                    elseif strcmp(sChildName, 'PixelDimsNano')
                        metaStruct.voxelSize.x = TrueSpotXML.getNumberAttribute(sChild, 'X', 0);
                        metaStruct.voxelSize.y = TrueSpotXML.getNumberAttribute(sChild, 'Y', 0);
                        metaStruct.voxelSize.z = 0;
                    elseif strcmp(sChildName, 'PointDimsNano')
                        metaStruct.pointSize.x = TrueSpotXML.getNumberAttribute(sChild, 'X', 0);
                        metaStruct.pointSize.y = TrueSpotXML.getNumberAttribute(sChild, 'Y', 0);
                        metaStruct.pointSize.z = TrueSpotXML.getNumberAttribute(sChild, 'Z', 0);
                    elseif strcmp(sChildName, 'TargetName')
                        metaStruct.targetName = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'ProbeName')
                        metaStruct.probeName = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'TargetType')
                        metaStruct.targetMolType = TrueSpotXML.getElementText(sChild);
                    end
                end
                sChild = getNextSibling(sChild);
            end

        end

        function pathInfo = readPathsNode(xmlNode)
            pathInfo = TrueSpotProjSettings.genPathsStruct();

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'ImageDir')
                        pathInfo.inputPath = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'OutputDir')
                        pathInfo.outputPath = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'Input')
                        pathInfo.inputPath = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'ControlPath')
                        pathInfo.controlPath = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'ExtCellMask')
                        pathInfo.extCellMaskStem = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'ExtNucMask')
                        pathInfo.extNucMaskStem = TrueSpotXML.getElementText(sChild);
                        pathInfo.extNucMaskZMin = TrueSpotXML.getNumberAttribute(sChild, 'ZMin', 1);
                    end
                end
                sChild = getNextSibling(sChild);
            end

        end

        function cellposeSubSettings = readCellposeSubSettingsNode(xmlNode, cellposeSubSettings)
            if nargin < 2; cellposeSettings = []; end
            if isempty(cellposeSubSettings)
                cellposeSubSettings = CellPoseTS.genCellposeSubParamStruct();
            end

            cellposeSubSettings.avg_dia = TrueSpotXML.getNumberAttribute(xmlNode, 'AvgDia', NaN);
            cellposeSubSettings.normalize_bool = TrueSpotXML.getBoolAttribute(xmlNode, 'Normalize', false);

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'Model')
                        cellposeSubSettings.model_name = char(getAttribute(sChild, 'Name'));
                        cellposeSubSettings.ensemble_bool = TrueSpotXML.getBoolAttribute(sChild, 'Ensemble', false);
                    elseif strcmp(sChildName, 'TuningThresholds')
                        cellposeSubSettings.cell_threshold = TrueSpotXML.getNumberAttribute(sChild, 'Cell', 0);
                        cellposeSubSettings.flow_threshold = TrueSpotXML.getNumberAttribute(sChild, 'Flow', 0.4);
                    end
                end
                sChild = getNextSibling(sChild);
            end
        end

        function cellposeSettings = readCellposeSettingsNode(xmlNode, cellposeSettings)
            if nargin < 2; cellposeSettings = []; end
            if isempty(cellposeSettings)
                cellposeSettings = CellPoseTS.genCellposeParamStruct();
            end

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'NucSettings')
                        cellposeSettings.nuc_params = TrueSpotXML.readCellposeSubSettingsNode(sChild, cellposeSettings.nuc_params);
                    elseif strcmp(sChildName, 'CytoSettings')
                        cellposeSettings.cyto_params = TrueSpotXML.readCellposeSubSettingsNode(sChild, cellposeSettings.cyto_params);
                    end
                end
                sChild = getNextSibling(sChild);
            end
        end

        function cellsegSettings = readCellSegSettingsNode(xmlNode, cellsegSettings)
            if nargin < 2; cellsegSettings = []; end
            if isempty(cellsegSettings)
                cellsegSettings = TrueSpotProjSettings.genCellsegSettingsStruct();
            end

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'PresetName')
                        cellsegSettings.presetName = TrueSpotXML.getElementText(sChild);
                    elseif strcmp(sChildName, 'TransZTrim')
                        cellsegSettings.lightZMin = TrueSpotXML.getNumberAttribute(sChild, 'Min', 0);
                        cellsegSettings.lightZMax = TrueSpotXML.getNumberAttribute(sChild, 'Max', 0);
                    elseif strcmp(sChildName, 'NucZTrim')
                        cellsegSettings.nucZMin = TrueSpotXML.getNumberAttribute(sChild, 'Min', 0);
                        cellsegSettings.nucZMax = TrueSpotXML.getNumberAttribute(sChild, 'Max', 0);
                    elseif strcmp(sChildName, 'CellSize')
                        cellsegSettings.cszmin = TrueSpotXML.getNumberAttribute(sChild, 'Min', 0);
                        cellsegSettings.cszmax = TrueSpotXML.getNumberAttribute(sChild, 'Max', 0);
                    elseif strcmp(sChildName, 'NucSize')
                        cellsegSettings.nszmin = TrueSpotXML.getNumberAttribute(sChild, 'Min', 0);
                        cellsegSettings.nszmax = TrueSpotXML.getNumberAttribute(sChild, 'Max', 0);
                    elseif strcmp(sChildName, 'XTrim')
                        cellsegSettings.xtrim = TrueSpotXML.getElementTextAsNumber(sChild, 0);
                    elseif strcmp(sChildName, 'YTrim')
                        cellsegSettings.ytrim = TrueSpotXML.getElementTextAsNumber(sChild, 0);
                    elseif strcmp(sChildName, 'NucZRange')
                        cellsegSettings.nzrange = TrueSpotXML.getElementTextAsNumber(sChild, 0);
                    elseif strcmp(sChildName, 'NucThSample')
                        cellsegSettings.nthsmpl = TrueSpotXML.getElementTextAsNumber(sChild, 0);
                    elseif strcmp(sChildName, 'NucCutoff')
                        cellsegSettings.ncutoff = TrueSpotXML.getElementTextAsNumber(sChild, 0);
                    elseif strcmp(sChildName, 'NucDXY')
                        cellsegSettings.ndxy = TrueSpotXML.getElementTextAsNumber(sChild);
                    elseif strcmp(sChildName, 'Options')
                        attrStr = char(getAttribute(sChild, 'ExportCellMaskToFormat'));
                        if ~isempty(attrStr)
                            if strcmp(attrStr, 'png')
                                cellsegSettings.outputCellMaskPNG = true;
                            elseif strcmp(attrStr, 'tif')
                                cellsegSettings.outputCellMaskTIF = true;
                            end
                        end

                        attrStr = char(getAttribute(sChild, 'ExportNucMaskToFormat'));
                        if ~isempty(attrStr)
                            if strcmp(attrStr, 'png')
                                cellsegSettings.outputNucMaskPNG = true;
                            elseif strcmp(attrStr, 'tif')
                                cellsegSettings.outputNucMaskTIF = true;
                            end
                        end

                        cellsegSettings.overwrite = TrueSpotXML.getBoolAttribute(sChild, 'Overwrite', true);
                        cellsegSettings.dumpSettings = TrueSpotXML.getBoolAttribute(sChild, 'DumpSettingsToText', false);
                    elseif strcmp(sChildName, 'CellposeSettings')
                        cellsegSettings.useCellposeNuc = TrueSpotXML.getBoolAttribute(sChild, 'UseCellposeNuc', true);
                        cellsegSettings.useCellposeCyto = TrueSpotXML.getBoolAttribute(sChild, 'UseCellposeCyto', true);
                        cellsegSettings.cellpose = TrueSpotXML.readCellposeSettingsNode(sChild, cellsegSettings.cellpose);
                    end
                end
                sChild = getNextSibling(sChild);
            end
        end

        function thresholdSettings = readThreshSettingsNode(xmlNode)
            thresholdSettings = TrueSpotChannelSettings.genThresholdSettingsStruct();

            thresholdSettings.preset = TrueSpotXML.getNumberAttribute(xmlNode, 'Preset', NaN);
            thresholdSettings.thMin = TrueSpotXML.getNumberAttribute(xmlNode, 'ScanMin', 0);
            thresholdSettings.thMax = TrueSpotXML.getNumberAttribute(xmlNode, 'ScanMax', 0);

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'WindowSettings')
                        winMin = TrueSpotXML.getNumberAttribute(sChild, 'Min', 0.006);
                        winMax = TrueSpotXML.getNumberAttribute(sChild, 'Max', 0.042);
                        winIncr = TrueSpotXML.getNumberAttribute(sChild, 'Increment', 0.006);
                        thresholdSettings.thParams.window_sizes = winMin:winIncr:winMax;
                        clear winMin winMax winIncr
                    elseif strcmp(sChildName, 'MADFactor')
                        thresholdSettings.thParams.mad_factor_min = TrueSpotXML.getNumberAttribute(sChild, 'Min', -1.0);
                        thresholdSettings.thParams.mad_factor_max = TrueSpotXML.getNumberAttribute(sChild, 'Max', 1.0);
                    elseif strcmp(sChildName, 'Weights')
                        thresholdSettings.thParams.fit_ri_weight = TrueSpotXML.getNumberAttribute(sChild, 'FitRightIntersect', 0.0);
                        thresholdSettings.thParams.madth_weight = TrueSpotXML.getNumberAttribute(sChild, 'MedMad', 0.0);
                        thresholdSettings.thParams.fit_weight = TrueSpotXML.getNumberAttribute(sChild, 'Fit', 1.0);
                    elseif strcmp(sChildName, 'MiscOptions')
                        thresholdSettings.thParams.test_data = TrueSpotXML.getBoolAttribute(sChild, 'IncludeRawCurve', false);
                        thresholdSettings.thParams.test_diff = TrueSpotXML.getBoolAttribute(sChild, 'IncludeDiffCurve', false);
                        thresholdSettings.thParams.std_factor = TrueSpotXML.getNumberAttribute(sChild, 'StDevFactor', 0.0);
                        attrStr = char(getAttribute(sChild, 'LogMode'));
                        if ~isempty(attrStr)
                            if strcmp(attrStr, 'None')
                                thresholdSettings.thParams.log_proj_mode = 0;
                            elseif strcmp(attrStr, 'All')
                                thresholdSettings.thParams.log_proj_mode = 1;
                            elseif strcmp(attrStr, 'FitOnly')
                                thresholdSettings.thParams.log_proj_mode = 2;
                            end
                        end
                    end
                end
                sChild = getNextSibling(sChild);
            end
        end

        function [spotSettings, thresholdSettings, optionsStruct] = readSpotDetectSettingsNode(xmlNode, spotSettings, thresholdSettings, optionsStruct)
            if nargin < 2; spotSettings = []; end
            if nargin < 3; thresholdSettings = []; end
            if nargin < 4; optionsStruct = []; end
            
            if isempty(spotSettings)
                spotSettings = TrueSpotChannelSettings.genCountSettingsStruct();
            end
            if isempty(thresholdSettings)
                thresholdSettings = TrueSpotChannelSettings.genThresholdSettingsStruct();
            end
            if isempty(optionsStruct)
                optionsStruct = TrueSpotChannelSettings.genOptionsStruct();
            end

            spotSettings.gaussRad = TrueSpotXML.getNumberAttribute(xmlNode, 'GaussRad', 7);
            spotSettings.useDPC = ~TrueSpotXML.getBoolAttribute(xmlNode, 'NoDPC', false);
            spotSettings.spotDetectThreads = TrueSpotXML.getNumberAttribute(xmlNode, 'Workers', 1);

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'ThresholdSettings')
                        thresholdSettings = TrueSpotXML.readThreshSettingsNode(sChild);
                    elseif strcmp(sChildName, 'ZTrim')
                        spotSettings.zMin = TrueSpotXML.getNumberAttribute(sChild, 'Min', 0);
                        spotSettings.zMax = TrueSpotXML.getNumberAttribute(sChild, 'Max', 0);
                    elseif strcmp(sChildName, 'Options')
                        optionsStruct.runparamtxt = TrueSpotXML.getBoolAttribute(sChild, 'DumpRunParamsToText', false);
                        optionsStruct.overwrite = TrueSpotXML.getBoolAttribute(sChild, 'Overwrite', true);
                        optionsStruct.csvzero = TrueSpotXML.getBoolAttribute(sChild, 'CsvZeroBased', false);
                        attrStr = char(getAttribute(sChild, 'CsvDump'));
                        if ~isempty(attrStr)
                            if strcmp(attrStr, 'All')
                                optionsStruct.csvfull = true;
                            elseif strcmp(attrStr, 'ThresholdRange')
                                optionsStruct.csvrange = true;
                            elseif strcmp(attrStr, 'SelectedThreshold')
                                optionsStruct.csvthonly = true;
                            end
                        end
                    end
                end
                sChild = getNextSibling(sChild);
            end
        end

        function [spotSettings, manualThresh] = readQuantSettingsNode(xmlNode, spotSettingsStruct)
            if ~isempty(spotSettingsStruct)
                spotSettings = spotSettingsStruct;
            else
                spotSettings = TrueSpotChannelSettings.genCountSettingsStruct();
            end

            spotSettings.quantNoClouds = ~TrueSpotXML.getBoolAttribute(xmlNode, 'DoClouds', false);
            spotSettings.quantCellZero = TrueSpotXML.getBoolAttribute(xmlNode, 'CellZero', false);
            spotSettings.quantThreads = TrueSpotXML.getNumberAttribute(xmlNode, 'Workers', 1);
            spotSettings.zAdj = TrueSpotXML.getNumberAttribute(xmlNode, 'ZtoXY', NaN);
            spotSettings.fitRadXY = TrueSpotXML.getNumberAttribute(xmlNode, 'FitRadXY', 4);
            spotSettings.fitRadZ = TrueSpotXML.getNumberAttribute(xmlNode, 'FitRadZ', 2);
            spotSettings.quantNoRefilter = ~TrueSpotXML.getBoolAttribute(xmlNode, 'DoRefilter', false);

            manualThresh = ~TrueSpotXML.getNumberAttribute(xmlNode, 'ManualThreshold', 0);
        end

        function channelSettings = readImageChannelNode(xmlNode, spotCommon, threshCommon, opsCommon)
            channelSettings = TrueSpotChannelSettings.newChannelSettings();
            
            if ~isempty(spotCommon)
                channelSettings.spotCountSettings = spotCommon;
            end
            if ~isempty(threshCommon)
                channelSettings.thresholdSettings = threshCommon;
            end
            if ~isempty(opsCommon)
                channelSettings.options = opsCommon;
            end

            channelSettings.channelIndex = TrueSpotXML.getNumberAttribute(xmlNode, 'ChannelNumber', 0);
            channelSettings.spotCountSettings.controlImageChannel = TrueSpotXML.getNumberAttribute(xmlNode, 'ControlChannelNumber', 0);

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'Meta')
                        channelSettings.metadata = TrueSpotXML.readMetaNode(sChild, channelSettings.metadata);
                    elseif strcmp(sChildName, 'SpotDetectSettings')
                        [channelSettings.spotCountSettings, channelSettings.thresholdSettings, channelSettings.options] = ...
                            TrueSpotXML.readSpotDetectSettingsNode(sChild, channelSettings.spotCountSettings, channelSettings.thresholdSettings, channelSettings.options);
                    elseif strcmp(sChildName, 'QuantSettings')
                        [channelSettings.spotCountSettings, manualThresh] = TrueSpotXML.readQuantSettingsNode(sChild, channelSettings.spotCountSettings);
                        if manualThresh > 0
                            channelSettings.thresholdSettings.manualTh = manualThresh;
                        end
                    end
                end
                sChild = getNextSibling(sChild);
            end
        end

        function [channels, channelInfo] = readBatchChannelInfo(xmlNode, channelInfo, spotCommon, threshCommon, opsCommon)
            if nargin < 2; channelInfo = []; end
            if nargin < 3; spotCommon = []; end
            if nargin < 4; threshCommon = []; end
            if nargin < 5; opsCommon = []; end

            channelInfo.channelCount = TrueSpotXML.getNumberAttribute(xmlNode, 'ChannelCount', 1);
            channelInfo.nucMarkerChannel = TrueSpotXML.getNumberAttribute(xmlNode, 'NucChannel', 0);
            channelInfo.lightChannel = TrueSpotXML.getNumberAttribute(xmlNode, 'TransChannel', 0);
            
            channelInfo.controlChannelCount = TrueSpotXML.getNumberAttribute(xmlNode, 'ControlChannelCount', 0);
            channelInfo.controlNucMarkerChannel = TrueSpotXML.getNumberAttribute(xmlNode, 'ControlNucChannel', 0);
            channelInfo.controlLightChannel = TrueSpotXML.getNumberAttribute(xmlNode, 'ControlTransChannel', 0);

            chElements = getElementsByTagName(xmlNode,'ImageChannel');
            chCount = chElements.getLength;

            channels = cell(1, chCount);
            for ci = 1:chCount
                ciz = ci-1;
                chNode = item(chElements, ciz);
                channelSettings = TrueSpotXML.readImageChannelNode(chNode, spotCommon, threshCommon, opsCommon);
                channels{ci} = channelSettings;
            end
        end

        function [batchSettings, channels] = readImageBatchNode(xmlNode, metaCommon, cellsegCommon, spotCommon, threshCommon, opsCommon)
            batchSettings = TrueSpotProjSettings.newBatchSettings();
            
            if ~isempty(metaCommon)
                batchSettings.metadata = metaCommon;
            end
            if ~isempty(cellsegCommon)
                batchSettings.cellsegSettings = cellsegCommon;
            end

            batchSettings.metadata.batchName = char(getAttribute(xmlNode, 'Name'));

            sChild = getFirstChild(xmlNode);
            while ~isempty(sChild)
                if sChild.getNodeType == sChild.ELEMENT_NODE
                    sChildName = char(sChild.getTagName);
                    if strcmp(sChildName, 'Meta')
                        batchSettings.metadata = TrueSpotXML.readMetaNode(sChild, batchSettings.metadata);
                    elseif strcmp(sChildName, 'Paths')
                        batchSettings.paths = TrueSpotXML.readPathsNode(sChild);
                    elseif strcmp(sChildName, 'CellSegSettings')
                        batchSettings.cellsegSettings = TrueSpotXML.readCellSegSettingsNode(sChild, batchSettings.cellsegSettings);
                    elseif strcmp(sChildName, 'ChannelInfo')
                        [channels, batchSettings.channelInfo] = ...
                            TrueSpotXML.readBatchChannelInfo(sChild, batchSettings.channelInfo, spotCommon, threshCommon, opsCommon);
                    end
                end
                sChild = getNextSibling(sChild);
            end
        end

        function batchList = readSettingsXML(xmlpath)
            rootnode = xmlread(xmlpath);
            %Extract common parameters from "ImageSet"
            childList = getElementsByTagName(rootnode,'ImageSet');
            childCount = childList.getLength;

            imgSetNode = item(childList, 0);
            batchElements = getElementsByTagName(imgSetNode,'ImageBatch');
            batchCount = batchElements.getLength;

            batchStruct = struct('batchSettings', [], 'channels', []);

            batchList(batchCount) = batchStruct;
            for i = (batchCount-1):-1:1
                batchList(i) = batchStruct;
            end
            clear batchStruct i batchElements

            metaCommon = TrueSpotProjSettings.genMetadataStruct();
            cellsegCommon = [];
            spotCommon = [];
            threshCommon = [];
            opsCommon = [];
            ii = 1;
            for i = 1:childCount
                iz = i-1;
                childNode = item(childList, iz);

                gChild = getFirstChild(childNode);
                while ~isempty(gChild)
                    if gChild.getNodeType == gChild.ELEMENT_NODE
                        gChildName = char(gChild.getTagName);
                        if strcmp(gChildName, 'ImageBatch')
                            batchStruct = batchList(ii);
                            [batchStruct.batchSettings, batchStruct.channels] = TrueSpotXML.readImageBatchNode(gChild, metaCommon, cellsegCommon, spotCommon, threshCommon, opsCommon);
                            batchList(ii) = batchStruct;
                            ii = ii + 1;
                            clear batchStruct
                        elseif strcmp(gChildName, 'CommonMeta')
                            metaCommon = TrueSpotXML.readMetaNode(gChild, metaCommon);
                        elseif strcmp(gChildName, 'CellSegSettings')
                            cellsegCommon = TrueSpotXML.readCellSegSettingsNode(gChild, cellsegCommon);
                        elseif strcmp(gChildName, 'SpotDetectSettings')
                            [spotCommon, threshCommon, opsCommon] = ...
                                TrueSpotXML.readSpotDetectSettingsNode(gChild, spotCommon, threshCommon, opsCommon);
                        elseif strcmp(gChildName, 'QuantSettings')
                            [spotCommon, ~] = TrueSpotXML.readQuantSettingsNode(gChild, spotCommon);
                        end
                    end
                    gChild = getNextSibling(gChild);
                end
            end

            %Detect any relative paths and convert to absolute
            [xmldir, ~, ~] = fileparts(xmlpath);
            for i = 1:batchCount
                batchStruct = batchList(i);
                if ~isempty(batchStruct.batchSettings.paths.inputPath)
                    if ~TrueSpotXML.isAbsolutePath(batchStruct.batchSettings.paths.inputPath)
                        batchStruct.batchSettings.paths.inputPath = TrueSpotXML.rel2absPath(batchStruct.batchSettings.paths.inputPath, xmldir);
                    end
                end

                if ~isempty(batchStruct.batchSettings.paths.controlPath)
                    if ~TrueSpotXML.isAbsolutePath(batchStruct.batchSettings.paths.controlPath)
                        batchStruct.batchSettings.paths.controlPath = TrueSpotXML.rel2absPath(batchStruct.batchSettings.paths.controlPath, xmldir);
                    end
                end

                if ~isempty(batchStruct.batchSettings.paths.outputPath)
                    if ~TrueSpotXML.isAbsolutePath(batchStruct.batchSettings.paths.outputPath)
                        batchStruct.batchSettings.paths.outputPath = TrueSpotXML.rel2absPath(batchStruct.batchSettings.paths.outputPath, xmldir);
                    end
                end

                if ~isempty(batchStruct.batchSettings.paths.extCellMaskStem)
                    if ~TrueSpotXML.isAbsolutePath(batchStruct.batchSettings.paths.extCellMaskStem)
                        batchStruct.batchSettings.paths.extCellMaskStem = TrueSpotXML.rel2absPath(batchStruct.batchSettings.paths.extCellMaskStem, xmldir);
                    end
                end

                if ~isempty(batchStruct.batchSettings.paths.extNucMaskStem)
                    if ~TrueSpotXML.isAbsolutePath(batchStruct.batchSettings.paths.extNucMaskStem)
                        batchStruct.batchSettings.paths.extNucMaskStem = TrueSpotXML.rel2absPath(batchStruct.batchSettings.paths.extNucMaskStem, xmldir);
                    end
                end

                batchList(i) = batchStruct;
            end

        end

        %% ========================== Write ==========================

        function strval = bool2str(boolVal)
            strval = 'True';
            if boolVal; return; end
            strval = 'False';
        end

        function writeBatchMetaBlock(xmlHandle, metaInfo, indent)
            fprintf(xmlHandle, '%s<CommonMeta>\n', indent);
            if ~isempty(metaInfo.species)
                fprintf(xmlHandle, '%s\t<Species>"%s"</Species>\n', indent, metaInfo.species);
            end
            if ~isempty(metaInfo.cellType)
                fprintf(xmlHandle, '%s\t<CellType>"%s"</CellType>\n', indent, metaInfo.cellType);
            end
            if ~isempty(metaInfo.voxelSize)
                fprintf(xmlHandle, '%s\t<VoxelDimsNano X="%d" Y="%d" Z="%d"/>\n',...
                    indent, metaInfo.voxelSize.x, metaInfo.voxelSize.y, metaInfo.voxelSize.z);
            end
            fprintf(xmlHandle, '%s</CommonMeta>\n', indent);
        end

        function writeChannelMetaBlock(xmlHandle, metaInfo, indent)
            fprintf(xmlHandle, '%s<Meta>\n', indent);
            if ~isempty(metaInfo.targetName)
                fprintf(xmlHandle, '%s\t<TargetName>"%s"</TargetName>\n', indent, metaInfo.targetName);
            end
            if ~isempty(metaInfo.probeName)
                fprintf(xmlHandle, '%s\t<ProbeName>"%s"</ProbeName>\n', indent, metaInfo.probeName);
            end
            if ~isempty(metaInfo.targetMolType)
                fprintf(xmlHandle, '%s\t<TargetType>"%s"</TargetType>\n', indent, metaInfo.targetMolType);
            end
            if ~isempty(metaInfo.pointSize)
                fprintf(xmlHandle, '%s\t<PointDimsNano X="%d" Y="%d" Z="%d"/>\n',...
                    indent, metaInfo.pointSize.x, metaInfo.pointSize.y, metaInfo.pointSize.z);
            end
            fprintf(xmlHandle, '%s</Meta>\n', indent);
        end

        function writePathsBlock(xmlHandle, pathsInfo, indent)
            fprintf(xmlHandle, '%s<Paths>\n', indent);
            if ~isempty(pathsInfo.inputPath)
                fprintf(xmlHandle, '%s\t<Input>"%s"</Input>\n', indent, pathsInfo.inputPath);
            end
            if ~isempty(pathsInfo.outputPath)
                fprintf(xmlHandle, '%s\t<OutputDir>"%s"</OutputDir>\n', indent, pathsInfo.outputPath);
            end
            if ~isempty(pathsInfo.controlPath)
                fprintf(xmlHandle, '%s\t<ControlPath>"%s"</ControlPath>\n', indent, pathsInfo.controlPath);
            end
            if ~isempty(pathsInfo.extCellMaskStem)
                fprintf(xmlHandle, '%s\t<ExtCellMask>"%s"</ExtCellMask>\n', indent, pathsInfo.extCellMaskStem);
            end
            if ~isempty(pathsInfo.extNucMaskStem)
                fprintf(xmlHandle, '%s\t<ExtNucMask ZMin="%d">"%s"</ExtNucMask>\n', indent, pathsInfo.extNucMaskStem, pathsInfo.extNucMaskZMin);
            end
            fprintf(xmlHandle, '%s</Paths>\n', indent);
        end

        function writeChannelInfoBlock(xmlHandle, channelInfo, channelSettingsList, indent)
            fprintf(xmlHandle, '%s<ChannelInfo', indent);
            fprintf(xmlHandle, ' ChannelCount="%d"', channelInfo.channelCount);
            fprintf(xmlHandle, ' NucChannel="%d"', channelInfo.nucMarkerChannel);
            fprintf(xmlHandle, ' TransChannel="%d"', channelInfo.lightChannel);
            if channelInfo.controlChannelCount > 0
                fprintf(xmlHandle, ' ControlChannelCount="%d"', channelInfo.controlChannelCount);
            end
            if channelInfo.controlNucMarkerChannel > 0
                fprintf(xmlHandle, ' ControlNucChannel="%d"', channelInfo.controlNucMarkerChannel);
            end
            if channelInfo.controlLightChannel > 0
                fprintf(xmlHandle, ' ControlTransChannel="%d"', channelInfo.controlLightChannel);
            end
            fprintf(xmlHandle, '>\n');

            cCount = size(channelSettingsList, 2);
            for c = 1:cCount
                if iscell(channelSettingsList)
                    myChannel = channelSettingsList{c};
                else
                    myChannel = channelSettingsList(c);
                end
                
                TrueSpotXML.writeImageChannelBlock(xmlHandle, myChannel, [indent '\t']);
            end

            fprintf(xmlHandle, '%s</ChannelInfo>\n', indent);
        end

        function writeCellposeBlock(xmlHandle, cellsegSettings, indent)
            fprintf(xmlHandle, '%s<CellposeSettings', indent);
            if cellsegSettings.useCellposeNuc
                fprintf(xmlHandle, ' UseCellposeNuc="True"');
            else
                fprintf(xmlHandle, ' UseCellposeNuc="False"');
            end
            if cellsegSettings.useCellposeCyto
                fprintf(xmlHandle, ' UseCellposeCyto="True"');
            else
                fprintf(xmlHandle, ' UseCellposeCyto="False"');
            end
            fprintf(xmlHandle, '>\n');

            if cellsegSettings.useCellposeNuc
                fprintf(xmlHandle, '%s\t<NucSettings', indent);
                fprintf(xmlHandle, '>\n');

                if ~isnan(cellsegSettings.cellpose.nuc_params.avg_dia)
                    fprintf(xmlHandle, ' AvgDia="%f"', cellsegSettings.cellpose.nuc_params.avg_dia);
                end

                if cellsegSettings.cellpose.nuc_params.normalize_bool
                    fprintf(xmlHandle, ' Normalize="True"');
                else
                    fprintf(xmlHandle, ' Normalize="False"');
                end

                fprintf(xmlHandle, '%s\t\t<Model Name="%s"', indent, cellsegSettings.cellpose.nuc_params.model_name);
                if cellsegSettings.cellpose.nuc_params.ensemble_bool
                    fprintf(xmlHandle, ' Ensemble="True"');
                else
                    fprintf(xmlHandle, ' Ensemble="False"');
                end
                fprintf(xmlHandle, '/>\n');

                fprintf(xmlHandle, '%s\t\t<TuningThresholds Cell="%f" Flow="%f"/>', ...
                    indent, cellsegSettings.cellpose.nuc_params.cell_threshold, cellsegSettings.cellpose.nuc_params.flow_threshold);

                fprintf(xmlHandle, '%s\t</NucSettings>\n', indent);
            end

            if cellsegSettings.useCellposeCyto
                fprintf(xmlHandle, '%s\t<CytoSettings', indent);
                fprintf(xmlHandle, '>\n');

                if ~isnan(cellsegSettings.cellpose.nuc_params.avg_dia)
                    fprintf(xmlHandle, ' AvgDia="%f"', cellsegSettings.cellpose.nuc_params.avg_dia);
                end

                if cellsegSettings.cellpose.nuc_params.normalize_bool
                    fprintf(xmlHandle, ' Normalize="True"');
                else
                    fprintf(xmlHandle, ' Normalize="False"');
                end

                fprintf(xmlHandle, '%s\t\t<Model Name="%s"', indent, cellsegSettings.cellpose.cyto_params.model_name);
                if cellsegSettings.cellpose.cyto_params.ensemble_bool
                    fprintf(xmlHandle, ' Ensemble="True"');
                else
                    fprintf(xmlHandle, ' Ensemble="False"');
                end
                fprintf(xmlHandle, '/>\n');

                fprintf(xmlHandle, '%s\t\t<TuningThresholds Cell="%f" Flow="%f"/>', ...
                    indent, cellsegSettings.cellpose.cyto_params.cell_threshold, cellsegSettings.cellpose.cyto_params.flow_threshold);

                fprintf(xmlHandle, '%s\t</CytoSettings>\n', indent);
            end

            fprintf(xmlHandle, '%s</CellposeSettings>\n', indent);
        end

        function writeCellSegBlock(xmlHandle, cellsegSettings, indent)
            fprintf(xmlHandle, '%s<CellSegSettings>\n', indent);
            if ~isempty(cellsegSettings.presetName)
                fprintf(xmlHandle, '%s\t<PresetName>"%s"</PresetName>\n', indent, cellsegSettings.presetName);
            end
            if (cellsegSettings.lightZMin > 0) | (cellsegSettings.lightZMax > 0)
                fprintf(xmlHandle, '%s\t<TransZTrim', indent);
                if (cellsegSettings.lightZMin > 0)
                    fprintf(xmlHandle, ' Min="%d"', cellsegSettings.lightZMin);
                end
                if (cellsegSettings.lightZMax > 0)
                    fprintf(xmlHandle, ' Max="%d"', cellsegSettings.lightZMax);
                end
                fprintf(xmlHandle, '/>\n');
            end
            if (cellsegSettings.nucZMin > 0) | (cellsegSettings.nucZMax > 0)
                fprintf(xmlHandle, '%s\t<NucZTrim', indent);
                if (cellsegSettings.nucZMin > 0)
                    fprintf(xmlHandle, ' Min="%d"', cellsegSettings.nucZMin);
                end
                if (cellsegSettings.nucZMax > 0)
                    fprintf(xmlHandle, ' Max="%d"', cellsegSettings.nucZMax);
                end
                fprintf(xmlHandle, '/>\n');
            end
            if (cellsegSettings.cszmin > 0) | (cellsegSettings.cszmax > 0)
                fprintf(xmlHandle, '%s\t<CellSize', indent);
                if (cellsegSettings.cszmin > 0)
                    fprintf(xmlHandle, ' Min="%d"', cellsegSettings.cszmin);
                end
                if (cellsegSettings.cszmax > 0)
                    fprintf(xmlHandle, ' Max="%d"', cellsegSettings.cszmax);
                end
                fprintf(xmlHandle, '/>\n');
            end
            if (cellsegSettings.nszmin > 0) | (cellsegSettings.nszmax > 0)
                fprintf(xmlHandle, '%s\t<NucSize', indent);
                if (cellsegSettings.nszmin > 0)
                    fprintf(xmlHandle, ' Min="%d"', cellsegSettings.nszmin);
                end
                if (cellsegSettings.nszmax > 0)
                    fprintf(xmlHandle, ' Max="%d"', cellsegSettings.nszmax);
                end
                fprintf(xmlHandle, '/>\n');
            end
            if cellsegSettings.xtrim > 0
                fprintf(xmlHandle, '%s\t<XTrim>"%d"</XTrim>\n', indent, cellsegSettings.xtrim);
            end
            if cellsegSettings.ytrim > 0
                fprintf(xmlHandle, '%s\t<YTrim>"%d"</YTrim>\n', indent, cellsegSettings.ytrim);
            end
            if cellsegSettings.nzrange > 0
                fprintf(xmlHandle, '%s\t<NucZRange>"%d"</NucZRange>\n', indent, cellsegSettings.nzrange);
            end
            if cellsegSettings.nthsmpl > 0
                fprintf(xmlHandle, '%s\t<NucThSample>"%d"</NucThSample>\n', indent, cellsegSettings.nthsmpl);
            end
            if cellsegSettings.ncutoff > 0
                fprintf(xmlHandle, '%s\t<NucCutoff>"%f"</NucCutoff>\n', indent, cellsegSettings.ncutoff);
            end
            if cellsegSettings.ndxy > 0
                fprintf(xmlHandle, '%s\t<NucDXY>"%f"</NucDXY>\n', indent, cellsegSettings.ndxy);
            end

            fprintf(xmlHandle, '%s\t<Options', indent);
            if cellsegSettings.outputCellMaskPNG
                fprintf(xmlHandle, ' ExportCellMaskToFormat="png"');
            end
            if cellsegSettings.outputCellMaskTIF
                fprintf(xmlHandle, ' ExportCellMaskToFormat="tif"');
            end
            if cellsegSettings.outputNucMaskPNG
                fprintf(xmlHandle, ' ExportNucMaskToFormat="png"');
            end
            if cellsegSettings.outputNucMaskTIF
                fprintf(xmlHandle, ' ExportNucMaskToFormat="tif"');
            end
            if cellsegSettings.overwrite
                fprintf(xmlHandle, ' Overwrite="True"');
            end
            if cellsegSettings.dumpSettings
                fprintf(xmlHandle, ' DumpSettingsToText="True"');
            end
            fprintf(xmlHandle, '/>\n');

            if cellsegSettings.useCellposeNuc | cellsegSettings.useCellposeCyto
                if ~isempty(cellsegSettings.cellpose)
                    TrueSpotXML.writeCellposeBlock(xmlHandle, cellsegSettings, [indent '\t']);
                end
            end

            fprintf(xmlHandle, '%s</CellSegSettings>\n', indent);
        end

        function writeThresholdSettingsBlock(xmlHandle, thSettings, indent)
            fprintf(xmlHandle, '%s<ThresholdSettings', indent);

            fprintf(xmlHandle, ' Preset="%d"', thSettings.preset);
            fprintf(xmlHandle, ' ScanMin="%d"', thSettings.thMin);
            fprintf(xmlHandle, ' ScanMax="%d"', thSettings.thMax);
            fprintf(xmlHandle, '>\n');

            if ~isempty(thSettings.thParams)
                if ~isempty(thSettings.thParams.window_sizes)
                    fprintf(xmlHandle, '%s\t<WindowSettings', indent);
                    fprintf(xmlHandle, ' Min="%d"', min(thSettings.thParams.window_sizes, [], 'all'));
                    fprintf(xmlHandle, ' Max="%d"', max(thSettings.thParams.window_sizes, [], 'all'));
                    fprintf(xmlHandle, ' Increment="%d"', min(diff(thSettings.thParams.window_sizes), [], 'all'));
                    fprintf(xmlHandle, '/>\n');
                end

                fprintf(xmlHandle, '%s\t<MADFactor', indent);
                fprintf(xmlHandle, ' Min="%d"', thSettings.thParams.mad_factor_min);
                fprintf(xmlHandle, ' Max="%d"', thSettings.thParams.mad_factor_max);
                fprintf(xmlHandle, '/>\n');

                fprintf(xmlHandle, '%s\t<Weights', indent);
                fprintf(xmlHandle, ' FitRightIntersect="%f"', thSettings.thParams.fit_ri_weight);
                fprintf(xmlHandle, ' MedMad="%f"', thSettings.thParams.madth_weight);
                fprintf(xmlHandle, ' Fit="%f"', thSettings.thParams.fit_weight);
                fprintf(xmlHandle, '/>\n');

                fprintf(xmlHandle, '%s\t<MiscOptions', indent);
                fprintf(xmlHandle, ' IncludeRawCurve="%s"', TrueSpotXML.bool2str(thSettings.thParams.test_data));
                fprintf(xmlHandle, ' IncludeDiffCurve="%s"', TrueSpotXML.bool2str(thSettings.thParams.test_diff));
                fprintf(xmlHandle, ' StDevFactor="%f"', thSettings.thParams.std_factor);
                fprintf(xmlHandle, ' LogMode=');
                if thSettings.thParams.log_proj_mode == 0
                    fprintf(xmlHandle, '"None"');
                elseif thSettings.thParams.log_proj_mode == 1
                    fprintf(xmlHandle, '"All"');
                elseif thSettings.thParams.log_proj_mode == 2
                    fprintf(xmlHandle, '"FitOnly"');
                end
                fprintf(xmlHandle, '/>\n');
            end

            fprintf(xmlHandle, '%s</ThresholdSettings>\n', indent);
        end

        function writeSpotsBlock(xmlHandle, spotsSettings, optionsStruct, thSettings, indent)
            fprintf(xmlHandle, '%s<SpotDetectSettings', indent);
            fprintf(xmlHandle, ' GaussRad="%d"', spotsSettings.gaussRad);
            fprintf(xmlHandle, ' Workers="%d"', spotsSettings.spotDetectThreads);
            if spotsSettings.useDPC
                fprintf(xmlHandle, ' NoDPC="False"');
            else
                fprintf(xmlHandle, ' NoDPC="True"');
            end

            if (spotsSettings.zMin > 0) | (spotsSettings.zMax > 0)
                fprintf(xmlHandle, '%s\t<ZTrim', indent);
                if (spotsSettings.zMin > 0)
                    fprintf(xmlHandle, ' Min="%d"', spotsSettings.zMin);
                end
                if (spotsSettings.zMax > 0)
                    fprintf(xmlHandle, ' Max="%d"', spotsSettings.zMax);
                end
                fprintf(xmlHandle, '/>\n');
            end

            fprintf(xmlHandle, '%s\t<Options', indent);
            if optionsStruct.runparamtxt
                fprintf(xmlHandle, ' DumpRunParamsToText="True"');
            else
                fprintf(xmlHandle, ' DumpRunParamsToText="False"');
            end
            if optionsStruct.overwrite
                fprintf(xmlHandle, ' Overwrite="True"');
            else
                fprintf(xmlHandle, ' Overwrite="False"');
            end
            if optionsStruct.csvrange
                fprintf(xmlHandle, ' CsvDump="ThresholdRange"');
            elseif optionsStruct.csvthonly
                fprintf(xmlHandle, ' CsvDump="SelectedThreshold"');
            elseif optionsStruct.csvfull
                fprintf(xmlHandle, ' CsvDump="All"');
            end
            if optionsStruct.csvzero
                fprintf(xmlHandle, ' CsvZeroBased="True"');
            end
            fprintf(xmlHandle, '/>\n');

            if ~isempty(thSettings)
                TrueSpotXML.writeThresholdSettingsBlock(xmlHandle, thSettings, [indent '\t']);
            end

            fprintf(xmlHandle, '%s</SpotDetectSettings>\n', indent);
        end

        function writeQuantBlock(xmlHandle, quantSettings, thSettings, indent)
            fprintf(xmlHandle, '%s<QuantSettings', indent);
            fprintf(xmlHandle, ' DoClouds="%s"', TrueSpotXML.bool2str(~quantSettings.quantNoClouds));
            fprintf(xmlHandle, ' DoRefilter="%s"', TrueSpotXML.bool2str(~quantSettings.quantNoRefilter));
            fprintf(xmlHandle, ' CellZero="%s"', TrueSpotXML.bool2str(quantSettings.quantCellZero));

            if ~isnan(quantSettings.zAdj)
                fprintf(xmlHandle, ' ZtoXY="%f"', quantSettings.zAdj);
            end
            if quantSettings.fitRadXY > 0
                fprintf(xmlHandle, ' FitRadXY="%d"', quantSettings.fitRadXY);
            end
            if quantSettings.fitRadZ > 0
                fprintf(xmlHandle, ' FitRadZ="%d"', quantSettings.fitRadZ);
            end

            if thSettings.manualTh > 0
                fprintf(xmlHandle, ' ManualThreshold="%d"', thSettings.manualTh);
            end

            fprintf(xmlHandle, '/>\n');
        end

        function writeImageChannelBlock(xmlHandle, channelSettings, indent)
            fprintf(xmlHandle, '%s<ImageChannel ChannelNumber="%d">\n', indent, channelSettings.channelIndex);
            if ~isempty(channelSettings.metadata)
                TrueSpotXML.writeChannelMetaBlock(xmlHandle, channelSettings.metadata, [indent '\t']);
            end
            if ~isempty(channelSettings.spotCountSettings)
                TrueSpotXML.writeSpotsBlock(xmlHandle, ...
                    channelSettings.spotCountSettings, channelSettings.options, channelSettings.thresholdSettings, [indent '\t']);
            end
            if ~isempty(channelSettings.spotCountSettings)
                TrueSpotXML.writeQuantBlock(xmlHandle, ...
                    channelSettings.spotCountSettings, channelSettings.thresholdSettings, [indent '\t']);
            end
            fprintf(xmlHandle, '%s</ImageChannel>\n', indent);
        end

        function writeBatchBlock(xmlHandle, batchSettings, channelSettingsList, indent)
            if ~isempty(batchSettings.metadata) & ~isempty(batchSettings.metadata.batchName)
                fprintf(xmlHandle, '%s<ImageBatch Name="%s">\n', indent, batchSettings.metadata.batchName);
            else
                fprintf(xmlHandle, '%s<ImageBatch>\n', indent);
            end

            if ~isempty(batchSettings.metadata)
                TrueSpotXML.writeBatchMetaBlock(xmlHandle, batchSettings.metadata, [indent '\t']);
            end
            if ~isempty(batchSettings.paths)
                TrueSpotXML.writePathsBlock(xmlHandle, batchSettings.paths, [indent '\t']);
            end
            if ~isempty(batchSettings.cellsegSettings)
                TrueSpotXML.writeCellSegBlock(xmlHandle, batchSettings.cellsegSettings, [indent '\t']);
            end
            if ~isempty(batchSettings.channelInfo)
                TrueSpotXML.writeChannelInfoBlock(xmlHandle, batchSettings.channelInfo, channelSettingsList, [indent '\t']);
            end

            fprintf(xmlHandle, '%s</ImageBatch>\n', indent);
        end

        function writeSettingsXML(xmlpath, batchSettings, channels, commonMeta)
            if nargin < 4; commonMeta = []; end

            fh = fopen(xmlpath, 'w');

            fprintf(fh, '<?xml version="1.0" encoding="UTF-8"?>\n');
            if isempty(commonMeta)
                fprintf(fh, '<ImageSet>\n');
            else
                fprintf(fh, '<ImageSet Name="%s">\n', commonMeta.batchName);
            end
            
            if ~isempty(commonMeta)
                fprintf(fh, '\t<CommonMeta>\n');
                if ~isempty(commonMeta.species)
                    fprintf(fh, '\t\t<Species>"%s"</Species>\n', commonMeta.species);
                end
                if ~isempty(commonMeta.cellType)
                    fprintf(fh, '\t\t<CellType>"%s"</CellType>\n', commonMeta.cellType);
                end
                if ~isempty(commonMeta.voxelSize)
                    fprintf(fh, '\t\t<VoxelDimsNano X="%d" Y="%d" Z="%d"/>\n',... 
                        commonMeta.voxelSize.x, commonMeta.voxelSize.y, commonMeta.voxelSize.z);
                end
                fprintf(fh, '\t</CommonMeta>\n');
            end

            batchCount = size(batchSettings, 2);
            for b = 1:batchCount
                if iscell(batchSettings)
                    myBatch = batchSettings{b};
                else
                    myBatch = batchSettings(b);
                end
                
                chList = channels{b};

                TrueSpotXML.writeBatchBlock(fh, myBatch, chList, '\t');
            end

            fprintf(fh, '</ImageSet>\n');
            fclose(fh);
        end

    end

end