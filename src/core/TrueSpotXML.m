%
%%
classdef TrueSpotXML

     %%
    methods (Static)

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
                        cellsegSettings.xtrim = TrueSpotXML.getElementTextAsNumber(sChild, 4);
                    elseif strcmp(sChildName, 'YTrim')
                        cellsegSettings.ytrim = TrueSpotXML.getElementTextAsNumber(sChild, 4);
                    elseif strcmp(sChildName, 'NucZRange')
                        cellsegSettings.nzrange = TrueSpotXML.getElementTextAsNumber(sChild, 3);
                    elseif strcmp(sChildName, 'NucThSample')
                        cellsegSettings.nthsmpl = TrueSpotXML.getElementTextAsNumber(sChild, 200);
                    elseif strcmp(sChildName, 'NucCutoff')
                        cellsegSettings.ncutoff = TrueSpotXML.getElementTextAsNumber(sChild, 0.05);
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
                spotSettings = TrueSpotChannelSettings.genThresholdSettingsStruct();
            end
            if isempty(optionsStruct)
                spotSettings = TrueSpotChannelSettings.genOptionsStruct();
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
        end

    end

end