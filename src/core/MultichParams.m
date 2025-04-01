%
%%
classdef MultichParams

    methods(Static)
    
        %%
        function sampleInfo = genEmptySampleStruct()
            sampleInfo = struct();
            sampleInfo.name = [];
            sampleInfo.chIdx = 0;
            sampleInfo.filters = [];
            sampleInfo.thLow = 0;
            sampleInfo.thMid = 0;
            sampleInfo.thHigh = 0;
        end

        %%
        function channelInfos = readChannelDefXML(xmlpath)
            %channelInfos = [];
            rootnode = xmlread(xmlpath);
            smplElements = getElementsByTagName(rootnode,'SampleChannel');
            smplCount = smplElements.getLength;

            channelInfos(smplCount) = MultichParams.genEmptySampleStruct();
            for i = (smplCount-1):-1:1
                channelInfos(i) = MultichParams.genEmptySampleStruct();
            end

            for i = 1:smplCount
                sInfo = channelInfos(i);

                ji = i-1;
                smplElement = item(smplElements, ji);

                %Get base attributes
                sInfo.name = char(getAttribute(smplElement, 'Name'));
                attrStr = char(getAttribute(smplElement, 'ChannelIndex'));
                if ~isempty(attrStr)
                    sInfo.chIdx = Force2Num(attrStr);
                end

                filtersPre = cell(1, 16);
                filtersUsed = 0;

                sChild = getFirstChild(smplElement);
                while ~isempty(sChild)
                    if sChild.getNodeType == sChild.ELEMENT_NODE
                        sChildName = char(sChild.getTagName);
                        if strcmp(sChildName, 'QueryParams')
                            gChild = getFirstChild(sChild);
                            while ~isempty(gChild)
                                if gChild.getNodeType == sChild.ELEMENT_NODE
                                    gChildName = char(gChild.getTagName);
                                    if strcmp(gChildName, 'QueryParam')
                                        myFilter = struct();
                                        myFilter.key = char(getAttribute(gChild, 'Key'));
                                        myFilter.value = char(getAttribute(gChild, 'Value'));
                                        filtersUsed = filtersUsed + 1;
                                        filtersPre{filtersUsed} = myFilter;
                                    end
                                end
                                gChild = getNextSibling(gChild);
                            end
                        elseif strcmp(sChildName, 'CountThresholds')
                            attrStr = char(getAttribute(sChild, 'ThMid'));
                            if ~isempty(attrStr)
                                sInfo.thMid = Force2Num(attrStr);
                            end
                            attrStr = char(getAttribute(sChild, 'ThLow'));
                            if ~isempty(attrStr)
                                sInfo.thLow = Force2Num(attrStr);
                            end
                            attrStr = char(getAttribute(sChild, 'ThHigh'));
                            if ~isempty(attrStr)
                                sInfo.thHigh = Force2Num(attrStr);
                            end
                        end
                    end
                    sChild = getNextSibling(sChild);
                end

                %Port over any found filters
                if filtersUsed > 0
                    %sInfo.filters = cell(1, filtersUsed);
                    sInfo.filters = filtersPre(1:filtersUsed); 
                end

                channelInfos(i) = sInfo;
            end

        end

        %%
        function sampleInfo = matchRunToChannel(spotsrun, channelInfos)
            sampleInfo = [];
            if isempty(spotsrun); return; end
            if isempty(channelInfos); return; end

            sampleCount = size(channelInfos, 2);
            for i = 1:sampleCount
                match = MultichParams.passesFilters(spotsrun, channelInfos(i));
                if match
                    sampleInfo = channelInfos(i);
                    return;
                end
            end
        end

        %%
        function boolres = passesFilters(spotsrun, sampleInfo)
            boolres = false;
            if isempty(spotsrun); return; end

            boolres = true;
            if isempty(sampleInfo); return; end
            if isempty(sampleInfo.filters); return; end

            filterCount = size(sampleInfo.filters, 2);
            for i = 1:filterCount
                myFilter = sampleInfo.filters{i};
                pass = MultichParams.passesFilter(spotsrun, myFilter);
                if ~pass
                    boolres = false;
                    return;
                end
            end
        end

        %%
        function boolres = passesFilter(spotsrun, filter)
            boolres = false;
            if isempty(spotsrun); return; end

            boolres = true;
            if isempty(filter); return; end

            if strcmp(filter.key, 'ChannelNo')
                val = Force2Num(filter.value);
                boolres = (spotsrun.channels.rna_ch == val);
            elseif strcmp(filter.key, 'INameContains')
                boolres = contains(spotsrun.img_name, filter.value);
            elseif strcmp(filter.key, 'ProbeName')
                boolres = contains(spotsrun.meta.type_probe, filter.value);
            elseif strcmp(filter.key, 'TargetName')
                boolres = contains(spotsrun.meta.type_target, filter.value);
            elseif strcmp(filter.key, 'CellType')
                boolres = contains(spotsrun.meta.type_cell, filter.value);
            elseif strcmp(filter.key, 'SpeciesName')
                boolres = contains(spotsrun.meta.type_species, filter.value);
            end
        end


    end

end