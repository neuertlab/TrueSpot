%%
classdef CountTables

     %%
    methods (Static)

        %Target -> Concentration -> Strain -> Repl -> Time

        function data = loadLiNeuertTable(path)
            trgNames = {'CTT1' 'STL1'};
            concNames = {'C200mM' 'C400mM'};
            rawTable = readtable(path, 'Format','%d%d%d%d%d%d%d');
            data =  struct();

            for ch = 1:2
                tname = trgNames{ch};
                chstr = struct();

                chData = rawTable(rawTable{:,'CH'} == ch,:);
                for exp = 1:2
                    expStr = struct();
                    strainStr = struct();
                    expData = chData(chData{:,'EXP'} == exp,:);

                    for rep = 1:3
                        if (rep == 3) & (exp == 1); break; end
                        repStr = struct();
                        repData = expData(expData{:,'REP'} == rep,:);

                        alltp = sort(unique(repData{:, 'TIME'}))';
                        tpCount = size(alltp, 2);
                        for tpi = 1:tpCount
                            tpVal = alltp(tpi);
                            tpName = CountTables.generateTimepointStructName(tpVal, 'min');
                            timeStr = struct();

                            timeStr.time = tpVal;
                            timeData = repData(repData{:,'TIME'} == tpVal,:);
                            %timeStr.total = timeData{:,'TOTAL'}';
                            timeStr.nuc = timeData{:,'NUC'}';
                            timeStr.nucCalls = timeData{:,'NUC'}';
                            timeStr.cyto = timeData{:,'TOTAL'}' - timeStr.nuc;

                            repStr.(tpName) = timeStr;
                        end

                        rName = ['R' num2str(rep)];
                        strainStr.(rName) = repStr;
                    end

                    expStr.BY4741 = strainStr;
                    eName = concNames{exp};
                    chstr.(eName) = expStr;
                end

                data.(tname) = chstr;
            end

        end

        function tableData = loadTSDumpTable_v1_2(path, timeColFunc, timeUnit, thSuffix, thSuffLo, thSuffHi)
            if nargin < 4; thSuffix = 'MID'; end
            if nargin < 5; thSuffLo = 'LO'; end
            if nargin < 6; thSuffHi = 'HI'; end

            cStatFmtStr = [repmat('%d', 1, 3) '%f' ];
            fmtStr = [repmat('%s', 1, 4) repmat('%d', 1, 2) repmat(cStatFmtStr, 1, 4) ...
                repmat('%d', 1, 7*3)];
            rawTable = readtable(path, 'FileType', 'text', 'Delimiter', '\t', 'Format', fmtStr);

            rawTable = timeColFunc(rawTable, 'TIME', timeUnit);

            tableData = struct();
            allTrg = unique(rawTable{:, 'TARGET'})';
            trgCount = size(allTrg, 2);

            %ismember('query', rawTable.Properties.VariableNames);
            useMedian = ismember('NUC_MED_INT', rawTable.Properties.VariableNames);

            for t = 1:trgCount
                trgName = allTrg{t};
                trgStr = struct();

                filtTable = rawTable(strcmp(rawTable{:,'TARGET'}, trgName),:);
                alltp = sort(unique(filtTable{:, 'TIME'}))';
                tpCount = size(alltp, 2);

                for tpi = 1:tpCount
                    tpVal = alltp(tpi);
                    tpName = CountTables.generateTimepointStructName(tpVal, timeUnit);
                    timeStr = struct();

                    timeData = filtTable(filtTable{:,'TIME'} == tpVal,:);

                    timeStr.time = tpVal;
                    timeStr.thVals = timeData{:,['THVAL' '_' thSuffix]}';
                    timeStr.cyto = timeData{:,['EST_COUNT_CYTO' '_' thSuffix]}';
                    timeStr.nuc = timeData{:,['EST_COUNT_NUC' '_' thSuffix]}';
                    %timeStr.nucCalls = timeData{:,['SPOTS_NUC' '_' thSuffix]}';
                    timeStr.nascent = timeData{:,['EST_NASCENT_COUNT_NUC' '_' thSuffix]}';
                    %timeStr.cytoSignal = timeData{:,['SIGNAL_CYTO' '_' thSuffix]}';
                    %timeStr.nucSignal = timeData{:,['SIGNAL_NUC' '_' thSuffix]}';
                    timeStr.cellAreaPix = timeData{:,'CELLAREA_PIX'}';
                    timeStr.nucVolVox = timeData{:,'NUCVOL_VOX'}';
                    timeStr.voxelDimsString = timeData{:,'VOXDIMS'}';

                    timeStr.nuc3Stats = struct();
                    timeStr.nuc3Stats.volume = timeData{:,'NUCVOL_VOX'}';
                    timeStr.nuc3Stats.totalIntensity = timeData{:,'NUC_TOT_INT'}';
                    if useMedian
                        timeStr.nuc3Stats.medianIntensity = timeData{:,'NUC_MED_INT'}';
                        timeStr.nuc3Stats.meanIntensity = NaN;
                    else
                        timeStr.nuc3Stats.medianIntensity = NaN;
                        timeStr.nuc3Stats.meanIntensity = timeData{:,'NUC_MEAN_INT'}';
                    end
                    timeStr.nuc3Stats.stdev = timeData{:,'NUC_INT_STDEV'}';

                    timeStr.nuc2Stats = struct();
                    timeStr.nuc2Stats.area = timeData{:,'NUC_MAXPROJ_AREA'}';
                    timeStr.nuc2Stats.totalIntensity = timeData{:,'NUC_MAXPROJ_TOT_INT'}';
                    if useMedian
                        timeStr.nuc2Stats.medianIntensity = timeData{:,'NUC_MAXPROJ_MED_INT'}';
                        timeStr.nuc2Stats.meanIntensity = NaN;
                    else
                        timeStr.nuc2Stats.medianIntensity = NaN;
                        timeStr.nuc2Stats.meanIntensity = timeData{:,'NUC_MAXPROJ_MEAN_INT'}';
                    end
                    timeStr.nuc2Stats.stdev = timeData{:,'NUC_MAXPROJ_INT_STDEV'}';

                    timeStr.cyto3Stats = struct();
                    timeStr.cyto3Stats.volume = timeData{:,'CYTOVOL_VOX'}';
                    timeStr.cyto3Stats.totalIntensity = timeData{:,'CYTO_TOT_INT'}';
                    if useMedian
                        timeStr.cyto3Stats.medianIntensity = timeData{:,'CYTO_MED_INT'}';
                        timeStr.cyto3Stats.meanIntensity = NaN;
                    else
                        timeStr.cyto3Stats.medianIntensity = NaN;
                        timeStr.cyto3Stats.meanIntensity = timeData{:,'CYTO_MEAN_INT'}';
                    end
                    timeStr.cyto3Stats.stdev = timeData{:,'CYTO_INT_STDEV'}';

                    timeStr.cyto2Stats = struct();
                    timeStr.cyto2Stats.volume = timeData{:,'CYTO_MAXPROJ_AREA'}';
                    timeStr.cyto2Stats.totalIntensity = timeData{:,'CYTO_MAXPROJ_TOT_INT'}';
                    if useMedian
                        timeStr.cyto2Stats.medianIntensity = timeData{:,'CYTO_MAXPROJ_MED_INT'}';
                        timeStr.cyto2Stats.meanIntensity = NaN;
                    else
                        timeStr.cyto2Stats.medianIntensity = NaN;
                        timeStr.cyto2Stats.meanIntensity = timeData{:,'CYTO_MAXPROJ_MEAN_INT'}';
                    end
                    timeStr.cyto2Stats.stdev = timeData{:,'CYTO_MAXPROJ_INT_STDEV'}';

                    timeStr.loThreshold = struct();
                    timeStr.loThreshold.thVals = timeData{:,['THVAL' '_' thSuffLo]}';
                    timeStr.loThreshold.cyto = timeData{:,['EST_COUNT_CYTO' '_' thSuffLo]}';
                    timeStr.loThreshold.nuc = timeData{:,['EST_COUNT_NUC' '_' thSuffLo]}';
                    %timeStr.loThreshold.nucCalls = timeData{:,['SPOTS_NUC' '_' thSuffLo]}';
                    timeStr.loThreshold.nascent = timeData{:,['EST_NASCENT_COUNT_NUC' '_' thSuffLo]}';
                    %timeStr.loThreshold.cytoSignal = timeData{:,['SIGNAL_CYTO' '_' thSuffLo]}';
                    %timeStr.loThreshold.nucSignal = timeData{:,['SIGNAL_NUC' '_' thSuffLo]}';

                    timeStr.hiThreshold = struct();
                    timeStr.hiThreshold.thVals = timeData{:,['THVAL' '_' thSuffHi]}';
                    timeStr.hiThreshold.cyto = timeData{:,['EST_COUNT_CYTO' '_' thSuffHi]}';
                    timeStr.hiThreshold.nuc = timeData{:,['EST_COUNT_NUC' '_' thSuffHi]}';
                    %timeStr.hiThreshold.nucCalls = timeData{:,['SPOTS_NUC' '_' thSuffHi]}';
                    timeStr.hiThreshold.nascent = timeData{:,['EST_NASCENT_COUNT_NUC' '_' thSuffHi]}';
                    %timeStr.hiThreshold.cytoSignal = timeData{:,['SIGNAL_CYTO' '_' thSuffHi]}';
                    %timeStr.hiThreshold.nucSignal = timeData{:,['SIGNAL_NUC' '_' thSuffHi]}';

                    trgStr.(tpName) = timeStr;
                end
                tableData.(trgName) = trgStr;
            end

        end

        function tableData = loadTSDumpTable(path, timeColFunc, timeUnit)
            %Target -> Time
            fmtStr = [repmat('%s', 1, 4) repmat('%d', 1, 13) repmat('%f', 1, 15)];
            rawTable = readtable(path, 'FileType', 'text', 'Delimiter', '\t', 'Format', fmtStr);
            rawTable(:, 'SRCIMGNAME') = rawTable(:, 'x_SRCIMGNAME');
            rawTable(:, 'x_SRCIMGNAME') = [];
            rawTable = timeColFunc(rawTable, 'TIME', timeUnit);

            tableData = struct();
            allTrg = unique(rawTable{:, 'TARGET'})';
            trgCount = size(allTrg, 2);

            for t = 1:trgCount
                trgName = allTrg{t};
                trgStr = struct();

                filtTable = rawTable(strcmp(rawTable{:,'TARGET'}, trgName),:);
                alltp = sort(unique(filtTable{:, 'TIME'}))';
                tpCount = size(alltp, 2);

                for tpi = 1:tpCount
                    tpVal = alltp(tpi);
                    tpName = CountTables.generateTimepointStructName(tpVal, timeUnit);
                    timeStr = struct();

                    timeData = filtTable(filtTable{:,'TIME'} == tpVal,:);

                    timeStr.time = tpVal;
                    timeStr.cyto = timeData{:,'EST_COUNT_CYTO'}';
                    timeStr.nuc = timeData{:,'EST_COUNT_NUC'}';
                    timeStr.nucCalls = timeData{:,'SPOTS_NUC'}';
                    timeStr.nascent = timeData{:,'EST_NASCENT_COUNT_NUC'}';
                    timeStr.cytoSignal = timeData{:,'SIGNAL_CYTO'}';
                    timeStr.nucSignal = timeData{:,'SIGNAL_NUC'}';
                    timeStr.cellAreaPix = timeData{:,'CELLAREA_PIX'}';
                    timeStr.nucVolVox = timeData{:,'NUCVOL_VOX'}';
                    timeStr.voxelDimsString = timeData{:,'VOXDIMS'}';

                    trgStr.(tpName) = timeStr;
                end
                tableData.(trgName) = trgStr;
            end
        end

        function tpStructName = generateTimepointStructName(tpValue, unitName)
            nstr = num2str(tpValue);
            nstr = replace(nstr, '.', 'p');
            tpStructName = ['T' nstr unitName];
        end

        function stringField = getTableFieldAsString(rawField)
            if iscell(rawField)
                stringField = strjoin(rawField,'');
            else
                stringField = string(rawField);
            end
            
        end

    end

end