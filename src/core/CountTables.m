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
                            tpName = ['T' num2str(tpVal) 'min'];
                            timeStr = struct();

                            timeStr.time = tpVal;
                            timeData = repData(repData{:,'TIME'} == tpVal,:);
                            %timeStr.total = timeData{:,'TOTAL'}';
                            timeStr.nuc = timeData{:,'NUC'}';
                            timeStr.cyto = timeData{:,'TOTAL'}' - timeStr.nuc;

                            repStr.(tpName) = timeStr;
                        end

                        rName = ['R' num2str(rep)];
                        strainStr.(rName) = repStr;
                    end

                    expStr.BY4147 = strainStr;
                    eName = concNames{exp};
                    chstr.(eName) = expStr;
                end

                data.(tname) = chstr;
            end

        end

        function tableData = loadTSDumpTable(path, timeColFunc)
            %Target -> Time
            fmtStr = [repmat('%s', 1, 4) repmat('%d', 1, 13)];
            rawTable = readtable(path, 'Format',fmtStr);
            rawTable = timeColFunc(rawTable, 'TIME');

            tableData = struct();
            allTrg = unique(rawTable{:, 'TARGET'});
            trgCount = size(allTrg, 2);

            for t = 1:trgCount
                trgName = allTrg{t};
                trgStr = struct();

                filtTable = rawTable(strcmp(rawTable{:,'TARGET'}, trgName),:);
                alltp = sort(unique(filtTable{:, 'TIME'}))';
                tpCount = size(alltp, 2);

                for tpi = 1:tpCount
                    tpVal = alltp(tpi);
                    tpName = ['T' num2str(tpVal) 'min'];
                    timeStr = struct();

                    timeData = filtTable(filtTable{:,'TIME'} == tpVal,:);

                    timeStr.time = tpVal;
                    timeStr.cyto = timeData{:,'EST_COUNT_CYTO'}';
                    timeStr.nuc = timeData{:,'EST_COUNT_NUC'}';
                    timeStr.nascent = timeData{:,'EST_NASCENT_COUNT_NUC'}';
                    timeStr.cytoSignal = timeData{:,'SIGNAL_CYTO'}';
                    timeStr.nucSignal = timeData{:,'SIGNAL_NUC'}';
                    timeStr.cellAreaPix = timeData{:,'CELLAREA(PIX)'}';
                    timeStr.nucVolVox = timeData{:,'NUCVOL(VOX)'}';
                    timeStr.voxelDimsString = timeData{:,'VOXDIMS'}';

                    trgStr.(tpName) = timeStr;
                end
                tableData.(trgName) = trgStr;
            end
        end

    end

end