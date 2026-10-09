%
%%
classdef ResViewerFunctions

    methods(Static)

        %---------------------------- Structs ---------------------------

        %%
        function rvsmStruct = genResViewerSampleMetaStuct()
            rvsmStruct = struct('structVer', 1);
            
            rvsmStruct.channelIndex = 0;
            rvsmStruct.fluorName = [];
            rvsmStruct.targetName = [];
        end

        %%
        function rviStruct = genResViewerImageStruct()
            rviStruct = struct('structVer', 1);
            
            rviStruct.imagePath = []; %TIF path
            rviStruct.cellSegPath = [];
            % rviStruct.totalCh = 0;
            % rviStruct.lightCh = 0;
            % rviStruct.nucCh = 0;

            rviStruct.sampleRunPaths = []; %Paths to spotsrun for each channel
            rviStruct.quantPaths = [];
        end

        %%
        function rvbStruct = genResViewerBatchStruct()
            rvbStruct = struct('structVer', 1);
            
            rvbStruct.batchName = [];
            rvbStruct.batchDir = [];
            rvbStruct.imageInfo = [];

            rvbStruct.totalCh = 0;
            rvbStruct.lightCh = 0;
            rvbStruct.nucCh = 0;
            rvbStruct.sampleInfo = [];
        end

        %---------------------------- Import ---------------------------

        %%
        function rvbStruct = importBatchDir(dirPath)
            %Scan recursively for spotsrun files
            rvbStruct = [];
            runFileList = FileUtils.scanForFilesEndingWith(dirPath, '_rnaspotsrun.mat', true);
            if isempty(runFileList)
                return;
            end

            %Determine which groups of spotsrun files belong to each image
            %Load first spotsrun to get some channel info
            %Channel info is assumed to be the same across batch
            srTestPath = runFileList{1};
            spotsRun = RNASpotsRun.loadFrom(srTestPath, true);
            rvbStruct = ResViewerFunctions.genResViewerBatchStruct();
            rvbStruct.batchDir = dirPath;
            [~, rvbStruct.batchName, ~] = fileparts(dirPath);
            rvbStruct.totalCh = spotsRun.channels.total_ch;
            rvbStruct.lightCh = spotsRun.channels.light_ch;
            rvbStruct.sampleInfo = cell(1, spotsRun.channels.total_ch);
            lsize = size(runFileList, 2);
            clear srTestPath spotsRun

            iPathLookup = dictionary('dummy', struct());
            for i = 1:lsize
                srPath = runFileList{i};
                spotsRun = RNASpotsRun.loadFrom(srPath, true);
                ogImgPath = spotsRun.paths.img_path;
                sampleCh = spotsRun.channels.rna_ch;

                if iPathLookup.isKey(ogImgPath)
                    rviStruct = iPathLookup.lookup(ogImgPath);
                else
                    rviStruct = ResViewerFunctions.genResViewerImageStruct();
                    rviStruct.imagePath = ogImgPath;
                    rviStruct.cellSegPath = spotsRun.paths.cellseg_path;

                    rviStruct.sampleRunPaths = cell(1, rvbStruct.totalCh);
                    rviStruct.quantPaths = cell(1, rvbStruct.totalCh);

                    %Check cellseg path and grab nuc channel
                    if ~isfile(rviStruct.cellSegPath)
                        [~, csName, csExt] = fileparts(rviStruct.cellSegPath);
                        csName = [csName csExt];
                        clear csExt

                        %First check the directory above the spots dir
                        [parentDir, ~, ~] = fileparts(srPath);
                        [parentDir, ~, ~] = fileparts(parentDir);
                        flist = FileUtils.scanForFilesWithName(parentDir, csName, false);
                        clear parentDir
                        if ~isempty(flist)
                            rviStruct.cellSegPath = flist{1};
                        else
                            %Then, scan entire batch directory
                            flist = FileUtils.scanForFilesWithName(dirPath, csName, true);
                            if ~isempty(flist)
                                rviStruct.cellSegPath = flist{1};
                            end
                        end
                        clear csName flist
                    end

                    if isfile(rviStruct.cellSegPath) & (rvbStruct.nucCh < 1)
                        load(rviStruct.cellSegPath, 'runMeta');
                        rvbStruct.nucCh = runMeta.srcImageChNuc;
                        clear runMeta
                    end
                end

                rviStruct.sampleRunPaths{sampleCh} = srPath;
                %Look for quant path in the same directory
                qpath = [];
                [srDir, ~, ~] = fileparts(srPath);
                dContents = dir(srDir);
                itemCount = size(dContents, 1);
                for j = 1:itemCount
                    dChild = dContents(j);
                    if ~dChild.isdir
                        if endsWith(dChild.name, '_quantData.mat')
                            qpath = [srDir filesep dChild.name];
                            break;
                        end
                    end
                end
                if ~isempty(qpath)
                    rviStruct.quantPaths{sampleCh} = qpath;
                end
                clear srDir qpath j itemCount dChild

                rvsmStruct = rvbStruct.sampleInfo{sampleCh};
                if isempty(rvsmStruct)
                    rvsmStruct = ResViewerFunctions.genResViewerSampleMetaStuct();
                    rvsmStruct.channelIndex = sampleCh;
                    rvsmStruct.fluorName = spotsRun.meta.type_probe;
                    rvsmStruct.targetName = spotsRun.meta.type_target;
                    rvbStruct.sampleInfo{sampleCh} = rvsmStruct;
                end

                iPathLookup = iPathLookup.insert(ogImgPath, rviStruct); %Overwrite
                clear srPath spotsRun ogImgPath sampleCh rvsmStruct rviStruct
            end

            %Copy image info structs to batch struct
            iPathLookup = iPathLookup.remove('dummy');
            rvbStruct.imageInfo = iPathLookup.values('cell');
            rvbStruct.imageInfo = rvbStruct.imageInfo';
        end

        %%
        function rvbStruct = updateImagePaths(rvbStruct, newImageDir)
            if isempty(rvbStruct)
                return;
            end

            rviCount = size(rvbStruct.imageInfo, 2);
            for i = 1:rviCount
                rviStruct = rvbStruct.imageInfo{i};
                [~,n,e] = fileparts(rviStruct.imagePath);
                ifileName = [n e];
                searchRes = FileUtils.scanForFilesWithName(newImageDir, ifileName, true);
                if ~isempty(searchRes)
                    rviStruct.imagePath = searchRes{1};
                end
                rvbStruct.imageInfo{i} = rviStruct;
                clear n e ifileName rviStruct searchRes
            end
        end

        %---------------------------- Save/Load ---------------------------

        %%
        function [rvbStruct, rvState] = loadResViewerBatch(filepath)
            load(filepath, 'rvbStruct', 'rvState');
        end

        %%
        function saveResViewerBatch(rvbStruct, rvState, filepath)
            save(filepath, 'rvbStruct', 'rvState');
        end

    end
end