%
%%
classdef VisInteractive

    %Call view modes:
            %   1: Results only
            %   2: Ref select
            %   3: Results versus other results
            %   4: Results versus reference

    %GUI Controls:
            %   space (or any key) - Ready listener
            %   s - Save (ref set)
            %   x - Exit


    %%
    properties
        spotsRun;

        imageCh; %Cache for loaded image channel.
        imageChFilt; %Cache for filtered channel.
        figNo = 615;
        figHandle;

        crosshairColor = [1.0 1.0 0.0];
        callViewMode = 1;
        currentSlice = 1;
        currentThresh = 0;
        otherThresh = 100;

        %Renderers/data containers
        visCommon = VisCommon;
        scVis = SpotCallVisualization;
        csVis = CellsegDrawer;
        qnVis = QuantVisualization;

        %View toggles
        useFilt = false; %Toggle between filtered and raw image channel
        maxProj = true; %Toggle between single slice view and max projection
        globalContrast = true; %(Slice view only) Toggle between contrast scale from full image or just that slice

        %Layer visibility
        callMode = 0; %Enum for mode of SpotCallVisualization
        imageLayerOn = true;
        cellLayerOn = false;
        nucLayerOn = false;
        spotCircleLayerOn = false;
        quantFitLayerOn = false;
        cloudLayerOn = false;

    end

    %%
    methods

        %% ================== Init ===============================

        %
        function obj = initCommon(obj, spotsrun, quantPath)
            if nargin < 3; quantPath = []; end %Try to guess from spotsrun

            obj.visCommon = obj.visCommon.initializeMe(spotsrun.dims.idims_sample);

            %Load whatever can be loaded (TIF, call table, quant, cellseg)
            fprintf('[VisInteractive.initCommon] Loading cell segmentation data...\n');
            obj.csVis = obj.csVis.initializeMe();
            obj.csVis.useMaxProj = obj.maxProj;
            csPath = spotsrun.paths.cellseg_path;
            if ~isempty(csPath)
                if isfile(csPath)
                    obj.csVis.cell_mask = CellSeg.openCellMask(csPath);
                    obj.csVis.nuc_mask = CellSeg.openNucMask(csPath);
                else
                    fprintf('[VisInteractive.initCommon] WARNING: Cell segmentation file "%s" not found. Cellseg will not be visible.\n', csPath);
                end
            else
                fprintf('[VisInteractive.initCommon] WARNING: Cell segmentation mask not provided. Cellseg will not be visible.\n');
            end

            %Call table
            fprintf('[VisInteractive.initCommon] Loading spot call data...\n');
            obj.scVis = obj.scVis.initializeMe();
            obj.scVis.visCommon = obj.visCommon;
            obj.scVis.zMode = 2;
            if obj.maxProj; obj.scVis.zMode = 1; end
            [spotsrun, obj.scVis.callTable] = spotsrun.loadCallTable();
            if isempty(obj.scVis.callTable)
                fprintf('[VisInteractive.initCommon] WARNING: Call table from "%s" could not be loaded!\n', [spotsrun.getFullOutStem() '_callTable.mat']);
            end

            %Quant results
            fprintf('[VisInteractive.initCommon] Loading quant/fit data...\n');
            obj.qnVis.maxProj = obj.maxProj;
            if isempty(quantPath)
                quantPath = [spotsrun.getFullOutStem() '_quantData.mat'];
            end

            if isfile(quantPath)
                load(quantPath, 'quant_results');

                if ~isfield(quant_results, 'cell_rna_data')
                    quant_results = RNAQuant.readResultsSavePackage(quant_results);
                end

                obj.qnVis.cells = quant_results.cell_rna_data;
            else
                fprintf('[VisInteractive.initCommon] WARNING: Quant data file "%s" not found. Quant results will not be visible.\n', quantPath);
            end

            %Source image
            fprintf('[VisInteractive.initCommon] Loading source image...\n');
            imgPath = spotsrun.paths.img_path;
            if ~isempty(imgPath)
                if isfile(imgPath)
                    [channels, ~] = LoadTif(imgPath, spotsrun.channels.total_ch, [spotsrun.channels.rna_ch], 1);
                    obj.imageCh = channels{spotsrun.channels.rna_ch, 1};
                    clear channels;

                    %Apply filter
                    [obj.imageChFilt] = RNA_Threshold_SpotDetector.run_spot_detection_pre(...
                        obj.imageCh, spotsrun.paths.out_dir, true, spotsrun.options.dtune_gaussrad, false);

                else
                    fprintf('[VisInteractive.initCommon] WARNING: Source image file "%s" not found. Cannot load source image.\n', imgPath);
                end
            else
                fprintf('[VisInteractive.initCommon] WARNING: Image path was not saved! Cannot load source image.\n');
            end

            %
            obj.currentThresh = spotsrun.intensity_threshold;
            obj.currentSlice = round(spotsrun.dims.idims_sample.z ./ 2);
            obj.spotsRun = spotsrun;
        end

        %
        function obj = initCallView(obj, spotsrun, quantPath)
            if nargin < 3; quantPath = []; end
            obj = obj.initCommon(obj, spotsrun, quantPath);
            obj.callViewMode = 1;
        end

        %
        function obj = initRefSelectView(obj, spotsrun, quantPath)
            if nargin < 3; quantPath = []; end
            obj = obj.initCommon(obj, spotsrun, quantPath);
            obj = obj.loadRefSet();
            obj.callViewMode = 2;
        end

        %
        function obj = initCallCompareView(obj, spotsrun, otherCallTable, quantPath)
            if nargin < 4; quantPath = []; end
            obj = obj.initCommon(obj, spotsrun, quantPath);
            obj.scVis.callTableCompare = otherCallTable;
            obj.callViewMode = 3;
        end

        %
        function obj = initRefCompareView(obj, spotsrun, quantPath)
            if nargin < 3; quantPath = []; end
            obj = obj.initCommon(obj, spotsrun, quantPath);
            obj = obj.loadRefSet();
            obj.callViewMode = 4;
        end

        %% ================== Save/Load ===============================
        %
        function path = getRefSavePath(obj)
            path = [obj.spotsRun.getFullOutStem() '_refset.mat'];
        end

        %
        function obj = loadRefSet(obj)
            refpath = obj.getRefSavePath();
            if isfile(refpath)
                load(refpath, 'refset');
                obj.scVis.referenceTable = refset;
            else
                %Start empty.
                obj.scVis.referenceTable = [];
            end
        end

        %
        function obj = saveRefSet(obj)
            refpath = obj.getRefSavePath();
            refset = obj.scVis.referenceTable;
            version = 1;
            timestamp = datetime();
            save(refpath, 'version', 'timestamp', 'refset');
        end

        %% ================== Rendering ===============================

        %
        function [obj, irender] = applyWorkAreaMask(obj, irender)
            %Generate xy mask to darken pixels outside selected region
            okayRegion = obj.scVis.workRegion;
            if ~isempty(okayRegion)
                Y = size(obj.imageCh,1);
                X = size(obj.imageCh,2);
                immul_mask = double(ones(Y,X));
                dval = 0.25;

                zout = false;
                if ~obj.maxProj
                    %Check if z in range.
                    if obj.currentSlice < okayRegion.z_min | obj.currentSlice > okayRegion.z_max
                        immul_mask(:,:) = dval;
                        zout = true;
                    end
                end

                if ~zout
                    immul_mask(:,1:okayRegion.x_min) = dval;
                    immul_mask(:,okayRegion.x_max:X) = dval;
                    immul_mask(1:okayRegion.y_min,:) = dval;
                    immul_mask(okayRegion.y_max:Y,:) = dval;
                end
                clear zout

                if ndims(irender) > 2
                    N = size(irender, 3);
                    for n = 1:N
                        irender(:,:,n) = immultiply(irender(:,:,n), immul_mask);
                    end
                else
                    irender = immultiply(irender, immul_mask);
                end

                clear immul_mask;
            end
        end

        %
        function [obj, irender] = render2dBase(obj)
            %Generate z projection.
            if obj.maxProj
                if obj.useFilt
                    irender = double(max(obj.imageChFilt,[],3));
                else
                    irender = double(max(obj.imageCh,[],3));
                end
            else
                if obj.useFilt
                    irender = double(obj.imageChFilt(:,:,obj.currentSlice));
                else
                    irender = double(obj.imageCh(:,:,obj.currentSlice));
                end
            end
            %mm = max(irender, [], 'all', 'omitnan');
            %irender = irender ./ mm;
            %imshow(irender, [obj.imshowMin obj.imshowMax]);
        end

        %
        function obj = updateRender(obj)
            if isempty(obj.figHandle)
                obj.figHandle = figure(obj.figNo);
                obj = updateImshowScaling(obj);
            else
                figure(obj.figHandle);
            end

            Y = size(obj.imageCh,1);
            X = size(obj.imageCh,2);

            clf;

            if obj.imageLayerOn
                [obj, irender] = obj.render2dBase();
            else
                irender = zeros(Y,X);
            end

            %Outputs a uint8
            if obj.cellLayerOn
                irender = obj.csVis.applyCellMask(irender, obj.currentSlice);
            else
                irender = uint8(zeros(Y, X, 3));
            end

            if obj.nucLayerOn
                irender = obj.csVis.applyNucMask(irender, obj.currentSlice);
            end

            %Quant?
            %RGB8 in, 0-1 double out
            if obj.quantFitLayerOn
                irender = obj.qnVis.applyToImage(irender, obj.currentSlice);
            else
                irender = double(irender) ./ 255.0;
            end

            %Apply workarea mask
            [obj, irender] = obj.applyWorkAreaMask(irender);

            %Rescale to uint8
            rgb = uint8(zeros(Y, X, 3));
            rgb(:,:,1) = uint8(irender(:,:,1) .* 255.0);
            rgb(:,:,2) = uint8(irender(:,:,2) .* 255.0);
            rgb(:,:,3) = uint8(irender(:,:,3) .* 255.0);

            clear irender
            imshow(rgb);

            %Spots
            if obj.spotCircleLayerOn
                if obj.callMode == 1
                    [obj.figHandle, ~] = obj.scVis.drawResultsBasic(obj.figHandle, obj.currentThresh, obj.currentSlice);
                elseif obj.callMode == 2
                    [obj.figHandle, ~] = obj.scVis.drawReferenceSet(obj.figHandle, obj.currentSlice);
                elseif obj.callMode == 3
                    [obj.figHandle, ~] = obj.scVis.drawResultsTTCompare(obj.figHandle, obj.currentThresh, obj.otherThresh, obj.currentSlice);
                elseif obj.callMode == 4
                    [obj.figHandle, ~] = obj.scVis.drawResultsTRCompare(obj.figHandle, obj.currentThresh, obj.currentSlice);
                end
            end

        end

    end

    %%
    methods(Static)
    end

end