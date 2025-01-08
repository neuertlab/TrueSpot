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
            %   c - Clear all ref spots
            %   3 - 3D snap to local maxima
            %   z - Z only snap to local maxima
            %   R - Increase xy snap radius by 1
            %   r - Decrease xy snap radius by 1

            % --- Threshold/Z Slice Movement ---
            %   + - Up one z slice
            %   - - Down one z slice
            %   = - Up ten z slices
            %   _ - Down ten z slices
            %   > - Increase threshold by 1
            %   < - Decrease threshold by 1
            %   ] - Increase threshold by 10
            %   [ - Decrease threshold by 10

            % --- Toggles ---
            %   f - Filtered/raw image toggle
            %   m - Maximum projection/indiv slice toggle
            %   g - Global (z) contrast toggle
            %   U - When snapping, removed any spots that do not snap to a maximum

            % --- Layers ---
            %   C - Cells on/off
            %   N - Nuclei on/off
            %   S - Spot calls on/off
            %   F - Spot fits on/off

            % --- Mask ---
            %   M - Enter XY mask selection mode (use two clicks to make rectangle)
            %   Z - Toggle Z min/max control
            %   ^ - Z boundary up
            %   v - Z boundary down


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
        zSetToggle = false; %False = min, true = max
        remUnsnapped = false;

        %Layer visibility
        %callMode = 0; %Enum for mode of SpotCallVisualization
        imageLayerOn = true;
        cellLayerOn = false;
        nucLayerOn = false;
        spotCircleLayerOn = false;
        quantFitLayerOn = false;
        cloudLayerOn = false;

        %GUI state
        loopBreaker = false;

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
                    obj.csVis.cell_mask = (obj.csVis.cell_mask > 0);
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
                clear quant_results;
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

            if ~isempty(obj.qnVis.cells)
                fprintf('[VisInteractive.initCommon] Prerendering fit data...\n');
                [obj.qnVis, ~] = obj.qnVis.prerenderCells();
            end

            %
            obj.currentThresh = spotsrun.intensity_threshold;
            obj.currentSlice = round(spotsrun.dims.idims_sample.z ./ 2);
            obj.spotsRun = spotsrun;
        end

        %
        function obj = initCallView(obj, spotsrun, quantPath)
            if nargin < 3; quantPath = []; end
            obj = obj.initCommon(spotsrun, quantPath);
            obj.callViewMode = 1;
        end

        %
        function obj = initRefSelectView(obj, spotsrun, quantPath)
            if nargin < 3; quantPath = []; end
            obj = obj.initCommon(spotsrun, quantPath);
            obj = obj.loadRefSet();
            obj.callViewMode = 2;
        end

        %
        function obj = initCallCompareView(obj, spotsrun, otherCallTable, quantPath)
            if nargin < 4; quantPath = []; end
            obj = obj.initCommon(spotsrun, quantPath);
            obj.scVis.callTableCompare = otherCallTable;
            obj.callViewMode = 3;
        end

        %
        function obj = initRefCompareView(obj, spotsrun, quantPath)
            if nargin < 3; quantPath = []; end
            obj = obj.initCommon(spotsrun, quantPath);
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
                load(refpath, 'refset', 'version', 'workRegion');
                obj.scVis.referenceTable = refset;
                if version > 1
                    obj.scVis.workRegion = workRegion;
                end
            else
                %Start empty.
                obj.scVis = obj.scVis.setRefSetEmpty();
            end
        end

        %
        function obj = saveRefSet(obj)
            refpath = obj.getRefSavePath();
            refset = obj.scVis.referenceTable;
            version = 2;
            workRegion = obj.scVis.workRegion;
            timestamp = datetime();
            save(refpath, 'version', 'timestamp', 'refset', 'workRegion');
        end

        %% ================== Rendering ===============================

        %
        function [obj, irender] = applyWorkAreaMask(obj, irender)
            %Generate xy mask to darken pixels outside selected region
            okayRegion = obj.scVis.workRegion;
            if ~isempty(okayRegion)
                Y = size(obj.imageCh,1);
                X = size(obj.imageCh,2);
                Z = size(obj.imageCh,3);
                immul_mask = double(ones(Y,X));
                dval = 0.25;

                x_min = max(1, okayRegion.x_min);
                x_max = min(X, okayRegion.x_max);
                if x_max < 1; x_max = X; end
                y_min = max(1, okayRegion.y_min);
                y_max = min(Y, okayRegion.y_max);
                if y_max < 1; y_max = Y; end
                z_min = max(1, okayRegion.z_min);
                z_max = min(Z, okayRegion.z_max);
                if z_max < 1; z_max = Z; end

                zout = false;
                if ~obj.maxProj
                    %Check if z in range.
                    if obj.currentSlice < z_min | obj.currentSlice > z_max
                        immul_mask(:,:) = dval;
                        zout = true;
                    end
                end

                if ~zout
                    if (x_min > 1); immul_mask(:,1:x_min) = dval; end
                    if (x_max < X); immul_mask(:,x_max:X) = dval; end
                    if (y_min > 1); immul_mask(1:y_min,:) = dval; end
                    if (y_max < Y); immul_mask(y_max:Y,:) = dval; end
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
                cScaleMin = min(irender(:));
                cScaleMax = median(irender(:)) + round(10 * std(irender(:)));
            else
                if obj.useFilt
                    irender = double(obj.imageChFilt(:,:,obj.currentSlice));
                else
                    irender = double(obj.imageCh(:,:,obj.currentSlice));
                end

                if obj.globalContrast
                    if obj.useFilt
                        allvals = double(obj.imageChFilt(:));
                        cScaleMin = min(allvals);
                        cScaleMax = median(allvals) + round(10 * std(allvals));
                        clear allvals
                    else
                        allvals = double(obj.imageCh(:));
                        cScaleMin = min(allvals);
                        cScaleMax = median(allvals) + round(10 * std(allvals));
                        clear allvals
                    end
                else
                    cScaleMin = min(irender(:));
                    cScaleMax = median(irender(:)) + round(10 * std(irender(:)));
                end
            end

            %Rescale
            irender = irender - cScaleMin;
            rr = cScaleMax - cScaleMin;
            irender = irender ./ rr;

            %mm = max(irender, [], 'all', 'omitnan');
            %irender = irender ./ mm;
            %imshow(irender, [obj.imshowMin obj.imshowMax]);
        end

        %
        function obj = updateRender(obj)
            if isempty(obj.figHandle)
                obj.figHandle = figure(obj.figNo);
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

            % ----- DEBUG
            imshow(irender);

            %Outputs a uint8
            if obj.cellLayerOn
                irender = obj.csVis.applyCellMask(irender, obj.currentSlice);
            else
                %Rescale input!!
                %inmax = max(irender, [], 'all', 'omitnan');
                %irender = irender ./ inmax;
                irender = irender .* 255.0;
                irender = round(irender);
                irender(irender > 255) = 255;
                %clear inmax

                irender8 = uint8(zeros(Y, X, 3));
                irender8(:,:,1) = uint8(irender(:,:));
                irender8(:,:,2) = irender8(:,:,1);
                irender8(:,:,3) = irender8(:,:,1);
                irender = irender8;
                clear irender8;
            end

            % ----- DEBUG
            clf;
            imshow(irender);

            if obj.nucLayerOn
                irender = obj.csVis.applyNucMask(irender, obj.currentSlice);
            end

            % ----- DEBUG
            clf;
            imshow(irender);

            %Quant?
            %RGB8 in, 0-1 double out
            if obj.quantFitLayerOn
                irender = obj.qnVis.applyToImage(irender, obj.currentSlice);
            else
                irender = double(irender) ./ 255.0;
            end

            % ----- DEBUG
            figure(obj.figHandle);
            clf;
            imshow(irender);

            %Apply workarea mask
            [obj, irender] = obj.applyWorkAreaMask(irender);

            % ----- DEBUG
            clf;
            imshow(irender);

            %Rescale to uint8
            rgb = uint8(zeros(Y, X, 3));
            rgb(:,:,1) = uint8(irender(:,:,1) .* 255.0);
            rgb(:,:,2) = uint8(irender(:,:,2) .* 255.0);
            rgb(:,:,3) = uint8(irender(:,:,3) .* 255.0);

            clear irender
            imshow(rgb);

            %Spots
            if obj.spotCircleLayerOn
                if obj.callViewMode == 1
                    [obj.figHandle, ~] = obj.scVis.drawResultsBasic(obj.figHandle, obj.currentThresh, obj.currentSlice);
                elseif obj.callViewMode == 2
                    [obj.figHandle, ~] = obj.scVis.drawReferenceSet(obj.figHandle, obj.currentSlice);
                elseif obj.callViewMode == 3
                    [obj.figHandle, ~] = obj.scVis.drawResultsTTCompare(obj.figHandle, obj.currentThresh, obj.otherThresh, obj.currentSlice);
                elseif obj.callViewMode == 4
                    [obj.figHandle, ~] = obj.scVis.drawResultsTRCompare(obj.figHandle, obj.currentThresh, obj.currentSlice);
                end
            end

            if ~obj.maxProj
                title(['z = ' num2str(obj.currentSlice)]);
            else
                title('');
            end

        end

        %% ================== GUI ===============================

        %
        function obj = launchFigureGUI(obj)
            %TODO
            obj = obj.updateRender();

            while ~obj.loopBreaker
                w = waitforbuttonpress;
                if w
                    obj = obj.onReadyKey();
                end
            end
        end

        %
        function obj = onReadyKey(obj)
            %TODO
            [x,y,btn] = ginput_color(1, obj.crosshairColor);

            if btn == 1 %Mouse click
                %Loop until any key pressed that isn't mouse
                obj = obj.whileMouseListening(x,y,btn);
            elseif  btn == 3 %Right click
                obj = obj.whileMouseListening(x,y,btn);
            elseif btn == 'x' %120 'x' - exit
                close(obj.figHandle);
                obj.loopBreaker = true;
            elseif btn == 's' %115 's' - save
                if(obj.callViewMode == 2)
                    obj = obj.saveRefSet();
                    refpath = obj.getRefSavePath();
                    fprintf('Ref set saved to: %s\n', refpath);
                else
                    fprintf('GUI is not in reference edit mode!\n');
                end
            elseif btn == 'c' %'c' - clear
                obj.scVis = obj.scVis.setRefSetEmpty();
                obj = obj.updateRender();
            elseif btn == '3' %'3' - 3D snap
                [obj.scVis, snappedCount, removedCount] =...
                    obj.scVis.refModeSnap(obj.imageChFilt, false, obj.maxProj, obj.remUnsnapped, true);
                fprintf('3D local snap done.\n');
                fprintf('\tSpots Snapped: %d\n', snappedCount);
                fprintf('\tSpots Removed: %d\n', removedCount);
                obj = obj.updateRender();
            elseif btn == 'z' %'z' - z only snap
                [obj.scVis, snappedCount, removedCount] =...
                    obj.scVis.refModeSnap(obj.imageChFilt, true, obj.maxProj, obj.remUnsnapped, true);
                fprintf('Z only local snap done.\n');
                fprintf('\tSpots Snapped: %d\n', snappedCount);
                fprintf('\tSpots Removed: %d\n', removedCount);
            elseif btn == 'U' %'U' - remove unsnapped spots toggle
                obj.remUnsnapped = ~obj.remUnsnapped;
                fprintf('Remove unsnapped spots toggle: %d\n', obj.remUnsnapped);
            elseif btn == 'f' %toggle filtered/raw
                obj.useFilt = ~obj.useFilt;
                fprintf('Filter toggle set: %d\n', obj.useFilt);
                obj = obj.updateRender();
            elseif btn == 'm' %toggle max projection
                obj.maxProj = ~obj.maxProj;
                obj.csVis.useMaxProj = obj.maxProj;
                if obj.maxProj
                    obj.scVis.zMode = 1;
                else
                    obj.scVis.zMode = 2;
                end
                obj.qnVis.maxProj = obj.maxProj;
                fprintf('Max projection toggle set: %d\n', obj.maxProj);
                obj = obj.updateRender();
            elseif btn == 'g' %toggle global contrast scale
                obj.globalContrast = ~obj.globalContrast;
                fprintf('Global contrast toggle set: %d\n', obj.globalContrast);
                obj = obj.updateRender();
            elseif btn == 'C' %toggle cell layer
                obj.cellLayerOn = ~obj.cellLayerOn;
                if(obj.cellLayerOn)
                    fprintf('Cell boundary estimate view: ON\n');
                else
                    fprintf('Cell boundary estimate view: OFF\n');
                end
                obj = obj.updateRender();
            elseif btn == 'N' %toggle nuc layer
                obj.nucLayerOn = ~obj.nucLayerOn;
                if(obj.nucLayerOn)
                    fprintf('Nucleus boundary estimate view: ON\n');
                else
                    fprintf('Nucleus boundary estimate view: OFF\n');
                end
                obj = obj.updateRender();
            elseif btn == 'S' %toggle spot call layer
                obj.spotCircleLayerOn = ~obj.spotCircleLayerOn;
                if(obj.spotCircleLayerOn)
                    fprintf('Spot call view: ON\n');
                else
                    fprintf('Spot call view: OFF\n');
                end
                obj = obj.updateRender();
            elseif btn == 'F' %toggle signal fit (quant) layer
                obj.quantFitLayerOn = ~obj.quantFitLayerOn;
                if(obj.quantFitLayerOn)
                    fprintf('Signal fit (quant) view: ON\n');
                else
                    fprintf('Signal fit (quant) view: OFF\n');
                end
                obj = obj.updateRender();
            elseif btn == 'M' %Mask select XY
                obj = obj.whileMouseListening_selectMask();
                obj = obj.updateRender();
            elseif btn == 'Z' %toggle selected z boundary top/bottom
                obj.zSetToggle = ~obj.zSetToggle;
                if obj.zSetToggle
                    fprintf('Z set mode: MAX\n');
                else
                    fprintf('Z set mode: MIN\n');
                end
            elseif btn == '^' %Z boundary up
                if obj.zSetToggle
                    Z = size(obj.imageCh,3);
                    zz = obj.scVis.workRegion.z_max;
                    zz = zz + 1;
                    if zz > Z
                        fprintf('Z MAX already at stack top!\n');
                    else
                        obj.scVis.workRegion.z_max = zz;
                        fprintf('Z MAX set to: %d\n', zz);
                        obj = obj.updateRender();
                    end
                else
                    zz = obj.scVis.workRegion.z_min;
                    zz = zz + 1;
                    if zz >  obj.scVis.workRegion.z_max
                        fprintf('Z MIN already at Z MAX!\n');
                    else
                        obj.scVis.workRegion.z_min = zz;
                        fprintf('Z MIN set to: %d\n', zz);
                        obj = obj.updateRender();
                    end
                end
            elseif btn == 'v' %Z boundary down
                if obj.zSetToggle
                    zz = obj.scVis.workRegion.z_max;
                    zz = zz - 1;
                    if zz < obj.scVis.workRegion.z_min
                        fprintf('Z MAX already at Z MIN!\n');
                    else
                        obj.scVis.workRegion.z_max = zz;
                        fprintf('Z MAX set to: %d\n', zz);
                        obj = obj.updateRender();
                    end
                else
                    zz = obj.scVis.workRegion.z_min;
                    zz = zz - 1;
                    if zz < 1
                        fprintf('Z MIN already at bottom of stack!\n');
                    else
                        obj.scVis.workRegion.z_min = zz;
                        fprintf('Z MIN set to: %d\n', zz);
                        obj = obj.updateRender();
                    end
                end
            elseif btn == '+' %Up 1 z slice
                Z = size(obj.imageCh,3);
                zz = obj.currentSlice + 1;
                if zz > Z
                    fprintf('Already at stack top!\n');
                else
                    obj.currentSlice = zz;
                    fprintf('Current slice set to: %d\n', zz);
                    obj = obj.updateRender();
                end
            elseif btn == '-' %Down 1 z slice
                zz = obj.currentSlice - 1;
                if zz < 1
                    fprintf('Already at stack bottom!\n');
                else
                    obj.currentSlice = zz;
                    fprintf('Current slice set to: %d\n', zz);
                    obj = obj.updateRender();
                end
            elseif btn == '=' %Up 10 z slices
                Z = size(obj.imageCh,3);
                if obj.currentSlice == Z
                    fprintf('Already at stack top!\n');
                else
                    zz = obj.currentSlice + 10;
                    if zz > Z; zz = Z; end
                    obj.currentSlice = zz;
                    fprintf('Current slice set to: %d\n', zz);
                    obj = obj.updateRender();
                end
            elseif btn == '_' %Down 10 z slices
                if obj.currentSlice == 1
                    fprintf('Already at stack bottom!\n');
                else
                    zz = obj.currentSlice - 10;
                    if zz < 1; zz = 1; end
                    obj.currentSlice = zz;
                    fprintf('Current slice set to: %d\n', zz);
                    obj = obj.updateRender();
                end
            elseif btn == '<' %Decrease threshold by 1
                %TODO
            elseif btn == '>' %Increase threshold by 1
                %TODO
            elseif btn == '[' %Decrease threshold by 10
                %TODO
            elseif btn == ']' %Increase threshold by 10
                %TODO
            elseif btn == 'R' %Increase xy snap radius by 1
                obj.scVis.matchRad_xy = obj.scVis.matchRad_xy + 1;
                fprintf('XY snap radius set to: %d\n', obj.scVis.matchRad_xy);
            elseif btn == 'r' %Decrease xy snap radius by 1
                obj.scVis.matchRad_xy = obj.scVis.matchRad_xy - 1;
                if obj.scVis.matchRad_xy > 0
                    fprintf('XY snap radius set to: %d\n', obj.scVis.matchRad_xy);
                else
                    obj.scVis.matchRad_xy = 1;
                    fprintf('XY snap radius cannot go below 1.\n');
                end
            end

        end

        %
        function obj = whileMouseListening(obj, x1, y1, b1)

            if ~obj.spotCircleLayerOn
                fprintf('Spot call layer is not visible!\n');
                return;
            end

            if (obj.callViewMode ~= 2) & (obj.callViewMode ~= 4)
                fprintf('Spot selection not available outside reference mode!\n');
                return;
            end

            if b1 == 1
                obj = obj.onLeftClick(x1, y1);
            elseif b1 == 3
                obj = obj.onRightClick(x1, y1);
            end

            loopy = 1;
            while loopy == 1
                [x,y,btn] = ginput_color(1, obj.crosshairColor);
                %fprintf("Click detected at %f,%f\n", x, y);
                if btn == 1
                    obj = obj.onLeftClick(x, y);
                elseif btn == 3
                    obj = obj.onRightClick(x, y);
                else
                    loopy = 0;
                end
            end

            %Save
            obj.scVis = obj.scVis.onRefModeAccept();
            obj = obj.updateRender();
        end

        %
        function obj = whileMouseListening_selectMask(obj)
            %Left clicks to set mask, right click to clear

            %Click 1
            [x,y,btn] = ginput_color(1, obj.crosshairColor);
            if btn == 1
                xPix0 = x;
                yPix0 = y;
            elseif btn == 3
                obj.scVis.workRegion.x_min = 1;
                obj.scVis.workRegion.x_max = size(obj.imageCh, 2);
                obj.scVis.workRegion.y_min = 1;
                obj.scVis.workRegion.y_max = size(obj.imageCh, 1);
                fprintf('XY mask cleared!\n');
                return;
            else
                return;
            end

            %Click 2
            [x,y,btn] = ginput_color(1, obj.crosshairColor);
            if btn == 1
                xPix1 = x;
                yPix1 = y;
            else
                return;
            end

            x0 = round(min(xPix0, xPix1));
            x1 = round(max(xPix0, xPix1));
            y0 = round(min(yPix0, yPix1));
            y1 = round(max(yPix0, yPix1));

            x0 = uint16(max(x0, 1));
            x1 = uint16(min(x1, size(obj.imageCh, 2)));
            y0 = uint16(max(y0, 1));
            y1 = uint16(min(y1, size(obj.imageCh, 1)));

            obj.scVis.workRegion.x_min = x0;
            obj.scVis.workRegion.x_max = x1;
            obj.scVis.workRegion.y_min = y0;
            obj.scVis.workRegion.y_max = y1;

            fprintf('XY Mask set: (%d,%d) to (%d,%d)\n',x0, y0, x1, y1);
        end

        %
        function obj = onLeftClick(obj, x, y)
            fh = obj.figHandle;
            [obj.scVis, fh] = obj.scVis.onRefModeClick(fh, x, y, obj.currentSlice, true, false);
            obj.figHandle = fh;
        end

        %
        function obj = onRightClick(obj, x, y)
            fh = obj.figHandle;
            [obj.scVis, fh] = obj.scVis.onRefModeClick(fh, x, y, obj.currentSlice, false, true);
            obj.figHandle = fh;
        end

    end

    %%
    methods(Static)
    end

end