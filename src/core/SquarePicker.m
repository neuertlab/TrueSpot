%
%%
classdef SquarePicker

    %%
    properties
        imageMaxproj = [];
        cellMask = [];
        nucMask = [];

        sampleLUT = SquarePicker.genGreyscaleLUT();
        cellMaskColor = [1 0 1];
        nucMaskColor = [1 1 0];

        figHandle = [];
        currentSelection = RNAUtils.genCoordRangeStruct(false);

        crosshairColor = [1 1 0];
        loopBreaker = false;
    end

     %%
    methods
        %% ================== GUI ===============================

        function obj = updateRender(obj)
            [bwImageScaled, ~, ~] = SquarePicker.rescaleSampleChannel(obj.imageMaxproj);
            iRender = SquarePicker.bw2rgb(bwImageScaled, obj.sampleLUT, false);
            if ~isempty(obj.cellMask)
                iRender = SquarePicker.compositeMaskOverlay(iRender, ...
                    obj.cellMask, obj.cellMaskColor, 0.25, true);
            end

            if ~isempty(obj.nucMask)
                iRender = SquarePicker.compositeMaskOverlay(iRender, ...
                    obj.nucMask, obj.nucMaskColor, 0.25, true);
            end

            figure(obj.figHandle);
            imshow(iRender); hold on;

            squareValid = SquarePicker.xyRegionValid(obj.currentSelection);
            if squareValid
                x0 = obj.currentSelection.x_min;
                x1 = obj.currentSelection.x_max;
                y0 = obj.currentSelection.y_min;
                y1 = obj.currentSelection.y_max;

                rectangle('Position',[x0 y0 (x1-x0) (y1-y0)],'EdgeColor','r','LineWidth',4);
            end
        end

        function obj = launchFigureGUI(obj)
            obj = obj.updateRender();
            obj.loopBreaker = false;

            while ~obj.loopBreaker
                w = waitforbuttonpress;
                if w
                    obj = obj.onReadyKey();
                end
            end
        end

        %
        function obj = onReadyKey(obj)
            [x,y,btn] = ginput_color(1, obj.crosshairColor);

            if btn == 1 %Mouse click
                %Loop until any key pressed that isn't mouse
                obj = obj.whileMouseListening_selectMask(x,y);
                obj = obj.updateRender();
            elseif btn == 3 %Right click
                obj.currentSelection.x_min = NaN;
                obj.currentSelection.x_max = NaN;
                obj.currentSelection.y_min = NaN;
                obj.currentSelection.y_max = NaN;
                obj = obj.updateRender();
            elseif btn == 'x'
                obj.loopBreaker = true;
            end

        end

        %
        function obj = whileMouseListening_selectMask(obj, x, y)
            %Left clicks to set mask, right click to clear

            %Click 1
            xPix0 = x;
            yPix0 = y;

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

            %Make it square.
            w = x1 - x0 + 1;
            cx = x0 + (w./2);
            h = y1 - y0 + 1;
            cy = y0 + (h./2);

            if w > h
                hw = w./2;
                y0 = round(cy - hw);
                y1 = round(cy + hw);
                clear hw;
            elseif h > w
                hh = h./2;
                x0 = round(cx - hh);
                x1 = round(cx + hh);
                clear hh;
            end

            x0 = uint16(max(x0, 1));
            x1 = uint16(min(x1, size(obj.imageMaxproj, 2)));
            y0 = uint16(max(y0, 1));
            y1 = uint16(min(y1, size(obj.imageMaxproj, 1)));

            obj.currentSelection.x_min = x0;
            obj.currentSelection.x_max = x1;
            obj.currentSelection.y_min = y0;
            obj.currentSelection.y_max = y1;

            fprintf('New Square picked: (%d,%d) to (%d,%d)\n',x0, y0, x1, y1);
        end
    end

    %%
    methods(Static)

        %
        function rgbImage = bw2rgb(bwImage, lut, doRescale)
            if nargin < 3; doRescale = true; end
            X = size(bwImage, 2);
            Y = size(bwImage, 1);
            rgbImage = zeros(Y,X,3);

            if doRescale
                bwImage = double(bwImage);
                bwMin = min(bwImage, [], 'all', 'omitnan');
                bwImage = bwImage - bwMin;
                bwMed = median(bwImage, 'all', 'omitnan');
                bwStd = std(bwImage, 0, 'all', 'omitnan');
                bwMax = bwMed + round(10 * bwStd);
                %bwMax = max(bwImage, [], 'all', 'omitnan');
                bwImage = bwImage ./ bwMax;
            end

            bwImage = round(bwImage .* 255.0);
            bwImage = bwImage + 1;
            bwImage = min(bwImage, 256);
            bwImage = max(bwImage, 1);

            for c = 1:3
                cmult = lut(bwImage, c);
                rgbImage(:,:,c) = reshape(cmult, Y, X);
            end

            rgbImage = rgbImage .* 255.0;
            rgbImage = uint8(rgbImage);
        end

        %
        function compImg = compositeNewChannel(baseImageRGB, overlayImage, overlayLUT, rescaleOverlay)
            if nargin < 4; rescaleOverlay = true; end
            rgbOverlay = bw2rgb(overlayImage, overlayLUT, rescaleOverlay);
            rgbOverlay = double(rgbOverlay) ./ 255.0;
            baseDbl = double(baseImageRGB) ./ 255.0;
            baseDbl = baseDbl .* (1.0 - rgbOverlay);
            baseDbl = baseDbl + rgbOverlay;
            compImg = uint8(round(baseDbl .* 255.0));
        end

        %
        function compImg = compositeMaskOverlay(baseImageRGB, mask, color, alpha, doOutline)
            baseDbl = double(baseImageRGB) ./ 255.0;
            if doOutline
                mask = bwperim(mask, 8);
                se = strel('disk',5);
                mask = imdilate(mask,se);
                clear se
            end
            Y = size(mask, 1);
            X = size(mask, 2);
            maskrgb = zeros(Y,X,3);
            maskrgb(:,:,1) = double(mask) .* color(1);
            maskrgb(:,:,2) = double(mask) .* color(2);
            maskrgb(:,:,3) = double(mask) .* color(3);
            maskrgb = maskrgb .* alpha;
            baseDbl = baseDbl .* (1.0 - maskrgb);
            baseDbl = baseDbl + maskrgb;
            compImg = uint8(round(baseDbl .* 255.0));
        end
   
        %
        function [bwImageScaled, bwShift, bwRange] = rescaleSampleChannel(bwImage)
            bwImage = double(bwImage);
            %bwMin = prctile(bwImage, 2, 'all');
            bwMin = min(bwImage, [], 'all', 'omitnan');
            bwImage = bwImage - bwMin;
            bwMed = median(bwImage, 'all', 'omitnan');
            bwStd = std(bwImage, 0, 'all', 'omitnan');
            bwMax = bwMed + round(10 * bwStd);
            %bwMax = max(bwImage, [], 'all', 'omitnan');
            bwImageScaled = bwImage ./ bwMax;
            bwShift = bwMin;
            bwRange = bwMax;
        end

        %
        function lut = genGreyscaleLUT()
            [~, lut] = meshgrid(1:3, 1:256);
            lut = lut - 1;
            lut = lut ./ 255.0;
        end

        %
        function bool = xyRegionValid(regionInfo)
            bool = false;
            if isnan(regionInfo.x_min); return; end
            if isnan(regionInfo.x_max); return; end
            if isnan(regionInfo.y_min); return; end
            if isnan(regionInfo.y_max); return; end
            if (regionInfo.x_min < 1); return; end
            if (regionInfo.x_max < 1); return; end
            if (regionInfo.y_min < 1); return; end
            if (regionInfo.y_max < 1); return; end

            bool = true;
        end
    end


end