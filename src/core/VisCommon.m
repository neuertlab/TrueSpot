%%
classdef VisCommon
    
    %%
    properties
        colortbl_red;
        colortbl_magenta;
        colortbl_yellow;
        colortbl_cyan;
        colortbl_white;

        colorBelowFar;
        colorBelowNear;
        colorAboveFar;
        colorAboveNear;

        idims;

        isInitialized = false;
    end

    %%
    methods
        
        %%
        function obj = initializeMe(obj, idims)
            if nargin < 2; idims = []; end

            if isempty(idims)
                idims = struct('x', 0, 'y', 0, 'z', 0);
            end

            obj = obj.generateColorTables(idims.z + 1);
            obj.idims = idims;

            obj.colorBelowFar = [0.60, 0.50, 0.0];
            obj.colorBelowNear = [1.0, 1.0, 0.0];
            obj.colorAboveFar = [0.0, 0.50, 0.20];
            obj.colorAboveNear = [0.0, 1.0, 0.0];

            obj.isInitialized = true;
        end

        %%
        function obj = generateColorTables(obj, levels)
            obj.colortbl_red = VisCommon.generateDarkColors(levels, [1.0, 0.0, 0.0]);
            obj.colortbl_magenta = VisCommon.generateDarkColors(levels, [1.0, 0.0, 1.0]);
            obj.colortbl_yellow = VisCommon.generateDarkColors(levels, [1.0, 1.0, 0.0]);
            obj.colortbl_cyan = VisCommon.generateDarkColors(levels, [0.0, 1.0, 1.0]);
            obj.colortbl_white = VisCommon.generateDarkColors(levels, [1.0, 1.0, 1.0]);
        end


    end

    %%
    methods(Static)

         %%
        function clrmtx = gen_z_color_matrix(color_tbl, coord_tbl, rel_z, range)
            maxidx = size(color_tbl,1);
            incr = 1;
            
            if range > 0
                incr = maxidx./range;
            end
            
            color_idxs = min((abs(rel_z - coord_tbl(:,3)) + 1) * incr, maxidx);
            clrmtx = color_tbl(color_idxs,:);
        end

        %%
        function color_tbl = generateDarkColors(levels, base_color)

            levels = double(levels+1);
            color_tbl = NaN(levels,3);
            konst = 1/log10(levels);

            r = base_color(1,1);
            g = base_color(1,2);
            b = base_color(1,3);

            color_tbl(1,1) = r;
            color_tbl(1,2) = g;
            color_tbl(1,3) = b;

            %r_div = r/levels;
            %g_div = g/levels;
            %b_div = b/levels;

            for l = 2:levels
                %r = r-r_div;
                r = (1 - (log10(l) * konst)) * base_color(1,1);
                color_tbl(l,1) = r;

                %g = g-g_div;
                g = (1 - (log10(l) * konst)) * base_color(1,2);
                color_tbl(l,2) = g;

                %b = b-b_div;
                b = (1 - (log10(l) * konst)) * base_color(1,3);
                color_tbl(l,3) = b;
            end

        end

        %%
        function color_tbl = generateLightColors(levels, base_color)

            levels = double(levels+1);
            color_tbl = NaN(levels,3);
            konst = 1/log10(levels);

            r = base_color(1,1);
            g = base_color(1,2);
            b = base_color(1,3);

            color_tbl(1,1) = r;
            color_tbl(1,2) = g;
            color_tbl(1,3) = b;

            for l = 2:levels
                %r = r-r_div;
                r = (log10(l) * konst) + base_color(1,1);
                color_tbl(l,1) = min(r, 1.0);

                %g = g-g_div;
                g = (log10(l) * konst) + base_color(1,2);
                color_tbl(l,2) = min(g, 1.0);

                %b = b-b_div;
                b = (log10(l) * konst) + base_color(1,3);
                color_tbl(l,3) = min(b, 1.0);
            end

        end

        %%
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

        %%
        function compImg = compositeNewChannel(baseImageRGB, overlayImageBW, overlayLUT, rescaleOverlay)
            if nargin < 4; rescaleOverlay = true; end
            rgbOverlay = VisCommon.bw2rgb(overlayImageBW, overlayLUT, rescaleOverlay);
            rgbOverlay = double(rgbOverlay) ./ 255.0;
            baseDbl = double(baseImageRGB) ./ 255.0;
            baseDbl = baseDbl .* (1.0 - rgbOverlay);
            baseDbl = baseDbl + rgbOverlay;
            compImg = uint8(round(baseDbl .* 255.0));
        end

        %%
        function compImg = compositeMaskOverlay(baseImageRGB, mask, color, alpha, doOutline, outlineDiskSize)
            if nargin < 6; outlineDiskSize = 3; end
            
            compImg = baseImageRGB;
            if isempty(mask); return; end
            if nnz(mask) < 1; return; end
            baseDbl = double(baseImageRGB) ./ 255.0;
            if doOutline
                mask = bwperim(mask, 8);
                se = strel('disk',outlineDiskSize);
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

        %%
        function lut = genGreyscaleLUT()
            [~, lut] = meshgrid(1:3, 1:256);
            lut = lut - 1;
            lut = lut ./ 255.0;
        end
    end
end