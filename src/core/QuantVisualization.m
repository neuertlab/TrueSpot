%
%%
classdef QuantVisualization

    %%
    properties
        maxProj = true;
        cells = SingleCell.empty();

        spotrad_xy = 4;
        spotrad_z = 2;

        spotColor = [1.000 1.000 0.000];
        cloudColor = [0.000 1.000 1.000];
        spotNascentColor = [0.000 1.000 0.000];
        cloudNascentColor = [0.300 1.000 1.000];

        countLabelColor = [1.000 0.000 0.000];
        showCountLabels = false;

        rgbCellRenders = [];
    end

    %%
    methods

        %%
        function [obj, okay] = prerenderCells(obj)
            cellCount = size(obj.cells, 2);
            
            %Determine scale intensity
            %TODO
            peakIntensity = 1;
            for c = 1:cellCount
                myCell = obj.cells(c);
                if isempty(myCell.spotTable) & ~isempty(myCell.spots)
                    myCell = myCell.convertSpotStorage();
                end

                maxi = 0;
                if ~isempty(myCell.spotZFits)
                    fits = myCell.spotZFits;
                    sCount = size(fits, 2);
                    for j = 1:sCount
                        myFits = fits{j};
                        mm = max([myFits.A], [], 'all', 'omitnan');
                        if mm > maxi; maxi = mm; end
                    end
                    clear mm fits
                else
                    if ~isempty(myCell.spotTable)
                        maxi = max(myCell.spotTable{:, 'expMInt'});
                    end
                end

                if maxi > peakIntensity; peakIntensity = maxi; end
                obj.cells(c) = myCell;
            end
            clear maxi;


            obj.rgbCellRenders = cell(1,cellCount);
            for c = 1:cellCount
                [obj, okay] = obj.renderCellRGB(obj.cells(c), peakIntensity);
                if ~okay; return; end
            end
        end
 
        %%
        function [obj, okay] = renderCellRGB(obj, myCell, peakIntensity)
            okay = false;
            if isempty(obj); return; end
            if isempty(myCell); return; end

            %cellRender = struct();
            %cellRender.cellBounds = myCell.cell_loc;

            r = zeros(myCell.cell_loc.height, myCell.cell_loc.width, myCell.cell_loc.depth);
            g = zeros(myCell.cell_loc.height, myCell.cell_loc.width, myCell.cell_loc.depth);
            b = zeros(myCell.cell_loc.height, myCell.cell_loc.width, myCell.cell_loc.depth);
            a = zeros(myCell.cell_loc.height, myCell.cell_loc.width, myCell.cell_loc.depth);

            [alphaMask, okay] = obj.renderCellCloudMask(myCell, true, peakIntensity);
            if ~okay; return; end
            if ~isempty(alphaMask)
                [r, g, b, a] = QuantVisualization.applyNewMask(alphaMask, obj.cloudNascentColor, r, g, b, a);
            end

            [alphaMask, okay] = obj.renderCellCloudMask(myCell, false, peakIntensity);
            if ~okay; return; end
            if ~isempty(alphaMask)
                [r, g, b, a] = QuantVisualization.applyNewMask(alphaMask, obj.cloudColor, r, g, b, a);
            end

            [alphaMask, okay] = obj.renderCellSpotsMask(myCell, true, peakIntensity);
            if ~okay; return; end
            if ~isempty(alphaMask)
                [r, g, b, a] = QuantVisualization.applyNewMask(alphaMask, obj.spotNascentColor, r, g, b, a);
            end

            [alphaMask, okay] = obj.renderCellSpotsMask(myCell, false, peakIntensity);
            if ~okay; return; end
            if ~isempty(alphaMask)
                [r, g, b, a] = QuantVisualization.applyNewMask(alphaMask, obj.spotColor, r, g, b, a);
            end

            cellRender = struct('cellNumber', myCell.cell_number);
            cellRender.cellBounds = myCell.cell_loc;
            cellRender.red = r;
            cellRender.green = g;
            cellRender.blue = b;
            cellRender.alpha = a;

            obj.rgbCellRenders{myCell.cell_number} = cellRender;
            okay = true;
        end

        %%
        function [alphaMask, okay] = renderCellSpotsMask(obj, myCell, nascentFlag, peakIntensity)
            %Generates a matrix the size of the cell region image with
            %values from 0-1
            okay = false;
            alphaMask = [];
            if isempty(obj); return; end
            if isempty(myCell); return; end

            %Get spots
            if isempty(myCell.spotTable) & ~isempty(myCell.spots)
                myCell = myCell.convertSpotStorage();
            end

            if isempty(myCell.spotTable)
                okay = true;
                return; %Return empty mask so don't waste time applying it
            end

            %Filter spots
            keepRows = (myCell.spotTable{:, 'nascent_flag'} == nascentFlag);
            if nnz(keepRows) < 1
                okay = true;
                return; %Return empty mask
            end
            spotList = myCell.spotTable(keepRows, :);
            zFits = myCell.spotZFits(keepRows);

            %Cell dims
            cX = myCell.cell_loc.width;
            cY = myCell.cell_loc.height;
            cZ = myCell.cell_loc.depth;
            alphaMask = zeros(cY, cX, cZ);

            %Spot locations
            spotCount = size(spotList, 1);
            sx = int32(spotList{:, 'xinit'}');
            sy = int32(spotList{:, 'yinit'}');
            sz = int32(spotList{:, 'zinit'}');

            %Determine render boxes
            [tx0, tx1, rx0, rx1] = QuantVisualization.getBoxBounds(sx, obj.spotrad_xy, cX);
            [ty0, ty1, ry0, ry1] = QuantVisualization.getBoxBounds(sy, obj.spotrad_xy, cY);
            [tz0, tz1, rz0, rz1] = QuantVisualization.getBoxBounds(sz, obj.spotrad_z, cZ);

            for i = 1:spotCount
                %subimg = spotList(i).generateSimSpotFromFit(obj.spotrad_xy);

                subimg = RNASpot.generateSimSpotFromFit_Table(spotList, i, zFits{i}, obj.spotrad_xy);
                alphaMask(ty0(i):ty1(i), tx0(i):tx1(i), tz0(i):tz1(i))...
                    = alphaMask(ty0(i):ty1(i), tx0(i):tx1(i), tz0(i):tz1(i))...
                    + subimg(ry0(i):ry1(i), rx0(i):rx1(i), rz0(i):rz1(i));
            end

            alphaMask = alphaMask ./ peakIntensity;
            okay = true;
        end

        %%
        function [alphaMask, okay] = renderCellCloudMask(obj, myCell, nascentFlag, peakIntensity)
            okay = false;
            alphaMask = [];
            if isempty(obj); return; end
            if isempty(myCell); return; end

            if nascentFlag
                nval = 1;
            else
                nval = 0;
            end

            cloudList = myCell.getAllClouds(-1, nval);
            if isempty(cloudList)
                okay = true;
                return; %Return empty mask so don't waste time applying it
            end

            %Cell dims
            cX = myCell.cell_loc.width;
            cY = myCell.cell_loc.height;
            cZ = myCell.cell_loc.depth;
            alphaMask = zeros(cY, cX, cZ);

            cloudCount = size(cloudList, 2);
            cloudLocs = [cloudList.cloud_box];
            clx0 = [cloudLocs.left]; clx1 = [cloudLocs.right];
            cly0 = [cloudLocs.top]; cly1 = [cloudLocs.bottom];
            clz0 = [cloudLocs.z_bottom]; clz1 = [cloudLocs.z_top];

            for i = 1:cloudCount
                subimg = double(cloudList(i).cloud_data);
                subimg(~cloudList(i).cloud_mask) = 0;
                alphaMask(cly0(i):cly1(i), clx0(i):clx1(i), clz0(i):clz1(i))...
                    = alphaMask(cly0(i):cly1(i), clx0(i):clx1(i), clz0(i):clz1(i))...
                    + subimg;
            end

            alphaMask = alphaMask ./ peakIntensity;
            okay = true;
        end

        %%
        function imgOut = applyToImage(obj, imgIn, z)
            %Don't forget to check if input is already RGB first
            if ndims(imgIn) > 2
                %See if it's already RGB
                if isa(imgIn, 'uint8') & (size(imgIn, 3) == 3)
                    red = double(imgIn(:,:,1)) ./ 255.0;
                    green = double(imgIn(:,:,2)) ./ 255.0;
                    blue = double(imgIn(:,:,3)) ./ 255.0;
                else
                    %Assume 3D BW
                    if obj.maxProj
                        bw = double(max(imgIn, [], 3));
                    else
                        bw = double(imgIn(:,:,z));
                    end
                    mm = max(bw, [], 'all', 'omitnan');
                    bw = bw ./ mm;
                    red = bw;
                    green = bw;
                    blue = bw;
                end
            else
                %Assume it's 2D BW
                mm = max(imgIn, [], 'all', 'omitnan');
                bw = imgIn ./ mm;
                red = bw;
                green = bw;
                blue = bw;
            end

            X = size(red, 2);
            Y = size(red, 1);
            imgOut = NaN(Y,X,3);
            imgOut(:,:,1) = red;
            imgOut(:,:,2) = green;
            imgOut(:,:,3) = blue;

            cellCount = size(obj.rgbCellRenders, 2);
            for c = 1:cellCount
                cellRender = obj.rgbCellRenders{c};
                x0 = cellRender.cellBounds.left;
                x1 = cellRender.cellBounds.right;
                y0 = cellRender.cellBounds.top;
                y1 = cellRender.cellBounds.bottom;
                z0 = cellRender.cellBounds.z_bottom;
                %z1 = cellRender.cellBounds.z_top;

                rr = red(y0:y1,x0:x1);
                gg = green(y0:y1,x0:x1);
                bb = blue(y0:y1,x0:x1);

                rc = cellRender.red;
                gc = cellRender.green;
                bc = cellRender.blue;
                aa = cellRender.alpha;

                if obj.maxProj
                    rc = max(rc, [], 3);
                    gc = max(gc, [], 3);
                    bc = max(bc, [], 3);
                    aa = max(aa, [], 3);
                else
                    zst = z - z0 + 1;
                    rc = rc(:,:,zst);
                    gc = gc(:,:,zst);
                    bc = bc(:,:,zst);
                    aa = aa(:,:,zst);
                end

                %DEBUG----
%                 figure(500);
%                 clf;
%                 imshow(rc, []);
%                 figure(501);
%                 clf;
%                 imshow(gc, []);
%                 figure(502);
%                 clf;
%                 imshow(bc, []);
%                 figure(503);
%                 clf;
%                 imshow(aa, []);

                aInv = 1.0 - aa;
                rr = (rr .* aInv) + (rc .* aa);
                gg = (gg .* aInv) + (gc .* aa);
                bb = (bb .* aInv) + (bc .* aa);

                imgOut(y0:y1,x0:x1,1) = rr;
                imgOut(y0:y1,x0:x1,2) = gg;
                imgOut(y0:y1,x0:x1,3) = bb;
            end

        end
    end

    %%
    methods(Static)

        %%
        function [main0, main1, sub0, sub1] = getBoxBounds(center, rad, maxDim)
            count = size(center, 2);
            subdim = (rad * 2) + 1;

            main0 = double(center - rad);
            main1 = double(center + rad);

            sub0 = ones(1, count);
            sub1 = repmat(subdim, 1, count);

            hangamt = max(1 - main0, 0);
            main0 = main0 + hangamt;
            sub0 = sub0 + hangamt;
            hangamt = max(main1 - maxDim, 0);
            main1 = main1 - hangamt;
            sub1 = sub1 - hangamt;

            main0 = uint16(round(main0));
            main1 = uint16(round(main1));
            sub0 = uint16(round(sub0));
            sub1 = uint16(round(sub1));
        end

        %%
        function [r, g, b, a] = applyNewMask(alphaMask, color, r, g, b, a)
            rr = alphaMask .* color(1);
            gg = alphaMask .* color(2);
            bb = alphaMask .* color(3);

            alphaInv = 1.0 - alphaMask;
            r = (r .* alphaInv) + rr;
            g = (g .* alphaInv) + gg;
            b = (b .* alphaInv) + bb;
            a = alphaMask + a; %I guess??
        end

    end

end