%
%%
classdef ImageUtils

    methods(Static)

        %-------------------------------------- Image --------------------------------------

        %%
        function [cmin, cmax] = autoPickRawContrastRange(channelStackData)
            maxproj = double(max(channelStackData, [], 3, 'omitnan'));
            cmin = min(maxproj(:));
            cmax = median(maxproj(:)) + round(10 * std(maxproj(:)));
        end

        %%
        function [cmin, cmax] = autoPickFiltContrastRange(channelStackData)
            maxproj = double(max(channelStackData, [], 3, 'omitnan'));
            cmin = median(maxproj(:)) - round(0 * std(maxproj(:)));
            cmax = median(maxproj(:)) + round(10 * std(maxproj(:)));
        end

        %%
        function image8 = linearRescaleToU8FullRange(imageIn)
            %Basically, minimum (whatever it is) is set to 0, and maximum
            %to 255
            valMin = min(imageIn, [], 'all', 'omitnan');
            imageIn = imageIn - valMin; %Normalize min to 0
            valMax = max(imageIn, [], 'all', 'omitnan');
            imageIn = imageIn ./ valMax; %Normalize max to 1

            image8 = uint8(round(imageIn .* 255.0));
        end
        

        %%
        function image8 = linearRescaleU16ToU8(image16)
            image = double(image16) ./ double(0xffff);
            image8 = uint8(round(image .* 255.0));
        end

        %%
        function imageDbl = linearRescaleTo01FullRange(imageIn)
            imageIn = double(imageIn);
            if nnz(imageIn) > 0
                valMin = min(imageIn, [], 'all', 'omitnan');
                imageIn = imageIn - valMin; %Normalize min to 0
                valMax = max(imageIn, [], 'all', 'omitnan');
                imageDbl = imageIn ./ valMax; %Normalize max to 1
            else
                imageDbl = imageIn;
            end
        end

        %%
        function imageDbl = linearRescaleTo01ContrastRange(imageIn, cMin, cMax)
            imageIn = double(imageIn);
            if cMin == cMax
                X = size(imageIn, 2);
                Y = size(imageIn, 1);
                %Z = size(imageIn, 3);
                if cMax == 0
                    imageDbl = zeros(Y,X);
                else
                    imageDbl = ones(Y,X);
                end
            else
                imageIn = imageIn - cMin; %Normalize min to 0
                imageDbl = imageIn ./ cMax; %Normalize max to 1
            end
        end

        %%
        function clutSet = loadUtilCLUTs()
            datPath = './utils/clut.mat';
            if isfile(datPath)
                load(datPath, 'clutSet');
            else
                clutSet = [];
            end
        end

        %-------------------------------------- Render --------------------------------------

        %%
        function outputRGBDbl = alphaLayerOn01Image(baseRGBDbl, layerDbl, layerColor, layerMasterAlpha)
            %Base and output are all 3 channel (RGB) 2D images with
            %intensity range of 0.0 - 1.0
            %Layer is single channel
            outputRGBDbl = baseRGBDbl;

            layerR = (layerDbl .* layerMasterAlpha) .* (layerColor(1));
            layerG = (layerDbl .* layerMasterAlpha) .* (layerColor(2));
            layerB = (layerDbl .* layerMasterAlpha) .* (layerColor(3));

            baseWeightR = 1.0 - layerR;
            baseWeightG = 1.0 - layerG;
            baseWeightB = 1.0 - layerB;

            r = (baseRGBDbl(:,:,1) .* baseWeightR) + layerR;
            g = (baseRGBDbl(:,:,2) .* baseWeightG) + layerG;
            b = (baseRGBDbl(:,:,3) .* baseWeightB) + layerB;

            outputRGBDbl(:,:,1) = r;
            outputRGBDbl(:,:,2) = g;
            outputRGBDbl(:,:,3) = b;
        end

        %%
        function figHandle = drawSpotCircles(figHandle, callTable, threshold, z, color)
            %If z is 0, use MIP (and all z)
            callTable = callTable(callTable{:, 'dropout_thresh'} >= threshold, :);
            if z > 0
                callTable = callTable(callTable{:, 'isnap_z'} == z, :);
            end

            xx = callTable{:, 'isnap_x'};
            yy = callTable{:, 'isnap_y'};
           
            figure(figHandle);
            plot(xx, yy, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none');
        end

    end
end