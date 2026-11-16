%
%%
classdef ImageUtils

    methods(Static)

        %-------------------------------------- Image --------------------------------------

        %%
        function [cmin, cmax] = autoPickRawContrastRange(channelStackData8)
            maxproj = double(max(channelStackData8, [], 3, 'omitnan'));
            cmin = min(maxproj(:));
            cmax = median(maxproj(:)) + round(10 * std(maxproj(:)));
        end

        %%
        function [cmin, cmax] = autoPickFiltContrastRange(channelStackData8)
            maxproj = double(max(channelStackData8, [], 3, 'omitnan'));
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

        %-------------------------------------- Render --------------------------------------

        %%
        function outputRGBDbl = alphaLayerOn01Image(baseRGBDbl, layerDbl, layerColor, layerMasterAlpha)
            %TODO
            %Base and output are all 3 channel (RGB) 2D images with
            %intensity range of 0.0 - 1.0
            %Layer is single channel

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