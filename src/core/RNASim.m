%
%%

classdef RNASim
    
    methods (Static)

        %----------------------------- Structs -----------------------------

        %%
        function paramStruct = genSpotGroupGenParamStruct()
            paramStruct = struct();
            paramStruct.spotCount = 0;

            %Intensity of a single spot
            paramStruct.meanPeak = NaN;
            paramStruct.stdPeak = NaN;

            %Gaussian rads (All in pixels)
            %(FWHM)
            paramStruct.gwMeanX = NaN;
            paramStruct.gwMeanY = NaN;
            paramStruct.gwMeanZ = NaN;
            paramStruct.gwStdX = NaN;
            paramStruct.gwStdY = NaN;
            paramStruct.gwStdZ = NaN;

            %Multi-spot parameters
            paramStruct.brightnessThreshold = NaN;
            paramStruct.singleIntensity = NaN;

            %From SingleCell
            paramStruct.spotTable = table.empty();
            paramStruct.spotZFits = [];

            %Misc
            paramStruct.randomSeed = NaN;
        end

        %----------------------------- Spot Image Gen -----------------------------

        %%
        function gauss_spot = generateGaussian2D(xdim, ydim, mu_x, mu_y, w_x, w_y, peak)
            [XX,YY] = meshgrid(1:1:xdim,1:1:ydim);
            
            x_factor = (XX - mu_x - 1).^2;
            y_factor = (YY - mu_y - 1).^2;
            xw_factor = 2 * (w_x^2);
            yw_factor = 2 * (w_y^2);
            
            gauss_spot = peak * exp(-((x_factor ./ xw_factor) + (y_factor ./ yw_factor)));
        end

        %%
        function gauss_spot = generateGaussian3D(xdim, ydim, zdim, mu_x, mu_y, mu_z, w_x, w_y, w_z, peak)
            [XX,YY,ZZ] = meshgrid(1:1:xdim,1:1:ydim,1:1:zdim);
            
            x_factor = (XX - mu_x - 1).^2;
            y_factor = (YY - mu_y - 1).^2;
            z_factor = (ZZ - mu_z - 1).^2;
            xw_factor = 2 * (w_x^2);
            yw_factor = 2 * (w_y^2);
            zw_factor = 2 * (w_z^2);
            
            gauss_spot = peak * exp(-((x_factor ./ xw_factor) + (y_factor ./ yw_factor) + (z_factor ./ zw_factor)));
        end

        %%
        function [genSpotTable, simData] = reconCloudSpots(cloudMask, cloudData, spotGenParams)
            genSpotTable = table.empty();
            if isempty(spotGenParams); return; end

            X = size(cloudMask, 2);
            Y = size(cloudMask, 1);
            Z = size(cloudMask, 3);

            remainingSpots = spotGenParams.spotCount;
            %Start with known spots
            %TODO
        end

    end
end