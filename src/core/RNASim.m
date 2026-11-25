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
        function imageOut = addSimSpotToImage(simSpot, imageIn, c_x, c_y, c_z)
            %TODO
        end

        %%
        function [genSpotTable, simData] = reconCloudSpots(cloudMask, cloudData, spotGenParams)
            genSpotTable = table.empty();
            if isempty(spotGenParams); return; end

            X = size(cloudMask, 2);
            Y = size(cloudMask, 1);
            Z = size(cloudMask, 3);

            %Prep output
            simData = zeros(Y,X,Z);
            %simMask = false(Y,X,Z);

            varNames = {'center_x', 'center_y', 'center_z', 'xgw', 'ygw', 'zgw', 'amplitude'};
            varTypes = {'single', 'single', 'single', 'single', 'single', 'single', 'single'};

            table_size = [spotGenParams.spotCount size(varNames,2)];
            genSpotTable = table('Size', table_size, 'VariableTypes', varTypes, 'VariableNames', varNames);

            remainingSpots = spotGenParams.spotCount;
            spotPos = 1;
            %Start with known spots
            if ~isempty(spotGenParams.spotTable)
                %Look for spots inside this cloud
                sx = spotGenParams.spotTable{:, 'xinit'};
                sy = spotGenParams.spotTable{:, 'yinit'};
                sz = spotGenParams.spotTable{:, 'zinit'};

                inCloud = false(1, spotGenParams.spotCount);
                for j = 1:spotGenParams.spotCount
                    inCloud(j) = cloudMask(sy, sx, sz);
                end
                clear j sx sy sz

                if nnz(inCloud) > 0
                    cloudSpots = spotGenParams.spotTable(inCloud, :);
                    cloudFits = spotGenParams.spotZFits(inCloud);
                    csCount = size(cloudSpots, 1);

                    for i = 1:csCount
                        %Determine how many sub-spots
                        sint = spotGenParams.spotTable{:, 'TotFitInt'};
                        if (sint > paramStruct.brightnessThreshold)
                            %Split into multiple
                            %TODO
                        else
                            %Render as single as-is
                            xy_rad = floor(cloudSpots{i, 'xFWHM'} .* 2);
                            sim_spot = RNASpot.generateSimSpotFromFit_Table(cloudSpots, i, cloudFits{i}, xy_rad);

                            c_x = cloudSpots{i, 'xinit'};
                            c_y = cloudSpots{i, 'yinit'};
                            c_z = cloudSpots{i, 'zinit'};
                            simData = RNASim.addSimSpotToImage(sim_spot, simData, c_x, c_y, c_z);
                            clear xy_rad sim_spot

                            %TODO
                            genSpotTable{spotPos, 'center_x'} = c_x;
                            genSpotTable{spotPos, 'center_y'} = c_y;
                            genSpotTable{spotPos, 'center_z'} = c_z;
                            genSpotTable{spotPos, 'xgw'} = cloudSpots{i, 'xgw'};
                            genSpotTable{spotPos, 'ygw'} = cloudSpots{i, 'ygw'};
                            %TODO - calculate zgw

                            clear c_x c_y c_z

                            remainingSpots = remainingSpots - 1;
                            spotPos = spotPos + 1;
                        end
                        clear sint
                    end
                    clear csCount cloudSpots cloudFits i
                end
                clear inCloud
            end

            %TODO
        end

    end
end