%
%%
classdef RNASpot
    
    properties
        %Position is relative to cell box!
        x;
        y;
        z;
        dropout_thresh = 0;
        
        gauss_fit;
        gfit_slices; %Gauss fit parameters for all z tested. Primary fit is for middle of these.
        adj_gauss_fit;
        
        in_cloud;
        fit_volume;
        nascent_flag = false;
    end
    
    methods
        
        %%
        function obj = initAdjGaussFit(obj)
            obj.adj_gauss_fit = obj.gauss_fit;
        end
        
        %%
        function obj = sumFromGaussFit(obj)
            if ~isempty(obj.gfit_slices)
                obj.fit_volume = 0.0;
                fZ = size(obj.gfit_slices,2);
                for k = 1:fZ
                    k_area = pi * obj.gfit_slices(k).xFWHM * obj.gfit_slices(k).yFWHM;
                    obj.fit_volume = obj.fit_volume + k_area;
                end
            else
                %Just use what we have of middle slice.
                if isempty(obj.gauss_fit)
                    obj.fit_volume = 0.0;
                    return;
                end
                obj.fit_volume = obj.gauss_fit.xFWHM * pi * obj.gauss_fit.yFWHM;
            end
        end
        
        %%
        function obj = saveFittedParamsFromVector(obj, fitted_params)
            fZ = size(fitted_params,2);
            obj.gfit_slices(fZ) = RNAQuant.genEmptyGaussFitParamStruct(true);
            for k = 1:fZ
                mu1 = fitted_params(1,k);
                mu2 = fitted_params(2,k);
                s1 = fitted_params(3,k);
                s2 = fitted_params(4,k);
                A = fitted_params(5,k);
                obj.gfit_slices(k) = RNAQuant.genGaussFitParamStruct(mu1, mu2, s1, s2, A);
            end
        end
        
        %%
        function sim_spot = generateSimSpotFromFit(obj, xy_rad)
            xydim = (xy_rad * 2) + 1;
            
            if ~isempty(obj.gfit_slices)
                zdim = size(obj.gfit_slices,2);
                sim_spot = NaN(xydim,xydim,zdim);
                fslices = obj.gfit_slices;

                mu1 = [fslices.mu1] - obj.x + xy_rad;
                mu2 = [fslices.mu2] - obj.y + xy_rad;
                s1 = [fslices.s1];
                s2 = [fslices.s2];
                A = [fslices.A];
                
                for k = 1:zdim
                    sim_spot(:,:,k) = RNAUtils.generateGaussian2D(...
                        xydim, xydim, mu1(k), mu2(k), ...
                        s1(k), s2(k), A(k));
                end
            else
                sim_spot = NaN(xydim,xydim,1);
                mu1 = obj.gauss_fit.xfit - obj.x + xy_rad;
                mu2 = obj.gauss_fit.yfit - obj.y + xy_rad;
                sim_spot(:,:,1) = RNAUtils.generateGaussian2D(...
                        xydim, xydim, mu1, mu2, ...
                        obj.gauss_fit.xgw, obj.gauss_fit.ygw, obj.gauss_fit.expMInt);
            end
        end
        
    end
    
    methods (Static)
        
        %%
        function [spot_mask, gaussian] = genSpotMask(dims, xw, yw, zw, gthresh)
            if nargin < 5; gthresh = 0.25; end

            spot_mask = false(dims);
            X = size(spot_mask, 2);
            Y = size(spot_mask, 1);
            Z = size(spot_mask, 3);

            mu_x = floor(X./2);
            mu_y = floor(Y./2);
            mu_z = floor(Z./2);

            gaussian = RNAUtils.generateGaussian3D(X, Y, Z, mu_x, mu_y, mu_z, xw, yw, zw, 1.0);
            spot_mask = (gaussian >= 0.45);
        end

        %%
        function sim_spot = generateSimSpotFromFit_Table(spotTable, row, gfit_slices, xy_rad)
            xy_rad_d = double(xy_rad);
            xydim = (xy_rad * 2) + 1;
            
            if ~isempty(gfit_slices)
                zdim = size(gfit_slices,2);
                sim_spot = NaN(xydim,xydim,zdim);

                mu1 = [gfit_slices.mu1] - double(spotTable{row, 'xinit'}) + xy_rad_d;
                mu2 = [gfit_slices.mu2] - double(spotTable{row, 'yinit'}) + xy_rad_d;
                s1 = [gfit_slices.s1];
                s2 = [gfit_slices.s2];
                A = [gfit_slices.A];
                
                for k = 1:zdim
                    sim_spot(:,:,k) = RNAUtils.generateGaussian2D(...
                        xydim, xydim, mu1(k), mu2(k), ...
                        s1(k), s2(k), A(k));
                end
            else
                sim_spot = NaN(xydim,xydim,1);
                mu1 = spotTable{row, 'xfit'} - spotTable{row, 'xinit'} + xy_rad;
                mu2 = spotTable{row, 'yfit'} - spotTable{row, 'yinit'} + xy_rad;
                sim_spot(:,:,1) = RNAUtils.generateGaussian2D(...
                        xydim, xydim, mu1, mu2, ...
                        spotTable{row, 'xgw'}, spotTable{row, 'ygw'}, spotTable{row, 'expMInt'});
            end
        end

        %%
        function rna_spot = newRNASpot()
            rna_spot = RNASpot;
            rna_spot.gfit_slices = RNAQuant.genEmptyGaussFitParamStruct(false);
            rna_spot.gauss_fit = RNAQuant.genGaussFitResultStruct();
            rna_spot.adj_gauss_fit = [];
            rna_spot.x = 0;
            rna_spot.y = 0;
            rna_spot.z = 0;
            rna_spot.in_cloud = false;
            rna_spot.fit_volume = 0.0;
        end
              
    end
    
end