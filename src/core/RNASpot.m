%
%%
classdef RNASpot
    
    properties
        %Position is relative to cell box!
        x;
        y;
        z;
        
        gauss_fit;
        gfit_slices; %Gauss fit parameters for all z tested. Primary fit is for middle of these.
        adj_gauss_fit;
        
        in_cloud;
        fit_volume;
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
                
                for k = 1:zdim
                    sim_spot(:,:,k) = RNAUtils.generateGaussian2D(...
                        xydim, xydim, obj.gfit_slices(k).mu1, obj.gfit_slices(k).mu2, ...
                        obj.gfit_slices(k).s1, obj.gfit_slices(k).s2, obj.gfit_slices(k).A);
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