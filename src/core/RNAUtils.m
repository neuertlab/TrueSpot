%
%%

classdef RNAUtils
    
    methods (Static)
        
        %%
        function border_mask = genBorderMask(dims, rads)
            dimcount = size(dims,2);
            
            if dimcount == 3
                [MY, MX, MZ] = meshgrid(1:dims(1),1:dims(2),1:dims(3));
                border_mask = (MX <= rads(2)) | (MX > (dims(2) - rads(2)));
                border_mask = border_mask | (MY <= rads(1)) | (MY > (dims(1) - rads(1)));
                border_mask = border_mask | (MZ <= rads(3)) | (MZ > (dims(3) - rads(3)));
            else
                [MY, MX] = meshgrid(1:dims(1),1:dims(2));
                border_mask = (MX <= rads(2)) | (MX > (dims(2) - rads(2)));
                border_mask = border_mask | (MY <= rads(1)) | (MY > (dims(1) - rads(1)));
            end
        end
        
        %%
        function img_filtered = medianifyBorder(img_filtered, rads)
            %TODO
        end
        
        %%
        function [d_min, d_max, dtrim_lo, dtrim_hi, needs_trim] = getDimSpotIsolationParams(d_coord, max_d, rad)
            needs_trim = false;
            d_min = d_coord - rad; d_max = d_coord + rad;
            dtrim_lo = 0; dtrim_hi = 0;
            if d_min < 1; dtrim_lo = 1 - d_min; d_min = 1; needs_trim = true; end
            if d_max > max_d; dtrim_hi = d_max - max_d; d_max = max_d; needs_trim = true; end
        end
        
        function spot_data = isolateSpotData2D(src_img, x, y, xrad, yrad)
            Y = size(src_img,1);
            X = size(src_img,2);
            
            [x_min, x_max, xtrim_lo, xtrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(x, X, xrad);
            [y_min, y_max, ytrim_lo, ytrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(y, Y, yrad);
            
            xdim = (xrad * 2) + 1;
            ydim = (yrad * 2) + 1;
            
            spot_data = zeros(ydim, xdim);
            spot_data(ytrim_lo+1:ydim-ytrim_hi, xtrim_lo+1:xdim-xtrim_hi) = ...
                src_img(y_min:y_max, x_min:x_max);
        end
        
        %%
        function spot_data = isolateSpotData(src_img, x, y, z, xyrad, zrad)
            Y = size(src_img,1);
            X = size(src_img,2);
            Z = size(src_img,3);
            
            [x_min, x_max, xtrim_lo, xtrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(x, X, xyrad);
            [y_min, y_max, ytrim_lo, ytrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(y, Y, xyrad);
            [z_min, z_max, ztrim_lo, ztrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(z, Z, zrad);
            
            xydim = (xyrad * 2) + 1;
            zdim = (zrad * 2) + 1;
            
            spot_data = zeros(xydim, xydim, zdim);
            spot_data(ytrim_lo+1:xydim-ytrim_hi, xtrim_lo+1:xydim-xtrim_hi, ztrim_lo+1:zdim-ztrim_hi) = ...
                src_img(y_min:y_max, x_min:x_max, z_min:z_max);
        end
        
    end
    
end