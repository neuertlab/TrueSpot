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
            rad_y = rads(1); Y = size(img_filtered,1);
            rad_x = rads(2); X = size(img_filtered,2);
            rad_z = rads(3); Z = size(img_filtered,3);
            mvalue = median(img_filtered(rad_y+1:Y-rad_y,rad_x+1:X-rad_x,rad_z+1:Z-rad_z), 'all');

            img_filtered(1:rad_y+1,:,:) = mvalue;
            img_filtered(Y-rad_y:Y,:,:) = mvalue;
            img_filtered(:,1:rad_x+1,:) = mvalue;
            img_filtered(:,X-rad_x:X,:) = mvalue;
            img_filtered(:,:,1:rad_z+1) = mvalue;
            img_filtered(:,:,Z-rad_z:Z) = mvalue;
        end
        
        %%
        function [d_min, d_max, dtrim_lo, dtrim_hi, needs_trim] = getDimSpotIsolationParams(d_coord, max_d, rad)
            d_min = d_coord - rad; d_max = d_coord + rad;
            dtrim_lo = max(1 - d_min,0);
            dtrim_hi = max(d_max - max_d,0);
            d_min = max(d_min, 1);
            d_max = min(d_max, max_d);
            needs_trim = isscalar(d_coord) & (dtrim_lo > 0 | dtrim_hi > 0);
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