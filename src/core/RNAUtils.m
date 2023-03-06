%
%%

classdef RNAUtils
    
    methods (Static)
        
        %%
        function thresh_idx = findThresholdIndex(thresh_value, thresh_x_tbl)
            %Give thresh_x_tbl as vector (single row)
            
            isge = thresh_x_tbl >= thresh_value;
            if nnz(isge) < 1
                %Nothing found
                thresh_idx = size(thresh_x_tbl,2);
                return;
            end
            
            thresh_idx = find(isge,1);
        end
        
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
        
        %%
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
        function auc_value = calculateAUC(x, y)
            pcount = size(x,2);
            data = NaN(pcount,2);
            data(:,1) = x;
            data(:,2) = y;
            data = unique(data, 'rows');

            %Remove duplicate x values (take higher y value)
            [~,uxi,uxic] = unique(data(:,1));
            uxi_count = size(uxi,2);
            if(uxi_count < pcount)
                keep_bool = true(1,pcount);
                uidiff = diff(uxic);
                zmode = false;
                for p = 1:pcount-1
                    %p-1 is index in original
                    if uidiff(p) == 0
                        keep_bool(p) = false;
                        if ~zmode
                            keep_bool(p-1) = false;
                            zmode = true;
                        end
                    else
                        if zmode
                            keep_bool(p-1) = true;
                            zmode = false;
                        end
                    end
                end

                data = data(find(keep_bool), :);
            end

            %Add point on left to bring to y axis
            data_adj = NaN(pcount+1,2);
                data_adj(2:pcount,:) = data(:,:);
            data_adj(1,1) = 0.0;
            data_adj(1,2) = data_adj(1,2);

            ddiff = NaN(pcount,2);
            ddiff(:,1) = diff(data_adj(:,1));
            ddiff(:,2) = diff(data_adj(:,2));
            
            area_a = 0.5 .* abs(ddiff(:,1)) .* abs(ddiff(:,1));
            area_b = abs(ddiff(:,1)) .* min(data_adj(1:pcount,2), data_adj(2:pcount+1,2));
            auc_value = sum(area_a) + sum(area_b);
        end
        
    end
    
end