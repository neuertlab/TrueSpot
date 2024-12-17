%
%%

classdef RNAUtils
    
    methods (Static)
        
        %%
        function bool = isTableVariable(myTable, varName)
            varNames = myTable.Properties.VariableNames;
            bool = ismember(varName, varNames);
        end

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
            if (thresh_x_tbl(thresh_idx) > thresh_value)
                if thresh_idx > 1; thresh_idx = thresh_idx - 1; end
            end
        end
        
        %%
        function [xx, yy] = spotCountFromCallTable(call_table, include_trimmed, min_th, max_th)
            if nargin < 2; include_trimmed = false; end
            if nargin < 3; min_th = 0; end
            if nargin < 4; max_th = 0; end

            xx = [];
            yy = [];
            if isempty(call_table); return; end

            %Get th range.
            allth = call_table{:, 'dropout_thresh'};
            allth = allth(allth ~= 0);
            allth = unique(allth);
            allth = sort(allth);
            
            if min_th < 1
                tmin = double(allth(1));
            else
                tmin = double(min_th);
            end

            if max_th < 1
                tmax = double(allth(size(allth, 1)));
            else
                tmax = double(max_th);
            end

            allth_d = diff(allth);
            allth_d = allth_d(allth_d ~= 0);
            tintr = double(min(allth_d, [], 'all', 'omitnan'));

            xx = [tmin:tintr:tmax];

            if include_trimmed | ~RNAUtils.isTableVariable(call_table, 'is_trimmed_out')
                yy = sum(call_table{:, 'dropout_thresh'} >= xx, 1);
            else
                yy = sum(call_table{~call_table{:,'is_trimmed_out'}, 'dropout_thresh'} >= xx, 1);
            end

        end

        %%
        function spot_table = spotTableFromCallTable(call_table, include_trimmed, min_th, max_th)
            if nargin < 2; include_trimmed = false; end
            if nargin < 3; min_th = 0; end
            if nargin < 4; max_th = 0; end

            spot_table = [];
            [xx, yy] = RNAUtils.spotCountFromCallTable(call_table, include_trimmed, min_th, max_th);
            
            if isempty(xx); return; end
            T = size(xx, 2);
            spot_table = NaN(T, 2);
            spot_table(:,1) = double(xx);
            spot_table(:,2) = double(yy);
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
            d_coord = int32(d_coord);
            max_d = int32(max_d);
            rad = int32(rad);
            
            d_min = d_coord - rad; d_max = d_coord + rad;
            dtrim_lo = max((1 - d_min), 0);
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
        function auc_value = calculateAUC(x, y)
            auc_value = NaN;
            if isempty(x); return; end
            if isempty(y); return; end
            
            dim1 = size(x,1);
            dim2 = size(x,2);
            
            if(dim1 > dim2)
                x = transpose(x);
                y = transpose(y);
            end
            
            %Remove any records where EITHER x or y is nan
            badrec_bool = (isnan(x) | isnan(y));
            badcount = nnz(badrec_bool);
            if badcount > 0
                if badcount >= size(x,2); return; end
                goodrecs = find(~badrec_bool);
                x = x(goodrecs);
                y = y(goodrecs);
            end
            
            %Sort by y, then by x
            [~, ysort_idx] = sort(y);
            x_sorted = x(ysort_idx);
            y_sorted = y(ysort_idx);

            [~, xsort_idx] = sort(x_sorted);
            x_sorted = x_sorted(xsort_idx);
            y_sorted = y_sorted(xsort_idx);

            %Remove duplicate x values
            ptcount = size(x_sorted,2);
            keep_bool = false(1,ptcount);
            keep_bool(1:(ptcount-1)) = (x_sorted(1:ptcount-1) ~= x_sorted(2:ptcount));
            keep_bool(ptcount) = true;
            
            keep_idx = find(keep_bool);
            x_sorted = x_sorted(keep_idx);
            y_sorted = y_sorted(keep_idx);

            %Add end points
            if(x_sorted(1) > 0.0)
                x_sorted = [0 0 x_sorted];
                y_sorted = [0 y_sorted(1) y_sorted];
            else
                x_sorted = [0 x_sorted];
                y_sorted = [0 y_sorted];
            end
            ptcount = size(x_sorted,2);

            if(y_sorted(ptcount) > 0.0)
                x_sorted = [x_sorted x_sorted(ptcount)];
                y_sorted = [y_sorted 0];
            end

            ply = polyshape(x_sorted, y_sorted);
            auc_value = area(ply);
            
            %DEBUG
%             figure(1);
%             plot(ply);
%             
        end
        
        %%
        function printVectorToTextFile(fhandle, vec, fmtstr, newline)
            vsize = size(vec, 2);
            fprintf(fhandle, '{');
            for i = 1:vsize
                if i ~= 1; fprintf(fhandle, ','); end
                fprintf(fhandle, fmtstr, vec(i));
            end
            
            if newline
                fprintf(fhandle, '}\n');
            else 
                fprintf(fhandle, '}');
            end
        end

        %%
        function boolRes = isInMask3(mask, x, y, z)
            if isvector(x)
                s1 = size(x,1);
                s2 = size(x,2);
                if s1 > s2
                    eCount = s1;
                    x = x';
                    y = y';
                    z = z';
                else
                    eCount = s2;
                end

                boolRes = false(eCount, 1);
                for i = 1:eCount; boolRes(i) = mask(y(i),x(i),z(i)); end

            else
                boolRes = mask(y,x,z);
            end
        end

        %%
        function boolRes = isInMask2(mask, x, y)
            if isvector(x)
                s1 = size(x,1);
                s2 = size(x,2);
                if s1 > s2
                    eCount = s1;
                    x = x';
                    y = y';
                else
                    eCount = s2;
                end

                boolRes = false(eCount, 1);
                for i = 1:eCount; boolRes(i) = mask(y(i),x(i)); end

            else
                boolRes = mask(y,x);
            end
        end

    end
    
end