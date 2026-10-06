%
%%

classdef RNAUtils
    
    methods (Static)
        
        %%
        function rangeStruct = genCoordRangeStruct(includeZ)
            rangeStruct = struct();
            rangeStruct.x_min = NaN;
            rangeStruct.x_max = NaN;
            rangeStruct.y_min = NaN;
            rangeStruct.y_max = NaN;
            if includeZ
                rangeStruct.z_min = NaN;
                rangeStruct.z_max = NaN;
            end
        end

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

        %%
        function imgname = imageNameFromFile(filePath)
            MAX_NAME_LEN = 48;
            [~, fname, ~] = fileparts(filePath);

            %Clean up so less likely to have future file naming issues
            imgname = replace(fname, '.', '');
            imgname = replace(imgname, ' ', '_');
            if length(imgname) > MAX_NAME_LEN
                imgname = imgname(1:MAX_NAME_LEN);
            end
        end

        %%
        function dead_pix_info = detectDeadPixels(in_img, verbose)
            %Adapted from RNA_Threshold_Common.saveDeadPixels
            %But without mandatory save...
            if nargin < 3; verbose = false; end

            Y = size(in_img,1);
            X = size(in_img,2);
            Z = size(in_img,3);
            A = Y * X;
            in_img = double(in_img);

            dead_pix_info = struct('recurring_all_count', 0, 'random_count', 0, 'recurring_count', 0);
            dead_pix_info.idims = struct('x', X, 'y', Y, 'z', Z);
            dead_pix_info.recurring_pixels = [];
            dead_pix_info.recurring_pixels_all = [];
            dead_pix_info.rand_recur_pix = [];

            slice_cut = 4/8;
            slice_check = round(min(Z,max(Z./4, 6)));
            hi_pixels = cell(slice_check, 1);
            hi_pixels_all = [];
            if verbose; fprintf("Determining recurring pixels...\n"); end
    
            w2 = [-1 -1 -1;...
                  -1 +8 -1;...
                  -1 -1 -1;];
      
            %For some subgroup of z slices, find pixels with unusually
            %   high values after edge filter is applied.
            if size(in_img,1) > 1
                for i = 1:slice_check
                    slice = imfilter(in_img(:,:,i),w2);
                    cutoff = mean(slice(:)) + (3 * std(slice(:)));
                    temp_pix = find(slice > cutoff);
                    hi_pixels{i,1} = temp_pix;
                    hi_pixels_all = cat(1, hi_pixels_all, temp_pix);
                end
            end
    
            %For each processed slice, pick the same number
            %   of pixels randomly
            rand_hi = [];
            for j = 1:slice_check
                temp_hi = hi_pixels{j,1};
                rand_hi_temp = randsample(A, size(temp_hi,2));
                rand_hi_temp = rand_hi_temp';
                rand_hi = cat(1, rand_hi, rand_hi_temp);
            end
            
            %See how often each pixel appears in random selection
            [hist_pix_rand, ~] = histcounts(double(rand_hi(:)), A);
            hist_pix_rand = transpose(hist_pix_rand);
            
            %Note which occur more than slice_cut proportion of the time.
            dead_pix_info.rand_recur_pix = find(hist_pix_rand >= slice_check*slice_cut);
            dead_pix_info.random_count = size(dead_pix_info.rand_recur_pix,1);
            if verbose; fprintf("%d randomly recurring pixels\n", dead_pix_info.random_count); end
            
            %Repeat with high value pixels
            [hist_pix, ~] = histcounts(double(hi_pixels_all(:)), A);
            hist_pix = transpose(hist_pix);
            dead_pix_info.recurring_pixels_all = find(hist_pix >= slice_check*slice_cut);
            dead_pix_info.recurring_all_count = size(dead_pix_info.recurring_pixels_all,1);
            if verbose; fprintf("%d total recurring pixels\n", dead_pix_info.recurring_all_count); end
            
            %Remove border pixels
            recmtx = NaN(dead_pix_info.recurring_all_count,3);
            recmtx(:,1) = dead_pix_info.recurring_pixels_all(:,1); %1D
            recmtx(:,2) = floor((recmtx(:,1) - 1)./Y) + 1; %x
            recmtx(:,3) = mod((recmtx(:,1) - 1), Y) + 1; %y
            testmtx = recmtx(:,2) > 1;
            testmtx = testmtx & (recmtx(:,2) < X);
            testmtx = testmtx & (recmtx(:,3) < Y);
            testmtx = testmtx & (recmtx(:,3) > 1);
            [keeprows, ~] = find(testmtx);
            dead_pix_info.recurring_pixels = recmtx(keeprows,1);
            dead_pix_info.recurring_count = size(recurring_pixels,1);
            if verbose; fprintf("%d non-border recurring pixels\n", dead_pix_info.recurring_count); end
        end

        %%
        function clean_img = cleanDeadPixels(in_img, dead_pix_info, verbose)
            %Adapted from RNA_Threshold_Common.cleanDeadPixels
            %But without mandatory save...
            if nargin < 3; verbose = false; end
            
            rp_count = size(dead_pix_info.recurring_pixels,1);
            cleaned_count = 0;
    
            %Get size(s) and save. Cleans up code.
            dim1 = size(in_img,1); %Height
            dim2 = size(in_img,2); %Width
            dim3 = size(in_img,3); %Depth
    
            clean_img = in_img;
                
            %The notice messages are kept from previous code
            %'Averaging recurring pixels'
            if verbose
                fprintf("Averaging recurring pixels...\t")
                tic
            end
            for k = 1:dim3
                slice = in_img(:,:,k); %[uint16[D][D]] 2D slice of image
                for c = 1:rp_count
                    p = dead_pix_info.recurring_pixels(c); %[int] 1D coordinate of bad pixel
                    
                    y = mod((p-1), dim1) + 1;
                    if y <= 1; continue; end
                    if y >= dim1; continue; end
                    
                    x = floor((p-1)/dim1) + 1;
                    if x <= 1; continue; end
                    if x >= dim2; continue; end
                    
                    west = p - dim1;
                    east = p + dim1;
                    surr_pix = [p-1, p+1, west, east, west-1, east-1, east+1, west+1];
                    %N S W E NW NE SE SW
                    slice(p) = mean(slice(surr_pix));
                    cleaned_count = cleaned_count+1;
                end
                clean_img(:,:,k) = slice;
            end
            if verbose
                toc; 
                fprintf("Recurring voxels cleaned: %d\n", cleaned_count);
            end
        end

        %%
        function [img_filtered, deadpix_info] = applyLoGFilter(in_img, gaussian_rad, dead_pix_detect_bool)
            if nargin < 2; gaussian_rad = 7; end
            if nargin < 3; dead_pix_detect_bool = true; end

            %Detect dead pixels...
            deadpix_info = RNAUtils.detectDeadPixels(in_img, true);

            %Pre-filtering...
            IMG3D = uint16(in_img);
            if (dead_pix_detect_bool)
                IMG3D = RNAUtils.cleanDeadPixels(IMG3D, deadpix_info, true);
            end
            
            %Actually, let's rescale the *raw* image.
            irmin = min(IMG3D, [], 'all');
            irmax = max(IMG3D, [], 'all');
            irrng = irmax - irmin;
            if irrng < 256
                fprintf("WARNING: Raw image has low dynamic range. Rescale triggered!\n");
                %Linear rescale to 0-511
                IMG3D = ((IMG3D - irmin) .* 512) ./ irrng;
            end

            img_filtered = RNA_Threshold_Common.applyGaussianFilter(IMG3D, gaussian_rad, 2);
            img_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(img_filtered);
            img_filtered = RNA_Threshold_Common.blackoutBorders(img_filtered, gaussian_rad+1, 0);
            
            %Rescale the filtered image if not enough range
            ifmin = min(img_filtered, [], 'all');
            ifmax = max(img_filtered, [], 'all');
            ifrng = ifmax - ifmin;
            if ifrng < 25
                fprintf("WARNING: Filtered image has low dynamic range. Rescale triggered!\n");
                %Linear rescale to 0-255
                img_filtered = ((img_filtered - ifmin) .* 255) ./ ifrng;
            end
            
            img_filtered = uint16(img_filtered); %To reduce memory usage. Note that on dim images this can have a dramatic effect.
        end

    end
    
end