%Methods for RNA cloud detection (mask) given image and set of points
%Blythe Hospelhorn
%Modified from code written by Ben Kesler & Gregor Neuert
%Version 1.1.0
%Updated Dec 6, 2022

%Update Log:
%   1.0.0 | 21.03.12
%       Init Doc
%   1.1.0 | 22.12.06
%       Outputs coordinates of individual cloud boxes when detected.
%       Also fixed a couple issues Ben brought up.

%TODO
%   1. Sort input coords by brightness, and don't bother looking for clouds
%   for spots already part of previous cloud.
%   2. Dead pix
%   3. ratio > obj.d_thresh ?

%%
classdef RNA_Clouds
    
    properties
        max_dist = 50; %Max distance (radius in pixels) to look for clouds around a point
        d_thresh = 0.97; %Deriv. threshold
        thresh_devs = 3; % # stdevs pixel intensity must be from median to be called in cloud
    end
    
    methods
        
        %%
        % Set the max pixel distance this RNA_Clouds instance should use as
        % a search radius relative to a start point when scanning for
        % clouds.
        %
        % ARGS
        %   value (int) - Value to set search radius to. Should be greater
        %       than or equal to 1.
        %
        function obj = setMaxDist(obj, value)
            obj.max_dist = value;
        end
        
        %%
        % Set the derivative threshold this RNA_Clouds instance should use
        % when looking for cloud edges using change in pixel intensity.
        %
        % ARGS
        %   value (double) - Value to set deriv threshold to.
        %
        function obj = setDerivThreshold(obj, value)
            obj.d_thresh = value;
        end
        
        %%
        % Set the threshold for number of standard deviations from the
        % median above which a pixel should be called as part of the cloud/
        %
        % ARGS
        %   value (int) - Value to set stdev threshold to. Should be greater
        %       than or equal to 0.
        %
        function obj = setThresholdStdevs(obj, value)
            obj.thresh_devs = value;
        end
        
        %%
        % Generate a boolean mask of the provided image where pixels with
        % "true" values represent pixels that are part of a detected signal
        % cloud, and pixels with "false" values are not.
        % This method variant automatically determines whether input is 2D
        % or 3D and calls appropriate function.
        %
        % ARGS
        %   img (num[Y][X]([Z])) - Image to detect cloud signals from
        %   nuc_coord_table (int[n][2+]) - Table of image coodinates to
        %       start cloud searches at.
        %   z_adj (num) - Z adjustment factor, or ratio of distance between
        %       adjacent pixels in the z direction relative to xy. This
        %       parameter is ignored if input is 2D.
        %
        % RETURN
        %   cloud_mask (bool[Y][X]([Z])) - Boolean matrix of the same
        %       dimensions of input image that can be used as a mask
        %       to determine which pixels are part of signal clouds.
        %
        function [cloud_mask, cloud_boxes] = detectClouds(obj, img, nuc_coord_table, z_adj)
            cloud_mask = [];
            cloud_boxes = [];

            if isempty(nuc_coord_table); return; end

            %Count columns in coord table to see if it's 2D or 3D
            coord_dims = size(nuc_coord_table, 2);
            if(coord_dims == 2)
                [cloud_mask, cloud_boxes] = obj.detectClouds2D(img, nuc_coord_table);
            elseif (coord_dims == 3)
                [cloud_mask, cloud_boxes] = obj.detectClouds3D(img, nuc_coord_table, z_adj);
            end
            
        end
        
        %%
        % Generate a boolean mask of the provided image where pixels with
        % "true" values represent pixels that are part of a detected signal
        % cloud, and pixels with "false" values are not.
        %
        % ARGS
        %   img (num[Y][X]) - 2D Image to detect cloud signals from
        %   nuc_coord_table (int[n][2]) - Table of 2D image coodinates to
        %       start cloud searches at.
        %
        % RETURN
        %   cloud_mask (bool[Y][X]) - 2D boolean matrix of the same
        %       dimensions of input image that can be used as a mask
        %       to determine which pixels are part of signal clouds.
        %
        function [cloud_mask, cloud_boxes] = detectClouds2D(obj, img, nuc_coord_table)
            cloud_boxes = [];
            X = size(img, 2);
            Y = size(img, 1);
            cloud_mask = logical(Y,X);
            
            spotcount = size(nuc_coord_table, 1);
            temp_mem = NaN(spotcount, 4);
            temp_next = 1;
            for s = 1:spotcount
                
                sx = nuc_coord_table(s, 1);
                sy = nuc_coord_table(s, 2);
                
                if cloud_mask(sy, sx); continue; end
                
                %Determine window radius
                
                %Make a version of the image where each pix value is dist
                %from spot
                %TODO - Maybe only perform on a subchunk of the image w/
                %size ~100x100x100? Would that be faster?
                
                [xx, yy] = meshgrid((1:X) - sx, (1:Y) - sy);
                s_dist = sqrt(xx.^2 + yy.^2);
                winrad = 0;

                last_avg = -1.0;
                for r = 1:obj.max_dist
                    %Mask
                    rmask = (s_dist <= r);
                    masked_img = immultiply(img, rmask);
                    masked_img(masked_img == 0) = NaN;
                    this_avg = nanmean(masked_img, 'all');
                    
                    if last_avg >= 0.0
                        ratio = this_avg/last_avg;
                        if ratio > obj.d_thresh
                            winrad = r;
                        end
                    end
                    
                    last_avg = this_avg;
                    if winrad > 0
                        break;
                    end
                end
                
                %The window radius is set
                %Get the median and stdev of the pix values in the window
                %rmask = (s_dist <= winrad);
                %windowed = immultiply(img_filter, rmask);
                %masked_img(masked_img == 0) = NaN;
                
                %Actually, let's see if we can just use the ones already
                %generated by the loop...
                w_med = nanmedian(masked_img, 'all');
                w_std = nanstd(masked_img, 0, 'all');
                c_thresh = w_med + (obj.thresh_devs * w_std);
                
                local_mask = masked_img > c_thresh;
                rprops = regionprops(local_mask,'BoundingBox','Area');
                cloud_box = rprops.BoundingBox;
                temp_mem(temp_next,:) = cloud_box(1,:);
                temp_next = temp_next + 1;
                cloud_mask = or(cloud_mask, local_mask);
            end
            if temp_next > 1
                cloud_boxes = uint16(temp_mem(1:temp_next-1,:));
            else
                cloud_boxes = [];
            end
            clear temp_mem;
        end
        
        %%
        % Generate a boolean mask of the provided image where pixels with
        % "true" values represent pixels that are part of a detected signal
        % cloud, and pixels with "false" values are not.
        %
        % ARGS
        %   img (num[Y][X][Z]) - Image to detect cloud signals from
        %   nuc_coord_table (int[n][3]) - Table of 3D image coodinates to
        %       start cloud searches at.
        %   z_adj (num) - Z adjustment factor, or ratio of distance between
        %       adjacent pixels in the z direction relative to xy.
        %
        % RETURN
        %   cloud_mask (bool[Y][X][Z]) - 3D boolean matrix of the same
        %       dimensions of input image that can be used as a mask
        %       to determine which pixels are part of signal clouds.
        %
        function [cloud_mask, cloud_boxes] = detectClouds3D(obj, img, nuc_coord_table, z_adj)
            %max_dist = 50;
            %d_thresh = 0.97;
            %thresh_devs = 3;
            
            X = size(img, 2);
            Y = size(img, 1);
            Z = size(img, 3);
            cloud_mask = false(Y,X,Z);
            
            spotcount = size(nuc_coord_table, 1);
            temp_mem = NaN(spotcount, 6);
            temp_next = 1;
            for s = 1:spotcount
                
                sx = double(nuc_coord_table(s, 1));
                sy = double(nuc_coord_table(s, 2));
                sz = double(nuc_coord_table(s, 3));
                
                %Skip if already in a cloud.
                if cloud_mask(sy, sx, sz); continue; end
                
                %Determine window radius
                
                %Make a version of the image where each pix value is dist
                %from spot
                
                [xx, yy, zz] = meshgrid((1:X) - sx, (1:Y) - sy, ((1:Z) - sz) * z_adj);
                s_dist = sqrt(xx.^2 + yy.^2 + zz.^2);
                winrad = 0;

                last_avg = -1.0;
                for r = 1:obj.max_dist
                    %Mask
                    rmask = (s_dist <= r);
                    masked_img = immultiply(img, rmask);
                    masked_img(masked_img == 0) = NaN;
                    this_avg = mean(masked_img, 'all', 'omitnan');
                    
                    if last_avg >= 0.0
                        ratio = this_avg/last_avg;
                        if ratio > obj.d_thresh
                            winrad = r;
                        end
                    end
                    
                    last_avg = this_avg;
                    if winrad > 0
                        break;
                    end
                end
                
                %The window radius is set
                %Get the median and stdev of the pix values in the window
                %rmask = (s_dist <= winrad);
                %windowed = immultiply(img_filter, rmask);
                %masked_img(masked_img == 0) = NaN;
                
                %Actually, let's see if we can just use the ones already
                %generated by the loop...
                w_med = median(masked_img, 'all', 'omitnan');
                w_std = std(masked_img, 0, 'all', 'omitnan');
                c_thresh = w_med + (obj.thresh_devs * w_std);
                
                local_mask = masked_img > c_thresh;
                if nnz(local_mask) <= 0
                    continue;
                end
                
                rprops = regionprops3(local_mask,'BoundingBox');
                cloud_box_all = rprops.BoundingBox;
                cloud_box_all(:,4) = cloud_box_all(:,4) + cloud_box_all(:,1);
                cloud_box_all(:,5) = cloud_box_all(:,5) + cloud_box_all(:,2);
                cloud_box_all(:,6) = cloud_box_all(:,6) + cloud_box_all(:,3);
                
                cloud_box = NaN(1,6);
                cloud_box(1,1) = max(floor(min(cloud_box_all(:,1))), 1);
                cloud_box(1,2) = max(floor(min(cloud_box_all(:,2))), 1);
                cloud_box(1,3) = max(floor(min(cloud_box_all(:,3))), 1);
                cloud_box(1,4) = min(ceil(max(cloud_box_all(:,4))),X) - cloud_box(1,1);
                cloud_box(1,5) = min(ceil(max(cloud_box_all(:,5))),Y) - cloud_box(1,2);
                cloud_box(1,6) = min(ceil(max(cloud_box_all(:,6))),Z) - cloud_box(1,3);
                
                temp_mem(temp_next,:) = cloud_box(1,:);
                temp_next = temp_next + 1;
                cloud_mask = or(cloud_mask, local_mask);
            end
            
            if temp_next > 1
                cloud_boxes = uint16(temp_mem(1:temp_next-1,:));
            else
                cloud_boxes = [];
            end
            
            clear temp_mem;
        end
           
    end
end