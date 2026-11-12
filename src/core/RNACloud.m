%
%%
classdef RNACloud
    
    properties
        cloud_box; %Coordinates of cloud, relative to its cell
        cloud_mask;
        cloud_data; 
        
        is_nuc;
        is_cyto;
        nascent_flag = false;
        
        cloud_volume; %Just the number of pix/vox cloud encompasses.
        max_intensity; %TODO Add
        total_intensity;
    end
    
    methods
        
        %%
        function obj = updateStats(obj)
            if isempty(obj.cloud_mask)
                obj.cloud_volume = 0.0;
                obj.max_intensity = 0;
                obj.total_intensity = 0;
                return; 
            end

            obj.cloud_volume = sum(obj.cloud_mask, 'all');

            if ~isempty(obj.cloud_data)
                maskedDat = immultiply(obj.cloud_data, obj.cloud_mask);
                obj.total_intensity = sum(maskedDat, 'all', 'omitnan');
                obj.max_intensity = max(maskedDat, [], 'all', 'omitnan');
            else
                obj.max_intensity = 0;
                obj.total_intensity = 0;
            end
        end

        %%
        function obj = recalculateVolume(obj)
            if isempty(obj.cloud_mask)
                obj.cloud_volume = 0.0;
                return; 
            end
            obj.cloud_volume = sum(obj.cloud_mask, 'all');
        end
        
        %%
        function pkg = packageForSave(obj)
            pkg = struct();
            pkg.cloud_box = obj.cloud_box;
            pkg.cloud_mask = obj.cloud_mask;
            pkg.cloud_data = obj.cloud_data;
            pkg.is_nuc = obj.is_nuc;
            pkg.is_cyto = obj.is_cyto;
            pkg.nascent_flag = obj.nascent_flag;
            pkg.cloud_volume = obj.cloud_volume;
            pkg.max_intensity = obj.max_intensity;
            pkg.total_intensity = obj.total_intensity;
        end
    
    end
    
    methods (Static)

        %---------------------------- Individual Clouds ---------------------------

        %%
        function rna_cloud = readFromSavePackage(pkg)
            rna_cloud = RNACloud;
            if isempty(pkg); return; end
            rna_cloud.cloud_box = pkg.cloud_box;
            rna_cloud.cloud_mask = pkg.cloud_mask;
            rna_cloud.cloud_data = pkg.cloud_data;
            rna_cloud.is_nuc = pkg.is_nuc;
            rna_cloud.is_cyto = pkg.is_cyto;
            rna_cloud.nascent_flag = pkg.nascent_flag;
            rna_cloud.cloud_volume = pkg.cloud_volume;
            rna_cloud.max_intensity = pkg.max_intensity;
            rna_cloud.total_intensity = pkg.total_intensity;
        end

        %%
        function rna_cloud = newRNACloud()
            rna_cloud = RNACloud;
            rna_cloud.cloud_box = SingleCell.generateRecPrismStruct(1,1,1,1,1,1);
            rna_cloud.cloud_mask = [];
            rna_cloud.is_nuc = false;
            rna_cloud.is_cyto = false;
            rna_cloud.cloud_volume = 0.0;
            rna_cloud.total_intensity = 0.0;
            rna_cloud.cloud_data = [];
        end

        %%
        function rna_cloud = cloudFromBoxStruct(cloudBoxStruct)
            rna_cloud = RNACloud.newRNACloud();
            rna_cloud.cloud_box.left = cloudBoxStruct.x0;
            rna_cloud.cloud_box.right = cloudBoxStruct.x1;
            rna_cloud.cloud_box.top = cloudBoxStruct.y0;
            rna_cloud.cloud_box.bottom = cloudBoxStruct.y1;
            rna_cloud.cloud_box.z_bottom = cloudBoxStruct.z0;
            rna_cloud.cloud_box.z_top = cloudBoxStruct.z1;

            rna_cloud.cloud_box.width = abs(cloudBoxStruct.x1 - cloudBoxStruct.x0) + 1;
            rna_cloud.cloud_box.height = abs(cloudBoxStruct.y1 - cloudBoxStruct.y0) + 1;
            rna_cloud.cloud_box.depth = abs(cloudBoxStruct.z1 - cloudBoxStruct.z0) + 1;
            rna_cloud.cloud_box.width_i = round(rna_cloud.cloud_box.width);
            rna_cloud.cloud_box.height_i = round(rna_cloud.cloud_box.height);
            rna_cloud.cloud_box.depth_i = round(rna_cloud.cloud_box.depth);
        end

        %---------------------------- Cloud Detection ---------------------------

        %%
        function cloudDetectSettings = genCloudDetectSettingsStruct()
            cloudDetectSettings = struct();
            cloudDetectSettings.xyMaxDist = 50; %Max distance (radius in pixels) to look for clouds around a point
            cloudDetectSettings.zMaxDist = 5; %Max distance (radius in pixels) to look for clouds around a point
            cloudDetectSettings.avgRatioThresh = 0.97; %TODO: This may need to be tweaked
            %cloudDetectSettings.zAvgRatioThresh = 0.97;
            cloudDetectSettings.stdevCount = 3; % # stdevs pixel intensity must be within from median to be called in cloud
            cloudDetectSettings.useMaxProj = true; %if false, use plane that spot is in
            cloudDetectSettings.zxRatio = 0.0; %if non-zero use 3D region (if applicable)
        end

        %%
        function cloudBox = genCloudBoxStruct()
            cloudBox = struct();
            cloudBox.x0 = 0;
            cloudBox.x1 = 0;
            cloudBox.y0 = 0;
            cloudBox.y1 = 0;
            cloudBox.z0 = 0;
            cloudBox.z1 = 0;
        end

        %%
        function cloudMask = detectCloudAroundLocalMax(imageRaw, settings, spot_x, spot_y, spot_z)
            if isempty(settings)
                settings = RNACloud.genCloudDetectSettingsStruct();
            end

            X = size(imageRaw, 2);
            Y = size(imageRaw, 1);
            Z = size(imageRaw, 3);

            if Z > 1
                if settings.useMaxProj
                    sliceRaw = max(imageRaw, [], 3, 'omitnan');
                else
                    sliceRaw = imageRaw(:,:,spot_z);
                end
            else
                sliceRaw = imageRaw;
            end

            spot_x = double(spot_x);
            spot_y = double(spot_y);
            spot_z = double(spot_z);
            %[xx, yy] = meshgrid((1:X) - spot_x, (1:Y) - spot_y);
            %s_dist = sqrt(xx.^2 + yy.^2);
            [xx, yy, zz] = meshgrid((1:X) - spot_x, (1:Y) - spot_y, ((1:Z) - spot_z) * settings.zxRatio);
            s_dist3 = sqrt(xx.^2 + yy.^2 + zz.^2);

            xx = xx(:,:,1);
            yy = yy(:,:,1);
            s_dist2 = sqrt(xx.^2 + yy.^2);
            winrad = 0;
            clear xx yy zz

            %!!! DEBUG
            % figure(1); clf;
            % imshow(sliceRaw, []);
            % hold on;
            % plot(spot_x, spot_y, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');

            last_avg = -1.0;
            for r = 1:settings.xyMaxDist
                if (settings.zxRatio > 0.0) & (Z > 1)
                    rmask3 = (s_dist3 <= r);
                    img_rmask3 = imageRaw;
                    img_rmask3(~rmask3) = NaN;

                    this_avg = mean(img_rmask3, 'all', 'omitnan');
                else
                    %Mask
                    rmask2 = (s_dist2 <= r);
                    img_rmask2 = sliceRaw;
                    img_rmask2(~rmask2) = NaN;

                    %Take the average of everything inside this radius
                    this_avg = mean(img_rmask2, 'all', 'omitnan');
                end

                %Take ratio to average of everything inside last tested
                %radius
                if last_avg >= 0.0
                    ratio = this_avg/last_avg;

                    %Look for leveloff
                    if ratio > settings.avgRatioThresh
                        winrad = r;
                    end

                    %Look for dropoff
                    % if ratio < settings.avgRatioThresh
                    %     winrad = r;
                    % end
                end

                last_avg = this_avg;
                if winrad > 0
                    break;
                end
            end
            clear r rmask2 img_rmask2 ratio last_avg this_avg
           
            %Find a cloud intensity threshold by finding the median and
            %standard dev within the xy window on the chosen slice/proj
            if winrad == 0; winrad = settings.xyMaxDist; end
            
            %!!! DEBUG
            % cx0 = max(spot_x - winrad, 1);
            % cx1 = min(spot_x + winrad, X);
            % cy0 = max(spot_y - winrad, 1);
            % cy1 = min(spot_y + winrad, Y);
            % rectangle('Position', [cx0 cy0 (cx1-cx0) (cy1-cy0)], 'Curvature', [1 1], 'EdgeColor', 'r');

            rmask2 = (s_dist2 <= winrad);
            img_rmask2 = sliceRaw;
            img_rmask2(~rmask2) = NaN;

            if (settings.zxRatio > 0.0) & (Z > 1)
                rmask3 = s_dist3 <= winrad;
                img_rmask3 = imageRaw;
                img_rmask3(~rmask3) = NaN;

                w_med = median(img_rmask3, 'all', 'omitnan');
                w_std = std(img_rmask3, 0, 'all', 'omitnan');
                c_thresh = w_med + (settings.stdevCount * w_std);
            else
                rmask3 = repmat(rmask2, [1 1 Z]);
                w_med = median(img_rmask2, 'all', 'omitnan');
                w_std = std(img_rmask2, 0, 'all', 'omitnan');
                c_thresh = w_med + (settings.stdevCount * w_std);
            end
            clear img_rmask2 w_med w_std sliceRaw s_dist2 s_dist3

            %Apply threshold in 3D to get mask
            cloudMask = and(rmask3, imageRaw >= c_thresh);

            if nnz(cloudMask) > 0
                %Keep only the region containing the original maximum
                cmlbl = bwlabeln(cloudMask);
                keepr = cmlbl(spot_y, spot_x, spot_z);
                if(keepr > 0)
                    cloudMask = (cmlbl == keepr);
                else
                    cloudMask(:) = 0;
                end
                clear cmlbl keepr
            end

            %!!! DEBUG
            % cmMax = max(cloudMask, [], 3, 'omitnan');
            % figure(2); clf;
            % imshow(cmMax, []);
        end

        %%
        function cloudBoxes = cloudBoxesFromMask(cloudMask)
            Z = size(cloudMask, 3);

            if Z > 1
                rprops = regionprops3(cloudMask, 'BoundingBox');
                rCount = size(rprops, 1);

                cloudBoxes(rCount) = RNACloud.genCloudBoxStruct();
                for r = 1:rCount
                    bb = rprops{r, 'BoundingBox'};
                    box = RNACloud.genCloudBoxStruct();
                    box.x0 = floor(bb(1));
                    box.y0 = floor(bb(2));
                    box.z0 = floor(bb(3));
                    box.x1 = ceil(bb(1) + bb(4) - 1);
                    box.y1 = ceil(bb(2) + bb(5) - 1);
                    box.z1 = ceil(bb(3) + bb(6) - 1);

                    cloudBoxes(r) = box;
                end
            else
                rprops = regionprops(cloudMask, 'BoundingBox');
                rCount = size(rprops, 1);
                boxCount = 0;
                for r = 1:rCount
                    rprop = rprops(r);
                    boxCount = boxCount + size(rprop.BoundingBox, 1);
                end
                
                if boxCount < 1
                    cloudBoxes = [];
                    return;
                end

                cloudBoxes(boxCount) = RNACloud.genCloudBoxStruct();
                b = 1;
                for r = 1:rCount
                    rprop = rprops(r);
                    bb = rprop.BoundingBox;
                    bc = size(bb, 1);
                    for i = 1:bc
                        box = RNACloud.genCloudBoxStruct();
                        box.x0 = floor(bb(1));
                        box.y0 = floor(bb(2));
                        box.x1 = ceil(bb(1) + bb(3) - 1);
                        box.y1 = ceil(bb(2) + bb(4) - 1);
                        box.z0 = 1;
                        box.z1 = 1;

                        cloudBoxes(b) = box;
                        b = b + 1;
                    end
                end
                clear rprops b box bb bc r rCount boxCount rprop
            end

        end
        
    end
    
end