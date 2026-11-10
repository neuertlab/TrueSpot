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
        end

        %%
        function cloudBox = genCloudBoxStruct()
            cloudBox = struct();
            cloudBox.x = 0;
            cloudBox.xSize = 0;
            cloudBox.y = 0;
            cloudBox.ySize = 0;
            cloudBox.z = 0;
            cloudBox.zSize = 0;
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

            [xx, yy] = meshgrid((1:X) - spot_x, (1:Y) - spot_y);
            s_dist = sqrt(xx.^2 + yy.^2);
            winrad = 0;

            last_avg = -1.0;
            for r = 1:settings.xyMaxDist
                %Mask
                rmask2 = (s_dist <= r);
                img_rmask2 = sliceRaw;
                img_rmask2(~rmask2) = NaN;

                %Take the average of everything inside this radius
                this_avg = mean(img_rmask2, 'all', 'omitnan');

                %Take ratio to average of everything inside last tested
                %radius
                %If it drops off (there is a dropoff in intensity), call
                %window
                if last_avg >= 0.0
                    ratio = this_avg/last_avg;
                    % if ratio > obj.d_thresh
                    %     winrad = r;
                    % end
                    if ratio < settings.avgRatioThresh
                        winrad = r;
                    end
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
            rmask2 = (s_dist <= winrad);
            img_rmask2 = sliceRaw;
            img_rmask2(~rmask2) = NaN;

            w_med = median(img_rmask2, 'all', 'omitnan');
            w_std = std(img_rmask2, 0, 'all', 'omitnan');
            c_thresh = w_med - (settings.stdevCount * w_std);
            clear img_rmask2 w_med w_std sliceRaw s_dist

            %Apply threshold in 3D to get mask
            rmask3 = repmat(rmask2, [1 1 Z]);
            cloudMask = and(rmask3, imageRaw >= c_thresh);
        end
        
    end
    
end