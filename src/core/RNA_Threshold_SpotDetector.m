%Master functions for running the RNA thresholding spot detection
%Blythe Hospelhorn
%Modified from code written by Ben Kesler
%Version 1.0.0
%Updated Oct 24, 2019

%Modified from ABs_Threshold3Dim

%Update Log:
%   1.0.0 | 19.10.24

%%

%%
classdef RNA_Threshold_SpotDetector
    
    %%
    methods (Static)
        
        %%
        %Run spot detection on the given 3D image channel (raw) with the
        %provided strategy and thesholds and save resulting tables to disk.
        %
        %ARGS
        %   img_channel(num[Y][X][Z]) - Input image channel.
        %   save_stem (string) - Path stem for saving tables and
        %                        intermediate images.
        %   strategy (string) - Algorithm to use for spot detection
        %                           'max_proj' - Use max projection
        %                           'max_avg' - Use the plane with max avg
        %                           'all_3d' - 3D analysis (default)
        %   th_min (int) - Minimum intensity threshold value to scan
        %   th_max (int) - Maximum intensity threshold value to scan
        %   dead_pix_detect (bool) - Whether to rerun the dead pixel detection procedure (if false, looks for file with dead pixel coords) 
        %
        function [auto_ztrim] = run_spot_detection(img_channel, save_stem, strategy, th_min, th_max, dead_pix_detect, verbose)
            
            if nargin < 7
                verbose = true;
            end
            
            %Detect dead pixels...
            dead_pix_path = [save_stem '_deadpix'];
            if (dead_pix_detect)
                RNA_Threshold_Common.saveDeadPixels(img_channel, dead_pix_path);
            end

            %Pre-filtering...
            IMG3D = uint16(img_channel);
            IMG3D = RNA_Threshold_Common.cleanDeadPixels(IMG3D, dead_pix_path);

            fs = 13;
            IMG_filtered = RNA_Threshold_Common.applyGaussianFilter(IMG3D, fs);
            IMG_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(IMG_filtered);
            IMG_filtered = RNA_Threshold_Common.blackoutBorders(IMG_filtered, fs, 0);

            %Generate threshold steps
            %th_min = thA(1);
            %th_step = thA(2);
            %th_max = thA(3);
            th_step = 1;
            thh = th_min:th_step:th_max;
            %t_count = size(thh,2);
    
            %Okay, this z border needs to be a liiiiittle more dynamic.
            %Can't have a 5 slice border on images that only have 8 slices,
            %you IDIOT
            z_trim = 5;
            Z = size(IMG_filtered,3);
            if(Z < 20)
                z_trim = 3;
                if (Z < 10)
                    z_trim = 0;
                end
            end
            auto_ztrim = z_trim;
            [spot_table, coord_table] = RNA_Threshold_Common.run_spotDetectOnThresholdList(IMG_filtered, thh, strategy, z_trim, false, [save_stem '_sfimg'], false, verbose);

            %Spot count table
            save([save_stem '_spotTable'], 'spot_table');
            %Spot coords
            save([save_stem '_coordTable'], 'coord_table');
            %Filtered image (before spot detect)
            img_filter = IMG_filtered;
            save([save_stem '_prefilteredIMG'], 'img_filter');
            
            %If 3D, collapse coord table to 2D and save
            t_count = size(thh,2);
            cdims = size(coord_table{1},2);
            if cdims == 3
                coord_table_2D = cell(t_count,1);
                spot_table_2D = zeros(t_count, 2);
                for thi = 1:t_count
                    spot_table_2D(thi,1) = thh(thi);
                    coord_table_2D{thi} = RNA_Threshold_Common.collapse3DCoordTable(coord_table{thi});
                    spot_table_2D(thi,2) = size(coord_table_2D{thi}, 1);
                end
                save([save_stem '_spotTable2d'], 'spot_table_2D');
                save([save_stem '_coordTable2d'], 'coord_table_2D');
            end
    
            max_proj = double(max(IMG3D,[],3));
            max_proj_f = double(max(IMG_filtered,[],3));
            Lmin = min(max_proj(:));
            Lmax = median(max_proj(:)) + round(10 * std(max_proj(:)));
            LminF = median(max_proj_f(:)) - round(0 * std(max_proj_f(:)));
            LmaxF = median(max_proj_f(:)) + round(10 * std(max_proj_f(:)));
    
            %Generate image structs for passing to visualization function
            my_images(1).image = max_proj_f;
            my_images(1).Lmin = LminF;
            my_images(1).Lmax = LmaxF;
    
            my_images(2).image = max_proj;
            my_images(2).Lmin = Lmin;
            my_images(2).Lmax = Lmax;

            %Save
            save([save_stem '_imgviewstructs'], 'my_images');
        end
       
        %%
        %Run spot detection on the given 3D image channel (raw) masked with the provided mask
        % using the provided strategy and thesholds and save resulting tables to disk.
        %
        %ARGS
        %   img_channel(num[Y][X][Z]) - Input image channel.
        %   bkg_mask (bool[Y][X]) - Mask. Only pixels from source channel
        %                           that pass through the masked are analyzed.
        %   save_stem (string) - Path stem for saving tables and
        %                        intermediate images.
        %   strategy (string) - Algorithm to use for spot detection
        %                           'max_proj' - Use max projection
        %                           'max_avg' - Use the plane with max avg
        %                           'all_3d' - 3D analysis (default)
        %   th_min (int) - Minimum intensity threshold value to scan
        %   th_max (int) - Maximum intensity threshold value to scan
        %   dead_pix_detect (bool) - Whether to rerun the dead pixel detection procedure (if false, looks for file with dead pixel coords) 
        %
        function run_spot_detection_on_background(img_channel, bkg_mask, save_stem, strategy, th_min, th_max, dead_pix_detect)
            
            Z = size(img_channel, 3);
            
            %Detect dead pixels...
            if (dead_pix_detect)
                RNA_Threshold_Common.saveDeadPixels(img_channel);
            end

            %Pre-filtering...
            IMG3D = uint16(img_channel);
            IMG3D = RNA_Threshold_Common.cleanDeadPixels(IMG3D);

            fs = 13;
            IMG_filtered = RNA_Threshold_Common.applyGaussianFilter(IMG3D, fs);
            IMG_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(IMG_filtered);
            IMG_filtered = RNA_Threshold_Common.blackoutBorders(IMG_filtered, fs, 0);
            
            MIMG3D = IMG3D;
            MIMG_filtered = IMG_filtered;
            for z = 1:Z
                MIMG_filtered(:,:,z) = immultiply(IMG_filtered(:,:,z), bkg_mask);
                MIMG3D(:,:,z) = immultiply(IMG3D(:,:,z), bkg_mask);
            end

            %Generate threshold steps
            %th_min = thA(1);
            %th_step = thA(2);
            %th_max = thA(3);
            th_step = 1;
            thh = th_min:th_step:th_max;
            %t_count = size(thh,2);
    
            [spot_table, coord_table] = RNA_Threshold_Common.run_spotDetectOnThresholdList(MIMG_filtered, thh, strategy, 5, 1, [save_stem '_sfimg'], false);

            %Spot count table
            save([save_stem '_spotTable'], 'spot_table');
            %Spot coords
            save([save_stem '_coordTable'], 'coord_table');
            %Filtered image (before spot detect)
            img_filter = MIMG_filtered;
            save([save_stem '_prefilteredIMG'], 'img_filter');
    
            max_proj = double(max(MIMG3D,[],3));
            max_proj_f = double(max(MIMG_filtered,[],3));
            Lmin = min(max_proj(:));
            Lmax = median(max_proj(:)) + round(10 * RNA_Threshold_Common.std_noZeros(max_proj));
            LminF = median(max_proj_f(:)) - round(0 * RNA_Threshold_Common.std_noZeros(max_proj_f));
            LmaxF = median(max_proj_f(:)) + round(10 * RNA_Threshold_Common.std_noZeros(max_proj_f));
    
            %Generate image structs for passing to visualization function
            my_images(1).image = max_proj_f;
            my_images(1).Lmin = LminF;
            my_images(1).Lmax = LmaxF;
    
            my_images(2).image = max_proj;
            my_images(2).Lmin = Lmin;
            my_images(2).Lmax = Lmax;

            %Save
            save([save_stem '_imgviewstructs'], 'my_images');
        
        end
        
    end
end