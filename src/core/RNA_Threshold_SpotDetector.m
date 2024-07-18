%Master functions for running the RNA thresholding spot detection
%Blythe Hospelhorn
%Modified from code written by Ben Kesler
%Version 2.3.0
%Updated March 31, 2023

%Modified from ABs_Threshold3Dim

%Update Log:
%   1.0.0 | 19.10.24
%   2.0.0 | 22.03.24
%       (There's definitely been more changes in between there that I
%       didn't log...)
%       - Split main detect function so double image can be cleared before
%       uint16 image is used for detect.
%       - Rearranged things within main detect for better memory efficiency
%       - Removed deprecated masked bkg detect function
%   2.0.1 | 22.10.17
%       Reduced border blackout length
%       Fixed params for Gaussian
%   2.0.2 | 22.10.27
%       Removed auto z-trim. 0 by default.
%   2.1.0 | 22.11.03
%       Added handling for coord table saves > 2GB.
%   2.2.0 | 23.01.23
%       Added auto-rescaling if filtered image doesn't have sufficient
%       dynamic range.
%       Also updated so filtered image remains as double instead of uint16
%       through filtering, since data were being lost for dim or noisy images.
%   2.3.0 | 23.03.31
%       Added auto-rescaling if raw image doesn't have sufficient
%       dynamic range.

%%

%%
classdef RNA_Threshold_SpotDetector
    
    %%
    methods (Static)
        
        %%
        %Run pre-detection operations (deadpix detect, filtering),
        %returning the 16-bit filtered image. Also generates and saves max
        %projections.
        %
        %ARGS
        %   img_channel(double[Y][X][Z]) - Input image channel.
        %   save_stem (string) - Path stem for saving tables and
        %                        intermediate images.
        %   dead_pix_detect (bool) - Whether to rerun the dead pixel detection procedure (if false, looks for file with dead pixel coords)
        %   gaussian_rad (int) - Radius for Gaussian filter (optional, defaults to 7)
        %
        function [IMG_filtered] = run_spot_detection_pre(img_channel, save_stem, dead_pix_detect, gaussian_rad, save_img_proj)
            if nargin < 4; gaussian_rad = 7; end
            if nargin < 5; save_img_proj = true; end
            
            %Detect dead pixels...
            dead_pix_path = [save_stem '_deadpix.mat'];
            if (dead_pix_detect)
                RNA_Threshold_Common.saveDeadPixels(img_channel, dead_pix_path, true);
            end

            %Pre-filtering...
            IMG3D = uint16(img_channel);
            if (dead_pix_detect)
                IMG3D = RNA_Threshold_Common.cleanDeadPixels(IMG3D, dead_pix_path, true);
                if ~save_img_proj
                    if isfile(dead_pix_path)
                        delete(dead_pix_path);
                    end
                end
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

            IMG_filtered = RNA_Threshold_Common.applyGaussianFilter(IMG3D, gaussian_rad, 2);
            IMG_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(IMG_filtered);
            IMG_filtered = RNA_Threshold_Common.blackoutBorders(IMG_filtered, gaussian_rad+1, 0);
            
            %Rescale the filtered image if not enough range
            ifmin = min(IMG_filtered, [], 'all');
            ifmax = max(IMG_filtered, [], 'all');
            ifrng = ifmax - ifmin;
            if ifrng < 25
                fprintf("WARNING: Filtered image has low dynamic range. Rescale triggered!\n");
                %Linear rescale to 0-255
                IMG_filtered = ((IMG_filtered - ifmin) .* 255) ./ ifrng;
            end
            
            %Generate image structs for passing to visualization function
            %This has been moved before detect so that unfiltered image can
            %be cleared.
            max_proj = double(max(IMG3D,[],3));
            max_proj_f = double(max(IMG_filtered,[],3));
            Lmin = min(max_proj(:));
            Lmax = median(max_proj(:)) + round(10 * std(max_proj(:)));
            LminF = median(max_proj_f(:)) - round(0 * std(max_proj_f(:)));
            LmaxF = median(max_proj_f(:)) + round(10 * std(max_proj_f(:)));
            clear IMG3D;
    
            %Save
            if save_img_proj
                my_images(1).image = max_proj_f;
                my_images(1).Lmin = LminF;
                my_images(1).Lmax = LmaxF;

                my_images(2).image = max_proj;
                my_images(2).Lmin = Lmin;
                my_images(2).Lmax = Lmax;

                save([save_stem '_imgviewstructs.mat'], 'my_images');
            end
            
            IMG_filtered = uint16(IMG_filtered); %To reduce memory usage. Note that on dim images this can have a dramatic effect.
        end
        
        %%
        function [auto_zmin, auto_zmax, call_table] = run_spot_detection_main(IMG_filtered, save_stem, strategy, th_min, th_max, z_min, z_max, verbose, thread_request)
            if nargin < 6
                z_min = 0;
            end
            if nargin < 7
                z_max = 0;
            end
            if nargin < 8
                verbose = true;
            end
            if nargin < 9
                thread_request = 1;
            end

            %Generate threshold steps
            th_step = 1;
            thh = th_min:th_step:th_max;

            %auto_ztrim = ztrim;
            detect_ctx = RNADetection.generateDetectContextStruct(IMG_filtered);
            detect_ctx.th_list = thh;
            %detect_ctx.zBorder = ztrim;
            detect_ctx.minZ = z_min;
            detect_ctx.maxZ = z_max;
            detect_ctx.save_stem = save_stem;
            detect_ctx.th_strategy = strategy;
            detect_ctx.verbose = verbose;
            detect_ctx.threads = thread_request;

            detect_ctx = RNADetection.run_spotDetect(detect_ctx);
            auto_zmin = detect_ctx.minZ;
            auto_zmax = detect_ctx.maxZ;
            spot_table = detect_ctx.spot_table;
            call_table = detect_ctx.call_table;
            clear detect_ctx;
        
            [spot_table, call_table] = RNASpotsRun.saveTables(spot_table, call_table, save_stem);

            %Clean up individual thresh coord table files generated by
            %parallel pipeline?
            if thread_request > 1
                T = size(spot_table,1);
                for i = 1:T
                    th_val = spot_table(i,1);
                    cpath = [save_stem '_coords_' sprintf('%04d.mat', th_val)];
                    if isfile(cpath)
                        delete(cpath);
                    end
                end
            end
        end
        
        %%
        %Run spot detection on the given 3D image channel (prefiltered) with the
        %provided strategy and thesholds and save resulting tables to disk.
        %
        %ARGS
        %   IMG_filtered(uint16[Y][X][Z]) - Image channel prefiltered by run_spot_detection_pre
        %   save_stem (string) - Path stem for saving tables and
        %                        intermediate images.
        %   strategy (string) - Algorithm to use for spot detection
        %                           'max_proj' - Use max projection
        %                           'max_avg' - Use the plane with max avg
        %                           'all_3d' - 3D analysis (default)
        %   th_min (int) - Minimum intensity threshold value to scan
        %   th_max (int) - Maximum intensity threshold value to scan
        %   verbose (bool) - [OPTIONAL] Whether to print messages regarding
        %           processing progress
        %
        function [auto_ztrim, new_th_min] = run_spot_detection_main_dbg(IMG_filtered, save_stem, strategy, th_min, th_max, ztrim, limit2gb, verbose, thread_request)
            if nargin < 6
                ztrim = 0;
            end
            if nargin < 7
                limit2gb = true;
            end
            if nargin < 8
                verbose = true;
            end
            if nargin < 9
                thread_request = 1;
            end
            
            %Filtered image (before spot detect)
            img_filter = IMG_filtered;
            save([save_stem '_prefilteredIMG'], 'img_filter');
            clear img_filter;
            
            %Generate threshold steps
            th_step = 1;
            thh = th_min:th_step:th_max;
            %t_count = size(thh,2);
    
            %Okay, this z border needs to be a liiiiittle more dynamic.
            %Can't have a 5 slice border on images that only have 8 slices,
            %you IDIOT
%             z_trim = 5;
%             Z = size(IMG_filtered,3);
%             if(Z < 20)
%                 z_trim = 3;
%                 if (Z <= 10)
%                     z_trim = 0;
%                 end
%             end
%             auto_ztrim = z_trim;
            auto_ztrim = ztrim;
            detect_ctx = RNA_Threshold_Common.generateThreshContextStruct(IMG_filtered);
            detect_ctx.th_list = thh;
            detect_ctx.zBorder = ztrim;
            detect_ctx.save_stem = save_stem;
            detect_ctx.th_strategy = strategy;
            detect_ctx.verbose = verbose;
            detect_ctx.threads = thread_request;
            detect_ctx = RNA_Threshold_Common.run_spotDetectOnThresholdList(detect_ctx);
            spot_table = detect_ctx.spot_table;
            coord_table = detect_ctx.coord_table;
            clear detect_ctx;

            %Check coord table size...
            T = size(spot_table,1);
            [spot_table, coord_table] = RNASpotsRun.saveTables_dbg(spot_table, coord_table, save_stem, limit2gb);
            new_T = size(spot_table,1);
            new_th_min = spot_table(1,1);
            if new_T < T
                fprintf("WARNING: Coordinate table size exceeded 2GB. Lower threshold boundary was trimmed to %d.\n", new_th_min);
            end
            T = new_T;
            
            %TODO Okay, so coord table is not being saved correctly? Why not?
            %I'm so tired.

            %Clean up individual thresh coord table files generated by
            %parallel pipeline?
            if thread_request > 1
                %Check if final coord_table was saved correctly.
                mem_ctbl = coord_table;
                load([save_stem '_coordTable.mat'], 'coord_table');
                if ~isempty(coord_table) 
                    fprintf("Saved coordinate table threshold count: %d\n", size(coord_table,1));
                    if(size(coord_table,1) == T)
                        for i = 1:T
                            th_val = spot_table(i,1);
                            cpath = [save_stem '_coords_' sprintf('%04d.mat', th_val)];
                            if isfile(cpath)
                                delete(cpath);
                            end
                        end
                    end
                end
                coord_table = mem_ctbl;
            end
            
            %If 3D, collapse coord table to 2D and save
            %t_count = size(thh,2);
            cdims = size(coord_table{1},2);
            if cdims == 3
                coord_table_2D = cell(T,1);
                spot_table_2D = zeros(T,2);
                for thi = 1:T
                    spot_table_2D(thi,1) = spot_table(thi,1);
                    coord_table_2D{thi} = RNA_Threshold_Common.collapse3DCoordTable(coord_table{thi});
                    spot_table_2D(thi,2) = size(coord_table_2D{thi}, 1);
                end
                save([save_stem '_spotTable2d.mat'], 'spot_table_2D');
                save([save_stem '_coordTable2d.mat'], 'coord_table_2D');
            end
        end
        
        %%
        %Run spot detection on the given 3D image channel (raw) with the
        %provided strategy and thesholds and save resulting tables to disk.
        %
        %ARGS
        %   img_channel(double[Y][X][Z]) - Input image channel.
        %   save_stem (string) - Path stem for saving tables and
        %                        intermediate images.
        %   strategy (string) - Algorithm to use for spot detection
        %                           'max_proj' - Use max projection
        %                           'max_avg' - Use the plane with max avg
        %                           'all_3d' - 3D analysis (default)
        %   th_min (int) - Minimum intensity threshold value to scan
        %   th_max (int) - Maximum intensity threshold value to scan
        %   dead_pix_detect (bool) - Whether to rerun the dead pixel detection procedure (if false, looks for file with dead pixel coords) 
        %   verbose (bool) - [OPTIONAL] Whether to print messages regarding
        %           processing progress
        %
        function [auto_ztrim] = run_spot_detection(img_channel, save_stem, strategy, th_min, th_max, dead_pix_detect, verbose)
            
            if nargin < 7
                verbose = true;
            end
            
            [IMG_filtered] = run_spot_detection_pre(img_channel, save_stem, dead_pix_detect);
            [auto_ztrim] = run_spot_detection_main(IMG_filtered, save_stem, strategy, th_min, th_max, verbose);
    
        end
        
    end
end