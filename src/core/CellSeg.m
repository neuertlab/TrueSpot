%
%%
classdef CellSeg

    methods(Static)

        %%
        function param_struct = genCellSegParameterStruct()
            param_struct = struct('focus_plane_strat', 'specify');
            param_struct.min_cell_size = 600;
            param_struct.max_cell_size = 1200;
            param_struct.min_plane = 1;
            param_struct.max_plane = inf;
            param_struct.focus_offset_min = 3;
            param_struct.focus_offset_max = 7;
            param_struct.start_plane = 0;
            param_struct.end_plane = 0;
            param_struct.x_trim = 4;
            param_struct.y_trim = 4;

            %These are hard trims - slices outside these are ignored from
            %the start.
            param_struct.z_min = 1;
            param_struct.z_max = 0;
        end

        %%
        function param_struct = genNucSegStruct()
            param_struct = struct();

            %Params
            param_struct.range = 3;
            param_struct.threshold_sampling = 200;
            param_struct.min_nucleus_size = 40;
            param_struct.max_nucleus_size = 200;
            param_struct.cutoff = 0.05;
            param_struct.dxy = NaN;
            param_struct.use_adding = true;
            param_struct.x_trim = 4;
            param_struct.y_trim = 4;
            param_struct.max_proj_z_trim_lo = 0.2;
            param_struct.max_proj_z_trim_hi = 0.0;
            param_struct.z_min = 1;
            param_struct.z_max = 0;
        end

        %%
        function res_struct = genNucSegResultsStruct()
            res_struct = struct();

            res_struct.counter = 0;
            res_struct.nuc_threshold = NaN;
            res_struct.nuc_label = [];
            res_struct.nuc_label_lo = [];
            res_struct.lbl_lo = [];
            res_struct.lbl_mid = [];
            res_struct.lbl_hi = [];
            res_struct.nuclei = [];
            res_struct.nuc_int = [];
            res_struct.nuc_vol = [];
            res_struct.nuc_axis_major = [];
            res_struct.nuc_axis_minor = [];
            res_struct.nuclei_num = [];
            res_struct.test_sum = [];
        end

        %%
        %Derived from C1_find_dapi_threshold
        function best_threshold = FindDAPIThreshold(dapi_max, strict_min_size, strict_max_size, step_size, step_count, min_dapi)
            %I had to dig back to 2014 to find these functions, so I will
            %ASSUME they are deprecated. - Blythe

            %{
                Runs through a series of thresholds on the dapi max projection, counting
                total objects.  The threshold that results in the maximum number of nuclei
                is returned.
            %}

            %[BH] These default to yeast 100X values
            if nargin < 2; strict_min_size = 40; end
            if nargin < 3; strict_max_size = 200; end
            if nargin < 4; step_size = 20; end
            if nargin < 5; step_count = 150; end
            if nargin < 5; min_dapi = 500; end

            %don't want to use mean/median because in images with very few cells, the
            %background intensity can be the false threshold.  Set this value 1000-2000
            %below minimum  nuclear DAPI signal

            best_threshold = 0;
            total_nuclei = 0;

            %% goes through a range of thresholds
            max_dapi = min_dapi + (step_size * step_count);
            for threshold = min_dapi:step_size:max_dapi
                bw = dapi_max > threshold;                                                % generate a binary image (0 or 1) if the pixel intensity is above the threshold
                bwlab = bwlabeln(bw);                                                   % Label connected components in binary image
                nuc_sizes = regionprops(bwlab,'Area');                                  % Measure properties of image regions such as the Area

                %only considering objects of reasonable size ; good size = difference between all objects above the minimum size substracted by objects above the maximum size
                num_good_size = ...
                    sum(cat(1,nuc_sizes.Area) > strict_min_size)- ...                % Objects larger then the minimum size
                    sum(cat(1,nuc_sizes.Area) > strict_max_size);                    % Objects larger then the maxium size

                if (total_nuclei<=num_good_size)                                        % if the total nuclei number is smaller or equal the number of nuclei with the right size;
                    best_threshold = threshold;                                         % update the threshold otherwise do nothing
                    total_nuclei = num_good_size;                                       % update the number of nuclei
                end

            end
        end

        %%
        %Derived from C2_find_trans_bestplane
        function best_plane = PickFocusPlane_MostCells(light_ch_data, nuclei_max, strict_min_size, strict_max_size)
            %{
            Runs through all stacks, segments, and finds the plane
            that gives the most cells.
            %}

            if nargin < 3; strict_min_size = 600; end
            if nargin < 4; strict_max_size = 1200; end

            best_plane = 1;                                                             % start value for the best plane
            total_cells = 0;                                                              % total number of cells;

            Z = 1;
            if ndims(light_ch_data) > 2; Z = size(light_ch_data,3); end
            for plane = 1:Z                                            % go through each trans image

                gradmag2 = imimposemin(light_ch_data(:,:,plane), nuclei_max);           % overlay transimage with nucleus to generate cells with nuclear signal
                L = watershed(gradmag2);                                            % segment cells
                cells = bwlabeln(L);                                                %

                cell_sizes = regionprops(cells,'Area');                             % determine area for each cell

                %only considering cells of reasonable size
                num_good_size = ...                                                 % determine cells that are larger then the "strict_min_size" and smaller then the "strict_max_size"
                    sum(cat(1,cell_sizes.Area) > strict_min_size)- ...
                    sum(cat(1,cell_sizes.Area) > strict_max_size);

                if (total_cells <= num_good_size)                                         % if "total_cells" are smaller then the number of cells that within the size range
                    total_cells = num_good_size;
                    best_plane = plane;                                                 % choose the plane of the cells with the right size
                end
            end
        end

        %%
        %Derived from C3_find_trans_midplane
        function mid_plane = PickFocusPlane_MidPlane(light_ch_data)
            %{
            Runs through all stacks and identifies all extended maxima.
            The plane with the least total maxima area is returned.

            The idea is that the cell-center trans plane has the 
            least contrast relative to the background and blends in
            the most.  Without the perimiter rings (maxima), this
            plane should be easy to identify with the
            imextendedmax function.
            %}

            mid_plane = 0;                                                              % start value for the mid plane
            min_area = inf;                                                               % start value for the minimum area

            Z = 1;
            if ndims(light_ch_data) > 2; Z = size(light_ch_data,3); end
            for plane = 1:Z                                           % go through each trans image
                I_sc = mat2gray(light_ch_data(:,:,plane));                              % convert into binary image
                mask_em = imextendedmax(I_sc, .1);                                  % generate mask
                if min_area > sum(mask_em(:))                                         % if minarea is larger then umber of objects in the image
                    min_area = sum(mask_em(:));                                       % update min area
                    mid_plane = plane;                                                % update mid plane
                end
            end
        end

        %%
        %Derived from C4_find_trans_plane
        function mid_plane = PickFocusPlane_LowStdev(light_ch_data)
            Z = 1;
            if ndims(light_ch_data) > 2; Z = size(light_ch_data,3); end

            transSTD = NaN(1,Z);
            for i = 1:Z                                               % go through each image in the stack
                ims = light_ch_data(:,:,i);
                transSTD(i) = std2(ims(:));                                             % compute the standard deviation (STD) in each image;
            end
            [~,mid_plane] = min(transSTD, [], 'omitnan');                                         % the image with the lowest contrast has the lowest STD

            %The original version applies the LoG filter and returns that
            %too.
        end

        %%
        function [trans_plane, param_struct] = GenBaseTransPlane(light_ch_data, nuc_label, param_struct)

            if ndims(light_ch_data) < 3
                trans_plane = light_ch_data;
                return;
            end
            Z = size(light_ch_data, 3);
            if param_struct.focus_offset_max > Z; param_struct.focus_offset_max = Z; end

            use_range = true;
            strat = param_struct.focus_plane_strat;
            if strcmp(strat,'max_cells')
                plane = CellSeg.PickFocusPlane_MostCells(light_ch_data,...
                    nuc_label, param_struct.min_cell_size, param_struct.max_cell_size);
                trans_plane = light_ch_data(:,:,plane);
                use_range = false;
            elseif strcmp(strat,'midplane')
                plane = CellSeg.PickFocusPlane_MidPlane(light_ch_data);
                start_plane = plane + param_struct.focus_offset_min;
                end_plane = plane + param_struct.focus_offset_max;
            elseif strcmp(strat,'midplane2')
                plane = CellSeg.PickFocusPlane_LowStdev(light_ch_data);
                start_plane = plane + param_struct.focus_offset_min;
                end_plane = plane + param_struct.focus_offset_max;
            elseif strcmp(strat,'first5')
                start_plane = 1;
                end_plane = 5;
            elseif strcmp(strat,'last5')
                start_plane = Z-4;
                end_plane = Z;
            else
                start_plane = param_struct.min_plane;
                end_plane = param_struct.max_plane;
            end

            if ~isfinite(start_plane); start_plane = 1; end
            if ~isfinite(end_plane); end_plane = Z; end
            if use_range
                if start_plane < param_struct.min_plane; start_plane = param_struct.min_plane; end
                if end_plane > param_struct.max_plane; end_plane = param_struct.max_plane; end
                if start_plane < 1; start_plane = 1; end
                if end_plane > Z; end_plane = Z; end
                if start_plane > end_plane; start_plane = end_plane; end
                trans_ring_planes = light_ch_data(:,:,start_plane:end_plane);
                trans_plane = max(trans_ring_planes,[],3);            % Maximum intensity z-projection

                param_struct.start_plane = start_plane;
                param_struct.end_plane = end_plane;
            end
        end

        %%
        %Derived from B1_autosegment_nuclei8_thres
        function [nucSegSpecs, nucSegRes] = InitSegNucleiThres(nuc_ch_data, nucSegSpecs, nucSegRes)
%             if ManTh
%                 threshold_sampling = 100
%             else
%                 threshold_sampling = 200
%             end
            
            Z = size(nuc_ch_data, 3);
            
            STD2D = NaN(1,Z);
            for z = 1:Z                                                               % Find the Image with the strongest DAPI signal (largest STD)
                STD2D(z) = std2(nuc_ch_data(:,:,z));
            end

            [~, ip] = max(STD2D);                                                            %use the 3 images around the max image
            range = nucSegSpecs.range;
            while (ip + range > Z) | (ip - range < 1)                        %reduce the range if this z stack is near the edges
                range = range - 1;
            end

            nuc_max_proj = max(nuc_ch_data(:,:,ip-range:ip+range),[],3);                               % maximum intensity projection in z-direction
            nuc_min = min(nuc_max_proj(:));
            nuc_max = max(nuc_max_proj(:));
            nuc_median = median(nuc_max_proj(:));

            % find the nuclei-maximizing dapi threshold
            threshold_sampling = nucSegSpecs.threshold_sampling;
            if (10 * nuc_median) < nuc_max
                dd = round((10 * nuc_median - nuc_min) / threshold_sampling);
                test_thresh = nuc_min:dd:(10 * nuc_median);
            else
                dd = round((nuc_max - nuc_min) / threshold_sampling);
                test_thresh = nuc_min:dd:nuc_max;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Ben Kesler 2/10/15 This goes through thresholds and determines how many
            % % nuclei would result from each. The threshold at which the max number of
            % % nuclei is obtained will be used
            % Check subsets of the image and
            T = size(test_thresh, 2);
            nucSegRes.nuclei_num = zeros(1,T);
            min_nucleus_size = nucSegSpecs.min_nucleus_size;
            max_nucleus_size = nucSegSpecs.max_nucleus_size;
            for i = (threshold_sampling / 10):T                    %see how many nuclei for every threshold
                nucSegRes.nuc_threshold = test_thresh(i);
                dapi_bw = DAPI_ims > nucSegRes.nuc_threshold;                                        % only take DAPI intensities above the identified threshold for 3D stack
                dapi_bw_max = max(dapi_bw,[],3);                                           % Maximum projection for the binary pixels above the threshold
                dapi_normal = bwareaopen(dapi_bw_max, min_nucleus_size);                        % remove DAPI signal that are too small
                dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
                dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
                dapi_bw_max = max(dapi_bw2,[],3);                                           % Maximum projection for the binary pixels above the threshold
                dapi_OK = bwareaopen(dapi_bw_max,4);                                       % segment all DAPI spots
                nucSegRes.nuc_label = bwlabeln(dapi_OK,8);                                          % label all segmented DAPI spots
                nucSegRes.nuclei_num(i) = max(nucSegRes.nuc_label(:));
            end
            clear test_thresh

            ik = find(nucSegRes.nuclei_num == max(nucSegRes.nuclei_num), 1, 'last');                              %set the threshold to the largest threshold that results in the maximum cells
            nucSegRes.nuc_threshold = dapi_threshold2(ik);

            dapi_bw = DAPI_ims > nucSegRes.nuc_threshold;                                        % only take DAPI intensities above the identified threshold for 3D stack
            dapi_bw_max = max(dapi_bw, [], 3);                % Maximum projection for the binary pixels above the threshold

            % remove nuclei > than max_nucleus and < than min_nucleus
            dapi_normal = bwareaopen(dapi_bw_max, min_nucleus_size);                        % remove DAPI signal that are too small
            dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
            dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
            dapi_bw_max = max(dapi_bw2,[],3);                                           % Maximum projection for the binary pixels above the threshold

            % Determine DAPI threshold for each individual cell
            % segment the DAPI signals in the image
            dapi_OK = bwareaopen(dapi_bw_max,4);                                       % segment all DAPI spots
            nucSegRes.nuc_label = bwlabeln(dapi_OK,8);                                          % label all segmented DAPI spots                                                % determine maximum nuber of DAPI stained nuclei
        end

        %%
        %Derived from B1_autosegment_nuclei8_adding
        function [nucSegSpecs, nucSegRes] = InitSegNucleiAdd(nuc_ch_data, nucSegSpecs, nucSegRes)
            % This segments the cells by adding together binary images at every threshold, then removing pixels that are not present in most, and then going to each cell to segment it further.
            X = size(nuc_ch_data,2);                                                       % image size in number of pixels
            Y = size(nuc_ch_data,1);
            Z = size(nuc_ch_data,3);                                                       % stack size in number of images

            %Z Trim
            if(nucSegSpecs.z_max < 1); nucSegSpecs.z_max = Z; end

%             nuc_ch_data_ztrim = nuc_ch_data;
%             if(nucSegSpecs.z_min > 1)
%                 nuc_ch_data_ztrim(:,:,1:(nucSegSpecs.z_min - 1)) = NaN;
%             end
%             if(nucSegSpecs.z_max < 1); nucSegSpecs.z_max = Z; end
%             if(nucSegSpecs.z_max < Z)
%                 nuc_ch_data_ztrim(:,:,(nucSegSpecs.z_max + 1):Z) = NaN;
%             end

            STD2D = NaN(1,Z);                                                           %BK 5/16, needed this so timepoints with more stacks wouldn't be left over
            for z = 1:Z                                                                % Find the Image with the strongest DAPI signal (largest STD)
                STD2D(z) = std2(nuc_ch_data(:,:,z));
            end
            clear z
            if(nucSegSpecs.z_min > 1)
                STD2D(1:(nucSegSpecs.z_min - 1)) = NaN;
            end
            if(nucSegSpecs.z_max < Z)
                STD2D((nucSegSpecs.z_max + 1):Z) = NaN;
            end

            [~,ip] = max(STD2D);

            range = nucSegSpecs.range;
            while (ip + range > Z) | (ip - range < 1)                        %reduce the range if this z stack is near the edges
                range = range - 1;
            end
            z_min = ip - range;
            z_max = ip + range;
            nuc_max_proj = max(nuc_ch_data(:,:,z_min:z_max), [], 3);                            % maximum intensity projection in z-direction

            nuc_min = min(nuc_max_proj(:));
            nuc_max = max(nuc_max_proj(:));
            nuc_median = median(nuc_max_proj(:));
            clear nuc_max_proj ip

%             if false;%Ywin
%                 threshold_sampling = 100
%             else
%                 threshold_sampling =200
%             end

            % find the nuclei-maximizing dapi threshold
            threshold_sampling = nucSegSpecs.threshold_sampling;
            med10 = 10 * nuc_median;
            if med10 < nuc_max                                         %BK 4/27/2016
                dd = round((med10 - nuc_min) / threshold_sampling);
                test_thresh = nuc_min:dd:med10;
            else
                dd = round((nuc_max - nuc_min) / threshold_sampling);
                test_thresh = nuc_min:dd:nuc_max;
            end
            clear nuc_min nuc_max nuc_median med10 dd

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Ben Kesler 2/10/15 This goes through thresholds and determines how many
            % % nuclei would result from each. The threshold at which the max number of
            % % nuclei is obtained will be used
            T = size(test_thresh, 2);
            nucSegRes.nuclei_num = zeros(1,T);  %added 9/7 BK
            nucSegRes.test_sum = zeros(Y,X); %added 4/29 BK
            min_nucleus_size = nucSegSpecs.min_nucleus_size;
            max_nucleus_size = nucSegSpecs.max_nucleus_size;
            min_slice = round(nucSegSpecs.max_proj_z_trim_lo * Z) + 1;
            max_slice = Z - round(nucSegSpecs.max_proj_z_trim_hi * Z);
            if min_slice < 1; min_slice = 1; end
            if max_slice > Z; max_slice = Z; end

%              figure(1);
%              clf;
            for i = (threshold_sampling/10):T                        %see how many nuclei for every threshold
                nucSegRes.nuc_threshold = test_thresh(i);
                dapi_bw = (nuc_ch_data > nucSegRes.nuc_threshold); % only take DAPI intensities above the identified threshold for 3D stack
                dapi_bw_max = max(dapi_bw(:,:,min_slice:max_slice), [], 3, 'omitnan');  % Maximum projection for the binary pixels above the threshold
%                  if i == T
%                      imshow(dapi_bw_max, []);
%                  end

                % remove nuclei > than max_nucleus and < than min_nucleus
                dapi_normal = bwareaopen(dapi_bw_max, min_nucleus_size);                        % remove DAPI signal that are too small
                dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
                dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
                dapi_bw_max = max(dapi_bw2, [], 3);                                           % Maximum projection for the binary pixels above the threshold
                dapi_OK = bwareaopen(dapi_bw_max, 4);                                       % segment all DAPI spots
                nucSegRes.nuc_label = bwlabeln(dapi_OK, 8);                                          % label all segmented DAPI spots
                nucSegRes.nuclei_num(i) = max(nucSegRes.nuc_label(:));

%                  if i == T
%                      imshow(dapi_normal, []);
%                      imshow(dapi_huge, []);
%                      imshow(dapi_bw2, []);
%                  end

                %added 4-29 BK
                nucSegRes.test_sum = nucSegRes.test_sum + dapi_OK;
            end
            clear test_thresh dapi_bw min_slice dapi_OK...
                dapi_bw_max dapi_bw_max

%             figure(1);
%             clf;
%             imshow(nucSegRes.test_sum, []);

            % Apply cutoff for added image
            cutoff = nucSegSpecs.cutoff;    %This is a kind of cutoff. It's a proportion of the maximum number of times a pixel is present
            DAPI_ims_final = nucSegRes.test_sum;
            if cutoff > 0
                DAPI_ims_final(find(DAPI_ims_final < (max(DAPI_ims_final(:)) * cutoff))) = 0;
                dapi_normal = bwareaopen(DAPI_ims_final, min_nucleus_size);                        % remove DAPI signal that are too small
                dapi_huge = bwareaopen(DAPI_ims_final, max_nucleus_size);                          % remove DAPI signal that are too large
                DAPI_ims_cut = dapi_normal - dapi_huge;
                DAPI_ims_final = immultiply(DAPI_ims_cut, DAPI_ims_final);
%                 imshow(DAPI_ims_cut, []);
                clear DAPI_ims_cut
            end

            nucSegRes.nuc_label = bwlabeln(DAPI_ims_final, 8);
            max_label = max(nucSegRes.nuc_label, [], 'all', 'omitnan');
            xs = zeros(size(nucSegRes.nuc_label(:)));
            ys = zeros(size(nucSegRes.nuc_label(:)));
            for j = 1:max_label
                %collect data about each spot. Find the center (not adjusted by intensity right now)
                [y_es, x_es] = find(nucSegRes.nuc_label == j);
                xs(j,1) = mean(x_es);
                ys(j,1) = mean(y_es);
            end
            clear x_es y_es

            % Separate individual DAPI spots
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dxy2 = round(nucSegSpecs.dxy * 2.5);                                               % determine maximum number of DAPI stained nuclei                                                              %counter for how many times the box had more than one nucleus
            nucSegRes.nuc_label_lo = DAPI_ims_final;
            for i = 1:max_label
                k1 = (nucSegRes.nuc_label == i);                                                  % Single cell per image
                k2a = immultiply(k1, DAPI_ims_final);

                xA = round(xs(i));                                                      %use the center found earlier as the starting point for drawing the box
                yA = round(ys(i));
                x0 = xA - dxy2;                                                            % generate corner points with are dxy1 pixels away from the nuclear center
                x1 = xA + dxy2;
                y0 = yA - dxy2;
                y1 = yA + dxy2;

                if x0 < 1; x0 = 1; end                                                           % if corner A is below 1 pixel or negative set equals 1
                if x1 > X; x1 = X; end                                                               % if corner B is above a pixels (1024 or 2048) set equals a
                if y0 < 1; y0 = 1; end                                                           % if corner C is below 1 pixel or negative set equals 1
                if y1 > Y; y1 = Y; end                                                             % if corner D is above a pixels (1024 or 2048) set equals a

                sqim = k2a(y0:y1,x0:x1);                                         % makes an image that consists only of the box around the nucleus                                          %label the image
                thresmat = 1:max(k2a(:));                                                % The array that will have the cell numb for diff thresholds
                %if there is more or less than one cell in the box, don't do the next steps
                for y = 1:max(k2a(:))
                    sqim_temp = sqim > 0;
                    sqim_temp(find(sqim < y)) = 0;
                    dapi_normal = bwareaopen(sqim_temp, min_nucleus_size);                        % remove DAPI signal that are too small
                    dapi_huge = bwareaopen(sqim_temp, max_nucleus_size);                          % remove DAPI signal that are too large
                    sqim_temp = dapi_normal - dapi_huge;
                    sqim_lbl1 = bwlabeln(sqim_temp, 8);
                    thresmat(2, y) = max(sqim_lbl1(:));
                end

                thres = thresmat(1, find(thresmat(2,:) == max(thresmat(2,:)), 1, 'last'));  %find the threshold that results in the most number of nuclei (maximum threhsold)
                thres2 = thresmat(1, find(thresmat(2,:) == max(thresmat(2,:)), 1, 'first'));  %find the threshold that results in the most number of nuclei (minimum threshold)
                k2a(k2a < thres) = 0;                                        %Change the scaled image so it is at the new threshold
                k2a_bin = k2a > 0;                                                   %Make binary of new scaled image
                subtr_im = k1 - k2a_bin;                                            %Determine pixels that need to be subtracted from scaled image
                DAPI_ims_final(find(subtr_im > 0)) = 0;                                     %Change final image to subtract other pixels
                k2a = immultiply(k1, nucSegRes.nuc_label_lo);                                %Reset the image
                k2a(k2a < thres2) = 0;                                        %Change the scaled image so it is at the new threshold (for minimum stacks removed)

                k2a_bin = k2a > 0;                                                   %Make binary of new scaled image
                subtr_im = k1 - k2a_bin;                                            %Determine pixels that need to be subtracted from scaled image
                nucSegRes.nuc_label_lo(find(subtr_im > 0)) = 0;                                     %Change final image to subtract other pixels
            end

            dapi_normal = bwareaopen(DAPI_ims_final, min_nucleus_size);                        % remove DAPI signal that are too small
            dapi_huge = bwareaopen(DAPI_ims_final, max_nucleus_size);                          % remove DAPI signal that are too large
            DAPI_ims_final = dapi_normal - dapi_huge;
            
%             figure(2);
%             clf;
%             imshow(max(dapi_normal, [], 3, 'omitnan'), []);
% 
%             figure(3);
%             clf;
%             imshow(max(DAPI_ims_final, [], 3, 'omitnan'), []);
            
            dapi_normal = bwareaopen(nucSegRes.nuc_label_lo, min_nucleus_size);                        % remove DAPI signal that are too small
            dapi_huge = bwareaopen(nucSegRes.nuc_label_lo, max_nucleus_size);                          % remove DAPI signal that are too large
            nucSegRes.nuc_label_lo = dapi_normal - dapi_huge;

            nucSegRes.nuc_label = bwlabeln(DAPI_ims_final, 8);
            nucSegRes.nuc_label_lo = bwlabeln(nucSegRes.nuc_label_lo, 8);
            xs = zeros(size(nucSegRes.nuc_label(:)));
            ys = zeros(size(nucSegRes.nuc_label(:)));
            for j = 1:max(nucSegRes.nuc_label(:))
                [y_es, x_es] = find(nucSegRes.nuc_label == j);
                xs(j,1) = mean(x_es);
                ys(j,1) = mean(y_es);
            end

            nucSegRes.counter = nucSegRes.counter + 1;
        end

        %%
        %Derived from B1_autosegment_nuclei5
        function [nucSegSpecs, nucSegRes] = SegmentThreshedNuclei(nuc_ch_data, nucSegSpecs, nucSegRes)
            X = size(nuc_ch_data,2);                                       % image size in number of pixels
            Y = size(nuc_ch_data,1);
            Z = size(nuc_ch_data,3);                                       % stack size in number of images

            dapi_label_low1 = nucSegRes.nuc_label_lo;
            
            dxy2 = round(nucSegSpecs.dxy * 1.3);
            nuc_count = max(dapi_label_low1(:));     % determine maximum number of DAPI stained nuclei

            min_slice = round(nucSegSpecs.max_proj_z_trim_lo * Z) + 1;
            max_slice = Z - round(nucSegSpecs.max_proj_z_trim_hi * Z);
            if min_slice < 1; min_slice = 1; end
            if max_slice > Z; max_slice = Z; end
            nuc_max_proj = max(nuc_ch_data(:,:,min_slice:max_slice), [], 3, 'omitnan');
            clear min_slice max_slice

%             figure(1);
%             clf;
%             imshow(nuc_max_proj, []);

            Label_low = zeros(Y,X,Z);                                      % generate zero 3D matrix of the size of the image stack
            Label_mid = zeros(Y,X,Z);
            Label_hi = zeros(Y,X,Z);
            Label_index = zeros(Y,X);
            Nuc_int = zeros(3,1);                                                     %Will store the integrated intensity of each nucleus
            Maj_axis = zeros(3,1);                                                     %Will store the length of the major axis of each nucleus
            Min_axis = zeros(3,1);                                                     %Will store the length of the minor axis of each nucleus
            Nuc_vol = zeros(3,1);                                                     %Will store the volume of each nucleus
            
            % go through all the DAPI nucleus labels one by one
            % Thresholding by individual cells
            for i = 1:nuc_count
                k1 = (dapi_label_low1 == i);                                                  % Single cell per image

%                 figure(9);
%                 clf;
%                 imshow(max(k1, [], 3), []);

                [ys, xs] = ind2sub(size(k1), find(k1));
                sumposx = sum(xs(:)); sumposy = sum(ys(:));
                sumfluor = sum(k1(:));
                xA = round(sumposx/sumfluor);
                yA = round(sumposy/sumfluor);

                clear k1 k2 sumposx sumfluor xs ys
                Imzero = uint16(zeros(Y,X));
                Imzero(yA, xA) = 1;

                xx0 = xA - dxy2;                                                            % generate corner points with are dxy2 pixels away from the nuclear center
                xx1 = xA + dxy2;
                yy0 = yA - dxy2;
                yy1 = yA + dxy2;

                xx0 = max(xx0, 1);      % if corner A is below 1 pixel or negative set equals 1
                xx1 = min(xx1, X);      % if corner B is above a pixels (1024 or 2048) set equals a
                yy0 = max(yy0, 1);      % if corner C is below 1 pixel or negative set equals 1
                yy1 = min(yy1, Y);      % if corner D is above a pixels (1024 or 2048) set equals a

                %Fast way of generating circle
%                 area1 = zeros(Y,X);
%                 for j = 1:X                                                                 %This sets the area by making a circle around the center
%                     for k = 1:Y
%                         if sqrt((j-xA)^2+(k-yA)^2) <= dxy2
%                             area1(j,k,1) = 1;
%                         end
%                     end
%                 end
%                 clear j k

                %Generate a circle around the center
                [xx,yy] = meshgrid(1:X, 1:Y);
                xx = (xx - xA) .^ 2;
                yy = (yy - yA) .^ 2;
                dd = sqrt(xx + yy);
                area1 = (dd <= dxy2);
                clear xx yy dd

                % Determine how much area inside circle (area1) is taken up by the dapi_label_low1 image
                total_area = sum(area1(:));       %Total area inside area1
                nuc_area = immultiply(area1(:,:), (dapi_label_low1 == i));
                nuc_area = sum(nuc_area(:));
                cumul_cutoff = 1 - (nuc_area/total_area); %Tentative cutoff for the cumulative distribution
                clear total_area nuc_area xA yA

                dapi_cell = immultiply(uint16(nuc_ch_data), uint16(repmat(area1,[1,1,Z])));                         % generate dapi image stack inside area1 and set the rest of the image to zero
                other_nuc =  immultiply((dapi_label_low1 ~= i), (dapi_label_low1 > 0));     %  Find where other nuclei are
                other_nuc =  immultiply(other_nuc, area1(:,:));        %find where other nuclei are in area
                dapi_cell_max = immultiply(nuc_max_proj, area1(:,:));                                    % determine maximum projection
                temp1 = dapi_cell_max(dapi_cell_max > 0);
                dapi_cell(repmat(other_nuc,[1,1,Z]) == 1) = min(temp1(:));
                temp1 = dapi_cell_max(dapi_cell_max > 0);
                dapi_cell_max(other_nuc == 1) = min(temp1(:));                                 %Set the intensity of the other nuclei to the min
                clear area1 other_nuc temp1

%                 figure(3);
%                 clf;
%                 imshow(dapi_cell_max, []);

                %TODO [BH] If corners (xx and yy) are vectors, will this be a
                %problem?
                dapi_cell2 = dapi_cell(yy0:yy1,xx0:xx1,:);
                dapi_cell2max = double(dapi_cell_max(yy0:yy1,xx0:xx1));                                 % cut out the maxiumum projection                                % cut out the maxiumum projection (smaller for max)
                dapi_cell2max(dapi_cell2max == 0) = NaN;  %This makes all zero elements NaN instead
                clear dapi_cell dapi_cell_max

%                 figure(8);
%                 clf;
%                 imshow(dapi_cell2max, []);
% 
%                 figure(12);
%                 clf;
%                 imshow(max(dapi_cell2, [], 3), []);

                % Determine thresholds based on cutoffs in cdf
                [ysss, xsss] = ecdf(dapi_cell2max(:));                                       %determine cumulative distribution                                                  %Find values closest to certain threshold in cumulative distribution function
                z2sss = abs(ysss - cumul_cutoff);                        %threshold at cutoff for cumulative distribution function
                mm5 = xsss(find((z2sss == min(z2sss)), 1, 'first'));                          %threshold at cutoff for cumulative distribution function                      %threshold at cutoff for cumulative distribution function
                mm4 = mm5 * 0.9;
                mm6 = mm5 * 1.1;

%                 figure(11);
%                 clf;
%                 plot(xsss, ysss);

                clear cumul_cutoff ysss xsss z2sss

                m1 = mm4;                                                       % set "LOW" threshold to first cumulative distribution cutoff
                m2 = mm5;                                                       % set "MID" threshold to second of cumulative distribution cutoff
                m3 = mm6;                                                       % set "HI" threshold to third of cumulative distribution cutoff
                filtxy_A = (dapi_cell2max > m1);
                filtxy_B = (dapi_cell2max > m2);
                filtxy_C = (dapi_cell2max > m3);
                dapi_label3A = dapi_cell2 > (m1);                                        % "LOW" generate 3D binary image above the thresholds
                dapi_label3B = dapi_cell2 > (m2);                                        % "MID" generate 3D binary image above the thresholds
                dapi_label3C = dapi_cell2 > (m3);                                        % "HI" generate 3D binary image above the thresholds
                clear m1 m2 m3 mm4 mm5 mm6 dapi_cell2max

                for z = 1:Z
                    dapi_label3A(:,:,z) = immultiply(dapi_label3A(:,:,z), filtxy_A);
                    dapi_label3B(:,:,z) = immultiply(dapi_label3B(:,:,z), filtxy_B);
                    dapi_label3C(:,:,z) = immultiply(dapi_label3C(:,:,z), filtxy_C);
                end

%                 figure(4);
%                 clf;
%                 imshow(max(dapi_label3A, [], 3), []);
% 
%                 figure(7);
%                 clf;
%                 imshow(max(dapi_label3B, [], 3), []);

                % This fills holes in the image in 3D, but it is very computationally intensive
                %%% Below fills in holes and takes the most connected area in 2D
                for k = 1:Z
                    dapi_label3A(:,:,k) = imfill(dapi_label3A(:,:,k),'holes');                                        % "LOW" fill holes in image
                    dapi_label3B(:,:,k) = imfill(dapi_label3B(:,:,k),'holes');                                        % "MID" fill holes in image
                    dapi_label3C(:,:,k) = imfill(dapi_label3C(:,:,k),'holes');                                        % "HI" fill holes in image                                 % "HI" fill holes in image
                end
                clear k

                %Find 3D connected area in 3D
                dapi_label_low = bwlabeln(dapi_label3A, 6);                                          % "LOW" label different connected neighborhoods
                dapi_label_mid = bwlabeln(dapi_label3B, 6);                                          % "MID" label different connected neighborhoods
                dapi_label_hi = bwlabeln(dapi_label3C, 6);                                          % "HI" label different connected neighborhoods
                dapi_label3A = (dapi_label_low == mode(dapi_label_low(dapi_label_low > 0)));                                        % "LOW" fill holes in image
                dapi_label3B = (dapi_label_mid == mode(dapi_label_mid(dapi_label_mid > 0)));                                        % "MID" fill holes in image
                dapi_label3C = (dapi_label_hi == mode(dapi_label_hi(dapi_label_hi > 0)));
                clear dapi_label_low dapi_label_mid dapi_label_hi

%                 figure(5);
%                 clf;
%                 imshow(max(dapi_label3A, [], 3), []);
% 
%                 figure(6);
%                 clf;
%                 imshow(max(dapi_label3B, [], 3), []);

                % Eliminate Cells at the border
                temp_props = regionprops3(dapi_label3A,'BoundingBox') ;              %Obtain axis lengths for nucleus with low threshold
                nucBounds = temp_props.BoundingBox;                                  %create the rectangular box around the cell
                shapeCount = size(nucBounds,1);
                trimOut = false;
                x_trim = nucSegSpecs.x_trim;
                y_trim = nucSegSpecs.y_trim;
                for ii = 1:shapeCount  %Go through each shape
                    X0=round(nucBounds(ii,2));
                    Y0=round(nucBounds(ii,1));
                    X1=round(nucBounds(ii,2)+ nucBounds(ii,5))-1;
                    Y1=round(nucBounds(ii,1)+ nucBounds(ii,4))-1;
                    if (X0+xx0 < x_trim) | (X1+xx0 > (X-x_trim)) | (Y0+yy0 < y_trim) | (Y1+yy0 > (Y-y_trim))
                        trimOut = true;
                    end
                end
                clear x_trim y_trim nucBounds shapeCount ii temp_props

                if ~trimOut
                    % Store aspects of nucleus
                    Nuc_int(1,i) = sum(sum(sum(immultiply(dapi_label3A,dapi_cell2))));          %Store integrated intensity of nucleus with lower threshold
                    Nuc_int(2,i) = sum(sum(sum(immultiply(dapi_label3B,dapi_cell2))));          %Store integrated intensity of nucleus with middle threshold
                    Nuc_int(3,i) = sum(sum(sum(immultiply(dapi_label3C,dapi_cell2))));          %Store integrated intensity of nucleus with lower threshold
                    Nuc_vol(1,i) = sum(dapi_label3A(:));                                        %Store volume of nucleus with low threshold
                    Nuc_vol(2,i) = sum(dapi_label3B(:));                                        %Store volume of nucleus with medium threshold
                    Nuc_vol(3,i) = sum(dapi_label3C(:));                                        %Store volume of nucleus with high threshold

                    Label_low(yy0:yy1,xx0:xx1,:) = Label_low(yy0:yy1,xx0:xx1,:) + dapi_label3A;                                   % Add to label matrix 'LOW'
                    Label_mid(yy0:yy1,xx0:xx1,:) = Label_mid(yy0:yy1,xx0:xx1,:) + dapi_label3B;                                   % Add to label matrix 'MID'
                    Label_hi(yy0:yy1,xx0:xx1,:) = Label_hi(yy0:yy1,xx0:xx1,:) + dapi_label3C;                                     % Add to label matrix 'HI'

%                     figure(10);
%                     clf;
%                     imshow(max(Label_low, [], 3), []);

                    Label_index = uint16(Label_index) + Imzero;
                else
                    %Too close to image edge, remove it.
                    nucSegRes.nuc_label(dapi_label_low1 == i) = 0;
                end
                clear trimOut Imzero xx0 xx1 yy0 yy1 
                clear dapi_cell2 dapi_label3A dapi_label3B dapi_label3C
            end

            nucSegRes.lbl_lo = Label_low > 0;                                   % Make sure it is binary image
            nucSegRes.lbl_mid = Label_mid > 0;                                   % Make sure it is binary image
            nucSegRes.lbl_hi = Label_hi > 0;                                     % Make sure it is binary image

%             figure(2);
%             clf;
%             imshow(max(nucSegRes.lbl_mid, [], 3), []);

            nucSegRes.nuclei = Label_index;                                                       % max projection of the segmented dapi signals using the 50% threshold
        
            %Save everything else to output struct
            nucSegRes.nuc_int = Nuc_int;
            nucSegRes.nuc_vol = Nuc_vol;
            nucSegRes.nuc_axis_major = Maj_axis;
            nucSegRes.nuc_axis_minor = Min_axis;
        end

        %%
        %Derived from B2_autosegment_cells_new &
        %B2_autosegment_cells_new_yeast
        function [Lab, cell_info, trans_plane, params] = AutosegmentCells(light_ch_data, nuc_label, params)

            %Grab dimensions
            X = size(light_ch_data, 2);
            Y = size(light_ch_data, 1);
            Z = size(light_ch_data, 3);

            nucmax = uint16(max(nuc_label, [], 3, 'omitnan'));

%             figure(6);
%             clf;
%             imshow(nucmax, []);

            %Z trim
            if(params.z_min > 1)
                light_ch_data(:,:,1:(params.z_min - 1)) = NaN;
                if(params.min_plane < params.z_min)
                    params.min_plane = params.z_min;
                end
            end
            if(params.z_max < 1); params.z_max = Z; end
            if(params.z_max < Z)
                light_ch_data(:,:,(params.z_max+1):Z) = NaN;
                if(params.max_plane > params.z_max)
                    params.max_plane = params.z_max;
                end
            end
            if(params.max_plane < 1)
                params.max_plane = params.z_max;
            end

            [trans_plane, params] = CellSeg.GenBaseTransPlane(light_ch_data, nuc_label, params);

            trans2 = imopen(trans_plane, strel('disk',10)); % 25     % smooth trans image with a disk filter of 50px to generate an background image
            trans3 = imsubtract(trans_plane,trans2);                 % subtract background for the real image to enhace signal
            trans_plane = double(trans3);                            % convert image into double format
            clear trans2 trans3

            gradmag2 = imimposemin(trans_plane, nuc_label, 26);      % generate superimposed image of the trans signal and the DAPI stained nucleus
            cells = watershed(gradmag2);                                                % segment cell using water shed algorithem
           
%             figure(1);
%             clf;
%             imshow(gradmag2, []);
% 
%             figure(2);
%             clf;
%             imshow(cells, []);

            cells_normal = bwareaopen(cells, params.min_cell_size);                     %filter out cells that are too small
            cells_huge = bwareaopen(cells, params.max_cell_size);                       %filter out cells that are too large
            cells_filtered = cells_normal - cells_huge;                                 % remaining cells have the right size
            cells = bwlabeln(cells_filtered);                                           % labels cells
            Label1 = uint16(cells);                                                     % convert segmented cell image to 16bit
            cell_count = max(Label1(:));                                                % determine the maximum number of segmented cells
            clear cells_filtered

%             figure(3);
%             clf;
%             imshow(Label1, []);

            %BH addition. Clean up using nuclei info
%             nonly = and(Label1 == 0, nucmax);
%             if nnz(nonly) >= 500
%                 %Merging into other cells. Need some kind of outlining to
%                 %put back in.
%                 nonly = imdilate(nonly, strel('disk',10));
%                 newcells = or(Label1 ~= 0, nonly);
%                 cells = bwlabeln(newcells);
%                 Label1 = uint16(cells);
%                 cell_count = max(Label1(:));
% 
% %                 figure(7);
% %                 clf;
% %                 imshow(nonly, []);
% 
%                 clear newcells
%             end

%             figure(4);
%             clf;
%             imshow(Label1, []);

            clear cells cells_normal cells_huge

            %Go through each cell and analyze the shape, remove cells and get information
            for j = 1:cell_count                                                                % m number of cells
                this_cell_mask = uint16(Label1 == j);                                                       % Single cell per image
                cell_box = regionprops(this_cell_mask,'BoundingBox','Area','Centroid');                     % cut out rectange with cell / increases post processing
                box_bounds = cell_box.BoundingBox;                                                    %create the rectangular box around the cell
                box_area = cell_box.Area;
                X0 = round(box_bounds(1));
                Y0 = round(box_bounds(2));
                X1 = round(box_bounds(1)+ box_bounds(3))-1;
                Y1 = round(box_bounds(2)+ box_bounds(4))-1;

                %remove cells that are too small or big
                if (box_area > params.max_cell_size) | (sum(box_area < params.min_cell_size))
                    Label1 = Label1 - (this_cell_mask * j);
                end

                %remove cells on the border of the image
                if (X0 < params.x_trim) | (X1 > (X - params.x_trim)) |...
                        (Y0 < params.y_trim) | (Y1 > (Y - params.y_trim))
                    Label1 = Label1 - (this_cell_mask * j);
                end

                clear this_cell_mask cell_box box_bounds box_area
            end

%             figure(5);
%             clf;
%             imshow(Label1, []);

            %generate cell label
            %[BH] This labelf appears to be for debug purposes, so I'll
            %   leave it.
            %Labelf = zeros(Y,X);
            Label1 = Label1 > 0;                                                        % generate a binary image of the segmented cells
            Label1 = imfill(Label1, 'holes');                                           % closes holes in the cells
            %Labelf = imadd(double(Labelf), double(Label1));                             % generate segmented cells with orginial intenity profile
            %Lab = Labelf > 0;                                        % make binary image
            [Lab, mm] = bwlabel(Label1, 8);                                                % index cells
            Lab = uint16(Lab);                                                          % convert image to 16bit format to reduce image size and save memory

            if (mm < 1)
                cell_info = [];
                return;
            end

            %Analyze cell shape for further analysis
            cell_info(mm) = struct('Centroid', NaN, 'MajorAxisLength', NaN, 'MinorAxisLength', NaN, 'FilledArea', NaN, 'Image', []);
            for j = 1 : mm                                                            % organise cell information
                this_cell_mask = uint16(Lab == j);  %sieve out dots in cell j
                
                % Try to cleanup edge to include full nucleus
                this_nuc_mask = immultiply(nucmax, this_cell_mask);
                %most common nonzero value
                all_nuc_nums = this_nuc_mask(:);
                all_nuc_nums = all_nuc_nums(all_nuc_nums ~= 0);
                if ~isempty(all_nuc_nums)
                    nuc_lbl_num = mode(all_nuc_nums, 'all');
                    trg_nuc_mask = uint16(nucmax == nuc_lbl_num);
                    not_in_cell = immultiply(trg_nuc_mask, ~this_nuc_mask);
                    total_nuc_pix = nnz(trg_nuc_mask);
                    outside_cell_pix = nnz(not_in_cell);
                    if outside_cell_pix > 0
                        if (outside_cell_pix <= (total_nuc_pix * 0.2))
                            Lab(outside_cell_pix) = j;
                            this_cell_mask = uint16(Lab == j);
                        end
                    end
                end

                % Determine cell properties
                try
                    cell_reg = regionprops(this_cell_mask,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea','Image');
                    cell_info(j).Centroid = cell_reg.Centroid;
                    cell_info(j).MajorAxisLength = cell_reg.MajorAxisLength;
                    cell_info(j).MinorAxisLength = cell_reg.MinorAxisLength;
                    cell_info(j).FilledArea = cell_reg.FilledArea;
                    cell_info(j).Image = cell_reg.Image;
                catch
                end
                clear this_cell_mask cell_reg
            end

        end

        %%
        function [nucSegSpecs, nucSegRes] = AutosegmentNuclei(nuc_ch_data, nucSegSpecs)
            if isempty(nuc_ch_data); return; end

            %Z Trim
%             Z = size(nuc_ch_data, 3);
%             if(nucSegSpecs.z_min > 1)
%                 nuc_ch_data(:,:,1:(nucSegSpecs.z_min - 1)) = NaN;
%             end
%             if(nucSegSpecs.z_max < 1); nucSegSpecs.z_max = Z; end
%             if(nucSegSpecs.z_max < Z)
%                 nuc_ch_data(:,:,(nucSegSpecs.z_max + 1):Z) = NaN;
%             end

            if isnan(nucSegSpecs.dxy)
                %From A0_segment_define_variables_streamlined_non_GUI
                nucSegSpecs.dxy = round(sqrt(nucSegSpecs.max_nucleus_size/pi)); 
            end
            nucSegRes = CellSeg.genNucSegResultsStruct();

            if nucSegSpecs.use_adding
                [nucSegSpecs, nucSegRes] = CellSeg.InitSegNucleiAdd(nuc_ch_data, nucSegSpecs, nucSegRes);
                nucSegRes.nuc_threshold = NaN;
            else
                [nucSegSpecs, nucSegRes] = CellSeg.InitSegNucleiThres(nuc_ch_data, nucSegSpecs, nucSegRes);
            end

            [nucSegSpecs, nucSegRes] = CellSeg.SegmentThreshedNuclei(nuc_ch_data, nucSegSpecs, nucSegRes);
        end

        %%
        function cell_mask = openCellMask(path)
            %For checking input format. Can open the CellSeg outputs, tsv, csv, or tif.
            cell_mask = [];
            if endsWith(path, '.mat')
                finfo = who('-file', path);
                if ~isempty(find(ismember(finfo, 'cellSeg'),1))
                    load(path, 'cellSeg');
                    if isfield(cellSeg, 'cell_mask')
                        cell_mask = cellSeg.cell_mask;
                    end
                    clear cellSeg;
                elseif ~isempty(find(ismember(finfo, 'cells'),1))
                    load(path, 'cells');
                    cell_mask = cells;
                    clear cells;
                end
            elseif endsWith(path, '.tif') | endsWith(path, '.tiff')
                %Assumes one channel, one plane.
                [channels, ~] = LoadTif(path, 1, [1], 0);
                cell_mask = channels{1,1};
                clear channels;
            elseif endsWith(path, '.tsv') | endsWith(path, '.csv')
                cell_mask = readmatrix(path);
            end
        end

        %%
        function nuc_mask = openNucMask(path, maskno, fill3)
            if nargin < 2
                maskno = 2; %lblmid
            end
            if nargin < 3
                fill3 = true;
            end

            nuc_mask = [];
            if endsWith(path, '.mat')
                finfo = who('-file', path);
                if ~isempty(find(ismember(finfo, 'nucleiSeg'),1))
                    load(path, 'nucleiSeg');
                    if isfield(nucleiSeg, 'results')
                        switch(maskno)
                            case 0
                                if isfield(nucleiSeg.results, 'nuc_label')
                                    nuc_mask = nucleiSeg.results.nuc_label;
                                end
                            case 1
                                if isfield(nucleiSeg.results, 'lbl_lo')
                                    nuc_mask = nucleiSeg.results.lbl_lo;
                                end
                            case 2
                                if isfield(nucleiSeg.results, 'lbl_mid')
                                    nuc_mask = nucleiSeg.results.lbl_mid;
                                end
                            case 3
                                if isfield(nucleiSeg.results, 'lbl_hi')
                                    nuc_mask = nucleiSeg.results.lbl_hi;
                                end
                        end
                    end
                    if fill3 & (ndims(nuc_mask) > 2)
                        if isfield(nucleiSeg, 'params')
                            Z = size(nuc_mask, 3);
                            z_min = 1;
                            z_max = Z;
                            if isfield(nucleiSeg.params, 'z_min')
                                z_min = nucleiSeg.params.z_min;
                            end
                            if isfield(nucleiSeg.params, 'z_max')
                                z_max = nucleiSeg.params.z_max;
                            end

                            if z_min > 1
                                for zz = 1:(z_min - 1)
                                    nuc_mask(:,:,zz) = nuc_mask(:,:,z_min);
                                end
                            end

                            if z_max < Z
                                for zz = (z_max + 1):Z
                                    nuc_mask(:,:,zz) = nuc_mask(:,:,z_max);
                                end
                            end
                        end
                    end
                    clear nucleiSeg;
                elseif ~isempty(find(ismember(finfo, 'nuclei'),1))
                    switch(maskno)
                        case 0
                            load(path, 'nuclei');
                            nuc_mask = nuclei;
                            clear nuclei;
                        case 1
                            load(path, 'Label_low');
                            nuc_mask = Label_low;
                            clear Label_low;
                        case 2
                            load(path, 'Label_mid');
                            nuc_mask = Label_mid;
                            clear Label_mid;
                        case 3
                            load(path, 'Label_hi');
                            nuc_mask = Label_hi;
                            clear Label_hi;
                    end
                end
            elseif endsWith(path, '.tif') | endsWith(path, '.tiff')
                %Assumes one channel, one plane.
                [channels, ~] = LoadTif(path, 1, [1], 0);
                nuc_mask = channels{1,1};
                clear channels;
            elseif endsWith(path, '.tsv') | endsWith(path, '.csv')
                nuc_mask = readmatrix(path);
            end


            
        end

    end

end