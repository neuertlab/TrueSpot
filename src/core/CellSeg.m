%
%%
classdef CellSeg

    methods(Static)

        function param_struct = genCellSegParameterStruct()
            param_struct = struct('focus_plane_strat', 'specify');
            param_struct.min_nuc_size = 40;
            param_struct.max_nuc_size = 200;
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
        end

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

        function [trans_plane, param_struct] = GenBaseTransPlane(light_ch_data, nuc_label, param_struct)

            if ndims(light_ch_data) < 3
                trans_plane = light_ch_data;
                return;
            end
            Z = size(light_ch_data, 3);
            if param_struct.focus_offset_max > Z; param_struct.focus_offset_max = Z; end

            use_range = true;
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
                start_plane = param_struct.focus_offset_min;
                end_plane = param_struct.focus_offset_max;
            end

            if use_range
                if start_plane < param_struct.focus_offset_min; start_plane = param_struct.focus_offset_min; end
                if end_plane > param_struct.focus_offset_max; end_plane = param_struct.focus_offset_max; end
                if start_plane > end_plane; start_plane = end_plane; end
                trans_ring_planes = light_ch_data(:,:,start_plane:end_plane);
                trans_plane = max(trans_ring_planes,[],3);            % Maximum intensity z-projection

                param_struct.start_plane = start_plane;
                param_struct.end_plane = end_plane;
            end
        end

        %Derived from B1_autosegment_nuclei8_thres
        function [dapi_threshold,dapi_label] = B1_autosegment_nuclei8_thres(DAPI_ims, Yim, ManTh, min_nucleus_size, max_nucleus_size)
            a = size(DAPI_ims,1);                                                       % image size in number of pixels
            z = size(DAPI_ims,3);                                                       % stack size in number of images

            for i = 1:z;                                                                % Find the Image with the strongest DAPI signal (largest STD)
                STD2D(i) = std2(DAPI_ims(:,:,i));
            end;
            [p,ip] = max(STD2D);
            range = 3                                                                   %use the 3 images around the max image
            while ip + range > size(DAPI_ims,3) | ip - range < 1                        %reduce the range if this z stack is near the edges
                range = range - 1
            end
            dapi_max = max(DAPI_ims(:,:,ip-range:ip+range),[],3);                               % maximum intensity projection in z-direction
            Dapi_Mean = mean(dapi_max(:))
            Dapi_Min = min(dapi_max(:))
            Dapi_Max = max(dapi_max(:))
            Dapi_Median = median(dapi_max(:))
            if ManTh
                threshold_sampling = 100
            else
                threshold_sampling = 200
            end
            %% find the nuclei-maximizing dapi threshold
            if 10*Dapi_Median < Dapi_Max
                dd = round((10*Dapi_Median-Dapi_Min)/threshold_sampling)
                dapi_threshold2 = Dapi_Min:dd:10*Dapi_Median;
            else
                dd = round((Dapi_Max-Dapi_Min)/threshold_sampling)
                dapi_threshold2 = Dapi_Min:dd:Dapi_Max;
            end
            for j = 1:size(dapi_threshold2,2)
                dapi_bw = DAPI_ims > dapi_threshold2(j);                                        % only take DAPI intensities above the identified threshold for 3D stack
                dapi_bw_max2(:,:,j) = max(dapi_bw,[],3);                                            % maxium z-direction
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Ben Kesler 2/10/15 This goes through thresholds and determines how many
            % % nuclei would result from each. The threshold at which the max number of
            % % nuclei is obtained will be used
            % Check subsets of the image and
            for i = threshold_sampling/10:size(dapi_threshold2,2);                          %see how many nuclei for every threshold
                i
                dapi_threshold = dapi_threshold2(i);
                dapi_bw = DAPI_ims > dapi_threshold;                                        % only take DAPI intensities above the identified threshold for 3D stack
                dapi_bw_max = max(dapi_bw,[],3);                                           % Maximum projection for the binary pixels above the threshold
                % remove nuclei > than max_nucleus and < than min_nucleus
                dapi_normal = bwareaopen(dapi_bw_max, min_nucleus_size);                        % remove DAPI signal that are too small
                dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
                dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
                dapi_bw_max = max(dapi_bw2,[],3);                                           % Maximum projection for the binary pixels above the threshold
                dapi_OK = bwareaopen(dapi_bw_max,4);                                       % segment all DAPI spots
                dapi_label = bwlabeln(dapi_OK,8);                                          % label all segmented DAPI spots
                nuclei_num(i) = max(dapi_label(:));
            end

            ik = find(nuclei_num == max(nuclei_num),1,'last');                              %set the threshold to the largest threshold that results in the maximum cells
            dapi_threshold = dapi_threshold2(ik)

            dapi_bw = DAPI_ims > dapi_threshold;                                        % only take DAPI intensities above the identified threshold for 3D stack
            dapi_bw_max = max(dapi_bw,[],3);                                           % Maximum projection for the binary pixels above the threshold

            %% remove nuclei > than max_nucleus and < than min_nucleus
            dapi_normal = bwareaopen(dapi_bw_max, min_nucleus_size);                        % remove DAPI signal that are too small
            dapi_huge = bwareaopen(dapi_bw_max, max_nucleus_size);                          % remove DAPI signal that are too large
            dapi_bw2 = dapi_normal - dapi_huge;                                         % DAPI signal of the right size
            dapi_bw_max = max(dapi_bw2,[],3);                                           % Maximum projection for the binary pixels above the threshold

            %% Determine DAPI threshold for each individual cell
            % segment the DAPI signals in the image
            dapi_OK = bwareaopen(dapi_bw_max,4);                                       % segment all DAPI spots
            dapi_label = bwlabeln(dapi_OK,8);                                          % label all segmented DAPI spots
            m22 = max(dapi_label(:));                                                   % determine maximum nuber of DAPI stained nuclei
        end

        %Derived from B2_autosegment_cells_new &
        %B2_autosegment_cells_new_yeast
        function [Lab, cell_info, trans_plane, params] = AutosegmentCells(light_ch_data, nuc_label, params)

            %Grab dimensions
            X = size(light_ch_data, 2);
            Y = size(light_ch_data, 1);

            [trans_plane, params] = CellSeg.GenBaseTransPlane(light_ch_data, nuc_label, params);

            trans2 = imopen(trans_plane, strel('disk',10)); % 25     % smooth trans image with a disk filter of 50px to generate an background image
            trans3 = imsubtract(trans_plane,trans2);                 % subtract background for the real image to enhace signal
            trans_plane = double(trans3);                            % convert image into double format
            clear trans2 trans3

            gradmag2 = imimposemin(trans_plane, nuc_label, 26);      % generate superimposed image of the trans signal and the DAPI stained nucleus
            cells = watershed(gradmag2);                                                % segment cell using water shed algorithem
            cells_normal = bwareaopen(cells, params.min_cell_size);                     %filter out cells that are too small
            cells_huge = bwareaopen(cells, params.max_cell_size);                       %filter out cells that are too large
            cells_filtered = cells_normal - cells_huge;                                 % remaining cells have the right size
            cells = bwlabeln(cells_filtered);                                           % labels cells
            Label1 = uint16(cells);                                                     % convert segmented cell image to 16bit
            cell_count = max(Label1(:));                                                % determine the maximum number of segmented cells

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

            %Analyze cell shape for further analysis
            cell_info(mm) = struct('Centroid', NaN, 'MajorAxisLength', NaN, 'MinorAxisLength', NaN, 'FilledArea', NaN, 'Image', []);
            for j = 1 : mm                                                            % organise cell information
                % Determine cell properties
                this_cell_mask = uint16(Lab == j);  %sieve out dots in cell j
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


    end

end