%
%%

classdef RNADetection
    methods (Static)

        %%
        function call_table = tempCalls2Table(temp_calls, img_filtered)
            %Cols: 1Dcoord, highest thresh passed

            varNames = {'coord_1d' 'isnap_x' 'isnap_y' 'isnap_z' 'intensity_f' 'intensity'...
                'dropout_thresh'};
            varTypes = {'uint32' 'uint16' 'uint16' 'uint16' 'single' 'single'...
                'single'};

            alloc = size(temp_calls, 1);
            table_size = [alloc size(varNames,2)];
            call_table = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);
            call_table{:,'coord_1d'} = uint32(temp_calls(:,1));
            call_table{:,'dropout_thresh'} = single(temp_calls(:,2));
            call_table{:,'intensity'} = NaN;

            call_table{:,'intensity_f'} = single(img_filtered(temp_calls(:,1)));

            [call_table{:,'isnap_y'}, call_table{:,'isnap_x'}, call_table{:,'isnap_z'}]...
                = ind2sub(size(img_filtered), call_table{:,'coord_1d'});

        end
        
        %%
        function callTable2csv(csv_path, call_table, th_min, zero_coords)
            %First, filter table to desired fields
            call_table = removevars(call_table, {'coord_1d'});
            call_table = renamevars(call_table,{'isnap_x' 'isnap_y' 'isnap_z'},{'x' 'y' 'z'});

            keep_idx = find(call_table{:,'dropout_thresh'} >= th_min);
            if isempty(keep_idx); return; end
            call_table = call_table(keep_idx,:);

            if zero_coords
                call_table{:,'x'} = call_table{:,'x'} - 1;
                call_table{:,'y'} = call_table{:,'y'} - 1;
                call_table{:,'z'} = call_table{:,'z'} - 1;
            end

            writetable(call_table, csv_path);
        end

        %%
        function callSet = condense3DCallResults(input_cell_vec, spot_count)
            Z = size(input_cell_vec,2);
            callSet = NaN(1, spot_count);
            tblpos = 1;
            for z = 1:Z
                ztbl = input_cell_vec{z};
                if isempty(ztbl); continue; end
                zcount = size(ztbl,1);
                endi = tblpos + zcount - 1;
                callSet(tblpos:endi) = ztbl(:,1);
                tblpos = endi + 1;
            end
        end

        %%
        function callSet = collapse3DCallResults(input_cell_vec, spot_count, idims)
            callSet3 = RNADetection.condense3DCallResults(input_cell_vec, spot_count);
            if isempty(callSet3)
                callSet = [];
                return;
            end

            [yy, xx, ~] = ind2sub(callSet3, idims);
            xy_table(:,2) = yy;
            xy_table(:,1) = xx;
            xy_table = unique(xy_table, 'rows');
            if isempty(xy_table)
                callSet = [];
                return;
            end
            callSet = sub2ind(idims(1:2), xy_table(:,2).', xy_table(:,1).');
        end

        %%
        function [img_filtered, calls, minVal, maxVal] = testThreshold_slice(in_slice, th1, th2)
            %Initialize outputs
            minVal = 0;
            maxVal = 0;

            %Initial filter
            IM = immultiply(in_slice, in_slice > th1); % Filters out pixels below threshold

            if sum(IM(:)) > 0 %Are there any spots at all?
                IM2 = imregionalmax(IM,8); %Filters down to only brightest pix in 2D for each spot (so can spot count)
                IM3 = immultiply(in_slice, IM2);
                minVal = min(IM3(:));
                maxVal = max(IM3(:))./2;

                calls = find(IM3 > round(th2));
                img_filtered = IM3;
            else
                calls = [];
                img_filtered = IM; %Nothing to see here
            end

        end

        %%
        function [img_filtered, calls, minVal, maxVal, spotsFound] = testThreshold_3D(in_img, th1, th2, zBorder, pretrimmed)

            if nargin < 5
                pretrimmed = false;
            end
            spotsFound = 0;

            %Initialize outputs
            %fprintf("Preparing spot testing...\t")
            %tic
            X = size(in_img, 2);
            Y = size(in_img, 1);
            Z = size(in_img, 3);

            calls = cell(1,Z); %Just to keep it easier for now.
            minVal = NaN(Z);
            maxVal = NaN(Z);

            %Dump unwanted slices
            minZ = zBorder+1;
            maxZ = Z-zBorder;
            %img_filtered = in_img(:,:,minZ:maxZ);
            img_filtered = in_img;
            if ~pretrimmed
                for z = 1:(minZ-1)
                    img_filtered(:,:,z) = 0;
                end
                for z = (maxZ+1):Z
                    img_filtered(:,:,z) = 0;
                end
            end
            %toc

            %Initial filter
            %fprintf("Filter 1\t")
            %tic
            IM = immultiply(img_filtered, img_filtered > th1); % Filters out pixels at or below threshold
            clear img_filtered; %Not using for now.
            %toc

            if sum(IM(:)) > 0 %Are there any spots at all?
                %fprintf("Filter 2\t")
                %tic
                IM2 = imregionalmax(IM,26); %Filters down to only brightest pix in each spot in 3D (so can spot count)
                clear IM;
                %toc
                %size(IM2)
                %size(in_img)

                %fprintf("Filter 3\t")
                %tic
                IM3 = immultiply(in_img, IM2);
                clear IM2; %It's not much, but maybe it'll help the memory spikes?
                %clear in_img; %Hopefully only deletes local copy of it?
                %toc

                %fprintf("Finding points...\t")
                %tic
                img_filtered = uint16(zeros(Y,X,Z));
                plane_size = X * Y;
                for z = minZ:maxZ
                    %IM3 = immultiply(in_img(:,:,z), IM2);
                    %idx = z - minZ + 1;
                    minVal(z) = min(IM3(:,:,z), [], [1, 2]);
                    maxVal(z) = max(IM3(:,:,z), [], [1, 2])./2;
                    call_coords = find(IM3(:,:,z) > round(th2(z)));
                    calls{z} = call_coords + ((z-1) * plane_size);
                    img_filtered(:,:,z) = IM3(:,:,z);
                    spotsFound = spotsFound + size(call_coords,1);
                end
                %toc
            else
                img_filtered = IM; %Nothing to see here
            end
        end

        %%
        function [img_filtered, calls, minVal, maxVal] = testThreshold_maxZ(in_img, threshold)

            %Prepare slice
            max_proj = max(double(in_img),[],3); % A 2D projection of the brightest pixels in each Z column
            p_std = std2(max_proj);

            %Run slice
            [img_filtered, calls, minVal, maxVal] = RNADetection.testThreshold_slice(max_proj, threshold, p_std);

        end

        %%
        function [img_filtered, calls, minVal, maxVal] = testThreshold_bestAvgZ(in_img, threshold, zBorder)
            %Little more than a wrapper of testThreshold_selectZ that recalculates
            %   the averages every time it is called.
            %TBH probably shouldn't use.
            Z = size(in_img, 3);
            plane_avgs = NaN(Z);

            minZ = zBorder+1;
            maxZ = Z-zBorder;

            for z = minZ:maxZ
                plane_avgs(z) = mean2(in_img(:,:,z));
            end

            [~,I] = nanmax(plane_avgs(:));
            [img_filtered, calls, minVal, maxVal] = RNADetection.testThreshold_selectZ(in_img, threshold, plane_avgs, I);

        end

        %%
        function [img_filtered, calls, minVal, maxVal] = testThreshold_selectZ(in_img, threshold, img_avg, use_planes)

            use_z = use_planes(1);
            fprintf("Using plane: %d\n", use_z)
            [img_filtered, calls, minVal, maxVal] = RNADetection.testThreshold_slice(in_img(:,:,use_z), threshold, img_avg(use_z));

        end

        %%
        function common_ctx = generateDetectContextStruct(img_filter)
            common_ctx = struct('img_filter', img_filter);
            %common_ctx.coord_table = [];
            common_ctx.call_table = table.empty();
            common_ctx.spot_table = [];
            common_ctx.plane_avgs = [];
            common_ctx.th_list = [];
            common_ctx.zBorder = 0;
            common_ctx.minZ = 0;
            common_ctx.maxZ = 0;
            common_ctx.save_stem = '';
            common_ctx.th_strategy = '';
            common_ctx.save_filtered = false;
            common_ctx.collapse3D = false;
            common_ctx.verbose = false;
            common_ctx.threads = 1;
            common_ctx.fimg_max_val = 0;
            common_ctx.valbin = [];
            common_ctx.img_size = [0 0 0]; %For quick usage
            common_ctx.temp_calls = []; %Cols: 1Dcoord, highest thresh passed
        end

        %%
        function common_ctx = do3DThresholdLoop_serial(common_ctx, th_idx)
            th = common_ctx.th_list(1,th_idx);
            if common_ctx.verbose
                tic;
                fprintf("Processing image using threshold = %f\n", th);
            end

            check_th = th+1;
            if ((check_th <= size(common_ctx.valbin,2)) & (common_ctx.valbin(1,check_th) ~= 0)) | (th_idx == 1)
                %Spot detect
                [f_img, calls, ~, ~, spot_count] = RNADetection.testThreshold_3D(common_ctx.img_filter, th, ...
                    common_ctx.plane_avgs, common_ctx.zBorder, true);

                if spot_count > 0

                    %Whittle down spots
                    if common_ctx.collapse3D
                        if common_ctx.verbose; fprintf("Determining unique XY spots...\n"); end
                        call_table = RNADetection.collapse3DCallResults(calls, spot_count, common_ctx.img_size);
                        spot_count = size(call_table, 2);
                        common_ctx.spot_table(th_idx,2) = spot_count;
                    else
                        call_table = RNADetection.condense3DCallResults(calls, spot_count);
                        common_ctx.spot_table(th_idx,2) = spot_count;
                        clear calls;
                    end

                    %Find in temp table
                    if spot_count > 1 %check again in case of collapse
                        match_mtx = ismember(common_ctx.temp_calls(:,1), call_table);
                        match_idx = find(match_mtx);
                        common_ctx.temp_calls(match_idx,2) = th;
                    end

                    if(common_ctx.save_filtered); save([common_ctx.save_stem '_sfimg_t' num2str(th)], 'f_img'); end
                    clear f_img;
                end
            else
                spot_count = common_ctx.spot_table(th_idx-1,2);
                common_ctx.spot_table(th_idx,2) = spot_count;
                if spot_count > 0
                    match_mtx = (common_ctx.temp_calls(:,2) == (th - 1));
                    match_idx = find(match_mtx);
                    common_ctx.temp_calls(match_idx,2) = th;
                end
            end
            if common_ctx.verbose; toc; end
        end

        %%
        function initParallel(threads, workers, workdir)
            mkdir(workdir);
            local_cluster = parcluster();
            local_cluster.NumThreads = threads;
            local_cluster.NumWorkers = workers;
            local_cluster.JobStorageLocation = workdir;
            parpool(local_cluster);
        end

        %%
        function do3DThresholdLoop_parallel(common_ctx, th_val)
            my_temp_file = [common_ctx.save_stem '_coords_' sprintf('%04d', th_val) '.mat'];
            if isfile(my_temp_file)
                if common_ctx.verbose
                    fprintf("Calls for threshold %d already found! Skipping...\n", th_val);
                end
                return;
            end

            %Spot detect
            [f_img, calls, ~, ~, spot_count] = RNADetection.testThreshold_3D(common_ctx.img_filter,...
                th_val, common_ctx.plane_avgs, common_ctx.zBorder, true);
            %common_ctx.img_filter = []; See if keeping it read-only helps?

            if(common_ctx.save_filtered); save([common_ctx.save_stem '_sfimg_t' num2str(th_val)], 'f_img'); end
            clear f_img;

            if spot_count > 0

                if common_ctx.collapse3D
                    t_coords = RNADetection.collapse3DCallResults(calls, spot_count, common_ctx.img_size);
                    spot_count = size(call_table, 2);
                else
                    t_coords = RNADetection.condense3DCallResults(calls, spot_count);
                end
                clear calls;

                if spot_count > 0
                    save(my_temp_file, 't_coords', '-v7.3');
                end

            end
        end

        %%
        function common_ctx = runThresholdList3D(common_ctx)

            T = size(common_ctx.th_list,2);

            %Do lowest threshold to establish base table.
            if common_ctx.verbose
                fprintf("Processing lowest threshold...\n");
            end

            th_val = common_ctx.th_list(1,1);
            [~, calls, ~, ~, spot_count] = RNADetection.testThreshold_3D(common_ctx.img_filter,...
                th_val, common_ctx.plane_avgs, common_ctx.zBorder, true);
            if common_ctx.collapse3D
                base_coords = RNADetection.collapse3DCallResults(calls, spot_count, common_ctx.img_size);
            else
                base_coords = RNADetection.condense3DCallResults(calls, spot_count);
            end
            common_ctx.spot_table(1,2) = spot_count;
            clear calls;

            if common_ctx.threads > 1
                ctx_clean.img_filter = common_ctx.img_filter;
                ctx_clean.plane_avgs = common_ctx.plane_avgs;
                ctx_clean.zBorder = common_ctx.zBorder;
                ctx_clean.minZ = common_ctx.minZ;
                ctx_clean.maxZ = common_ctx.maxZ;
                ctx_clean.save_stem = common_ctx.save_stem;
                ctx_clean.collapse3D = common_ctx.collapse3D;
                ctx_clean.save_filtered = common_ctx.save_filtered;
                ctx_clean.img_size = common_ctx.img_size;
                ctx_clean.verbose = common_ctx.verbose;
                tlist = common_ctx.th_list;

                if common_ctx.verbose
                    fprintf("Now testing thresholds using %d workers...\n", common_ctx.threads);
                    tic;
                end

                %parpool(common_ctx.threads);
                [outdir, ~, ~] = fileparts(common_ctx.save_stem);
                pardir = [outdir filesep 'parallel'];
                RNADetection.initParallel(1, common_ctx.threads, pardir);
                parfor c = 2:T
                    RNADetection.do3DThresholdLoop_parallel(ctx_clean, tlist(1,c));
                end
                delete(gcp('nocreate'));
                rmres = rmdir(pardir, 's');
                if common_ctx.verbose; toc; end

                %Rejoin into tables...
                common_ctx.temp_calls = NaN(spot_count, 2);
                common_ctx.temp_calls(:,1) = base_coords(:);
                common_ctx.temp_calls(:,2) = common_ctx.th_list(1,1);

                for c = 2:T
                    th_val = tlist(1,c);
                    partemp_path = [common_ctx.save_stem '_coords_' sprintf('%04d', th_val) '.mat'];
                    if isfile(partemp_path)
                        load(partemp_path, 't_coords');
                        match_mtx = ismember(common_ctx.temp_calls(:,1), t_coords);
                        match_idx = find(match_mtx);
                        common_ctx.temp_calls(match_idx,2) = th_val;
                        common_ctx.spot_table(c,2) = size(t_coords,2);
                    else
                        common_ctx.spot_table(c,2) = 0;
                    end
                end
            else
                common_ctx.temp_calls = NaN(spot_count, 2);
                common_ctx.temp_calls(:,1) = base_coords(:);
                common_ctx.temp_calls(:,2) = common_ctx.th_list(1,1);

                %Skip threshes with no spots at that intensity.
                common_ctx.fimg_max_val = max(common_ctx.img_filter, [], 'all', 'omitnan');
                common_ctx.valbin = histcounts(common_ctx.img_filter, common_ctx.fimg_max_val+1);

                %fprintf("breakpoint\n");
                for c = 2:T
                    common_ctx = RNADetection.do3DThresholdLoop_serial(common_ctx, c);
                end
            end
        end

        %%
        function common_ctx = run_spotDetect(common_ctx)

            %Save the basics
            Z = size(common_ctx.img_filter, 3);
            T = size(common_ctx.th_list, 2);

            if common_ctx.minZ < 1
                common_ctx.minZ = common_ctx.zBorder+1;
            end
            if common_ctx.maxZ < 1
                common_ctx.maxZ = Z-common_ctx.zBorder;
            end

            if common_ctx.maxZ > Z
                common_ctx.maxZ = Z;
            end

            minZ = common_ctx.minZ;
            maxZ = common_ctx.maxZ;

            %Pre-allocate some vars we gonna use later...
            common_ctx.plane_avgs = NaN(1,Z);
            for z = minZ:maxZ
                %plane_avgs(z) = mean2(img_filter(:,:,z));
                common_ctx.plane_avgs(z) = RNA_Threshold_Common.mean_noZeros(common_ctx.img_filter(:,:,z));
            end

            if minZ > 1
                common_ctx.img_filter(:,:,1:(minZ-1)) = 0;
            end
            if maxZ < Z
                common_ctx.img_filter(:,:,(maxZ+1):Z) = 0;
            end

            common_ctx.spot_table = NaN(T,2); %Columns are threshold, spot count. Rows are entries.

            %Move this list population up here
            common_ctx.spot_table(:,1) = common_ctx.th_list(:);

            %Generate coord & spot tables (do spot detection)
            if strcmp(common_ctx.th_strategy, 'max_proj')
                if common_ctx.verbose
                    fprintf("Running on max projection!\n");
                end
                common_ctx.maxViewMinZ = common_ctx.minZ;
                common_ctx.maxViewMaxZ = common_ctx.maxZ;
                common_ctx.minZ = 1;
                common_ctx.maxZ = 1;
                common_ctx.zBorder = 0;

                %Use only the max projection slice for spot detection
                for c = 1:T
                    th = common_ctx.th_list(1,c);
                    if common_ctx.verbose
                        fprintf("Testing threshold %d...\n", th);
                        tic;
                    end

                    %Spot detect
                    [f_img, calls, ~, ~] = RNADetection.testThreshold_maxZ(common_ctx.img_filter, th);
                    %Spot count
                    if ~isempty(calls)
                        spot_count = size(calls,1);
                        common_ctx.spot_table(c,2) = spot_count;
                        %Save coordinates
                        if isempty(common_ctx.temp_calls)
                            %Alloc and copy
                            common_ctx.temp_calls = zeros(spot_count, 2);
                            common_ctx.temp_calls(:,1) = calls(:);
                        end
                        %Update max ths
                        matches = ismember(common_ctx.temp_calls(:,1), calls);
                        common_ctx.temp_calls(matches,2) = th;
                    else
                        common_ctx.spot_table(c,2) = 0;
                    end

                    if common_ctx.save_filtered
                        savepath = sprintf("%s%s%s", common_ctx.save_stem, '_t', num2str(c));
                        save(savepath, 'f_img');
                    end
                    if common_ctx.verbose
                        toc;
                    end
                end
            else
                if strcmp(common_ctx.th_strategy, 'max_avg')
                    %Use only the slice w/ highest avg intensity for spot detection
                    [~,I] = nanmax(common_ctx.plane_avgs(:));
                    common_ctx.usedPlane = I;
                    common_ctx.minZ = 1;
                    common_ctx.maxZ = 1;
                    common_ctx.zBorder = 0;
                    for c = 1:T
                        th = common_ctx.th_list(1,c);
                        %Spot detect
                        [f_img, calls, ~, ~] = RNADetection.testThreshold_selectZ(common_ctx.img_filter, th, common_ctx.plane_avgs, I);
                        if ~isempty(calls)
                            %Spot count
                            spot_count = size(calls,1);
                            common_ctx.spot_table(c,2) = spot_count;
                            %Save coordinates
                            if isempty(common_ctx.temp_calls)
                                %Alloc and copy
                                common_ctx.temp_calls = zeros(spot_count, 2);
                                common_ctx.temp_calls(:,1) = calls(:);
                            end
                            %Update max ths
                            matches = ismember(common_ctx.temp_calls(:,1), calls);
                            common_ctx.temp_calls(matches,2) = th;
                        else
                            common_ctx.spot_table(c,2) = 0;
                        end
                    end
                else
                    %All 3D! The slow one!
                    common_ctx = RNADetection.runThresholdList3D(common_ctx);
                end
            end

            %Convert temp table.
            common_ctx.call_table = RNADetection.tempCalls2Table(common_ctx.temp_calls, common_ctx.img_filter);
            common_ctx.temp_calls = [];

            %Correct intensity values for 2D projections
            if strcmp(common_ctx.th_strategy, 'max_proj')
                max_proj = max(double(common_ctx.img_filter),[],3);
                common_ctx.call_table{:,'intensity_f'} = single(max_proj(common_ctx.call_table{:,'coord_1d'}));
                clear max_proj;
            end
            if strcmp(common_ctx.th_strategy, 'max_avg')
                my_slice = common_ctx.img_filter(:,:,common_ctx.usedPlane);
                common_ctx.call_table{:,'intensity_f'} = single(my_slice(common_ctx.call_table{:,'coord_1d'}));
                clear my_slice;
            end

        end
    end
end