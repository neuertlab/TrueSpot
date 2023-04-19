%
%%

classdef RNACoords

methods (Static)

    function call_table = convertOldCoordTable(spot_table, coord_table_old, image_filtered, image_raw, fit_incl)
        if isempty(coord_table_old)
            call_table = table.empty();
            return;
        end

        %Grab lowest threshold results (to determine alloc)
        tcoords = coord_table_old{1,1};
        alloc = size(tcoords,1);
        call_table = RNACoords.genNewCoordTable(alloc, fit_incl);

        %Fill in initial data
        is3d = (ndims(image_raw) >= 3);
        Y = size(image_raw,1);
        X = size(image_raw,2);
        if is3d
            Z = size(image_raw,3);
        else
            Z = 1;
        end
        idim = [Y X Z];
        icoords = sub2ind(idim, tcoords(:,2), tcoords(:,1), tcoords(:,3));

        call_table(:,'coord_1d') = array2table(uint32(icoords));
        call_table(:,'isnap_x') = array2table(uint16(tcoords(:,1)));
        call_table(:,'isnap_y') = array2table(uint16(tcoords(:,2)));
        call_table(:,'isnap_z') = array2table(uint16(tcoords(:,3)));
        if ~isempty(image_filtered)
            call_table(:,'intensity_f') = array2table(single(image_filtered(icoords)));
        end
        call_table(:,'intensity') = array2table(single(image_raw(icoords)));
        call_table(:,'dropout_thresh') = array2table(single(spot_table(1,1)));

        T = size(spot_table,1);
        for t = 1:T
            th = spot_table(t,1);
            tcoords = coord_table_old{t,1};
            if isempty(tcoords)
                break;
            end

            %Find in full call table.
            ticoords = sub2ind(idim, tcoords(:,2), tcoords(:,1), tcoords(:,3));
            rowidxs = find(ismember(icoords, ticoords));

            call_table(rowidxs,'dropout_thresh') = array2table(single(th));
        end
    end

    function app_table = genAppendableTable(og_call_table, alloc)
        fitlvl = 0;
        colcount = size(og_call_table,2);
        if colcount > 16; fitlvl = 1; end
        if colcount > 20; fitlvl = 2; end
        app_table = RNACoords.genNewCoordTable(alloc, fitlvl);
    end

    function call_table = genNewCoordTable(alloc, fit_incl)
        if nargin < 2
            fit_incl = 0;
        end

        varNames = {'coord_1d' 'isnap_x' 'isnap_y' 'isnap_z' 'intensity_f' 'intensity'...
            'dropout_thresh' 'is_true' 'is_trimmed_out' 'in_truth_region' 'cell'...
            'xdist_ref' 'ydist_ref' 'zdist_ref' 'xydist_ref' 'xyzdist_ref'};
        varTypes = {'uint32' 'uint16' 'uint16' 'uint16' 'single' 'single'...
            'single' 'logical' 'logical' 'logical' 'uint16'...
            'single' 'single' 'single' 'single' 'single'};

        if fit_incl >= 1
            varNames = cat(2, varNames, {'fit_x', 'fit_y', 'fit_z', 'fit_intensity'});
            varTypes = cat(2, varTypes, {'single', 'single', 'single', 'single'});
        end

        if fit_incl >= 2
            varNames = cat(2, varNames, {'xFWHM', 'yFWHM', 'fit_total_intensity'});
            varTypes = cat(2, varTypes, {'double', 'double', 'double'});
        end

        table_size = [alloc size(varNames,2)];
        call_table = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);

        call_table{:,'isnap_x'} = 0;
        call_table{:,'isnap_y'} = 0;
        call_table{:,'isnap_z'} = 0;
        call_table{:,'cell'} = 0;
        call_table{:,'coord_1d'} = 0;

        call_table{:,'intensity_f'} = NaN;
        call_table{:,'intensity'} = NaN;
        call_table{:,'dropout_thresh'} = NaN;
        call_table{:,'xdist_ref'} = NaN;
        call_table{:,'ydist_ref'} = NaN;
        call_table{:,'zdist_ref'} = NaN;
        call_table{:,'xydist_ref'} = NaN;
        call_table{:,'xyzdist_ref'} = NaN;
        if fit_incl >= 1
            call_table{:,'fit_x'} = NaN;
            call_table{:,'fit_y'} = NaN;
            call_table{:,'fit_z'} = NaN;
            call_table{:,'fit_intensity'} = NaN;
        end

        if fit_incl >= 2
            call_table{:,'xFWHM'} = NaN;
            call_table{:,'yFWHM'} = NaN;
            call_table{:,'fit_total_intensity'} = NaN;
        end

        call_table{:,'is_true'} = false;
        call_table{:,'is_trimmed_out'} = false;
        call_table{:,'in_truth_region'} = false;
    end

    function ref_assign = matchALotOfSpots(callmtx, ref_set, snaprad_3, snaprad_z, snap_minth)
        %TODO
        %For big tables
        BLOCK_SIZE = 1024;
        rcount = size(callmtx, 1);
        ccount = size(ref_set, 1);
        R = ceil(rcount/BLOCK_SIZE);
        C = ceil(ccount/BLOCK_SIZE);

        alloc = ceil(rcount * ccount * 0.001);
        combo_table = NaN(alloc,3); %dist3 row col

        %Go through blocks to get combo candidates
        rr = 1;
        for rb = 1:R
            cc = 1;
            r_end = rr + BLOCK_SIZE - 1;
            if r_end > rcount
                r_end = rcount;
            end
        
            tempmtx_rx = single(repmat(callmtx(rr:r_end,1), [1 ccount]));
            tempmtx_ry = single(repmat(callmtx(rr:r_end,2), [1 ccount]));
            tempmtx_rz = single(repmat(callmtx(rr:r_end,3), [1 ccount]));
            for cb = 1:C
                c_end = cc + BLOCK_SIZE - 1;
                if c_end > ccount
                    c_end = ccount;
                end
                tempmtx_cols = single(repmat(ref_set(cc:c_end,1).', [rcount 1]));
                xdist = tempmtx_rx - tempmtx_cols;

                tempmtx_cols = single(repmat(ref_set(cc:c_end,2).', [rcount 1]));
                ydist = tempmtx_ry - tempmtx_cols;

                tempmtx_cols = single(repmat(ref_set(cc:c_end,3).', [rcount 1]));
                zdist = tempmtx_rz - tempmtx_cols;
                %TODO

                cc = c_end + 1;
            end
            rr = r_end + 1;
        end
        clear tempmtx_rx tempmtx_ry tempmtx_rz

        %Go through combos

    end

    function ref_assign =  matchSpots(call_table, ref_set, snaprad_3, snaprad_z, snap_minth)

        %Convert callxyz to a matrix
        ref_spots = size(ref_set,1);
        call_spots = size(call_table,1);
        callmtx = zeros(call_spots,3);
        callmtx(:,1) = double(table2array(call_table(:,'isnap_x')));
        callmtx(:,2) = double(table2array(call_table(:,'isnap_y')));
        callmtx(:,3) = double(table2array(call_table(:,'isnap_z')));

        %NaN out any rows where dropout th is below minimum
        dropth = table2array(call_table(:,'dropout_thresh'));
        callmtx(dropth < snap_minth, :) = NaN;

        tempmtx_r = double(repmat(ref_set(:,1).', [call_spots 1]));
        tempmtx_c = double(repmat(callmtx(:,1), [1 ref_spots]));
        xdist = tempmtx_r - tempmtx_c;

        tempmtx_r = double(repmat(ref_set(:,2).', [call_spots 1]));
        tempmtx_c = double(repmat(callmtx(:,2), [1 ref_spots]));
        ydist = tempmtx_r - tempmtx_c;

        tempmtx_r = double(repmat(ref_set(:,3).', [call_spots 1]));
        tempmtx_c = double(repmat(callmtx(:,3), [1 ref_spots]));
        zdist = tempmtx_r - tempmtx_c;
        clear tempmtx_r tempmtx_c

        dist3 = sqrt((xdist.^2) + (ydist.^2) + (zdist.^2));
        clear xdist;
        clear ydist;
        zdist = abs(zdist);

        %Eliminate any combos above radius thresholds
        dist3(zdist > snaprad_z) = NaN;
        dist3(dist3 > snaprad_3) = NaN;
        clear zdist

        %Try to make matches
        ref_assign = zeros(1,ref_spots); %Index of call spot assigned to ref spot
        while nnz(~isnan(dist3)) > 0
            [~, minidx] = min(dist3, [], 'all', 'omitnan');
            [r, c] = ind2sub([call_spots ref_spots], minidx);
            dist3(:,c) = NaN;
            dist3(r,:) = NaN;
            ref_assign(c) = r;
        end
        clear dist3;

    end

    function [call_table, ref_assign] = updateTFCalls(call_table, ref_set, snaprad_3, snaprad_z, snap_minth)
        if nargin < 3
            snaprad_3 = 4;
        end
        if nargin < 4
            snaprad_z = 2;
        end
        if nargin < 5
            snap_minth = 25;
        end

        call_table{:,'is_true'} = false;

        ref_assign = RNACoords.matchSpots(call_table, ref_set, snaprad_3, snaprad_z, snap_minth);

        %Mark trues in table
        if nnz(ref_assign) > 0
            fidxs = find(ref_assign > 0);
            cidxs = ref_assign(fidxs);
            call_table{cidxs,'is_true'} = true;

            %Record distances from match
            temparr = double(table2array(call_table(cidxs, 'isnap_x')));
            xdist = abs(double(ref_set(fidxs,1)) - temparr);
            temparr = double(table2array(call_table(cidxs, 'isnap_y')));
            ydist = abs(double(ref_set(fidxs,2)) - temparr);
            temparr = double(table2array(call_table(cidxs, 'isnap_z')));
            zdist = abs(double(ref_set(fidxs,3)) - temparr);

            xydist = sqrt((xdist .^ 2) + (ydist .^ 2));
            xyzdist = sqrt((xdist .^ 2) + (ydist .^ 2) + (zdist .^ 2));

            call_table(cidxs, 'xdist_ref') = array2table(single(xdist));
            call_table(cidxs, 'ydist_ref') = array2table(single(ydist));
            call_table(cidxs, 'zdist_ref') = array2table(single(zdist));
            call_table(cidxs, 'xydist_ref') = array2table(single(xydist));
            call_table(cidxs, 'xyzdist_ref') = array2table(single(xyzdist));

            clear xdist ydist zdist xydist xyzdist temparr
        end

        %Import fnegs
        no_assign = (ref_assign == 0);
        fn_count = nnz(no_assign);
        if fn_count > 0
            ridxs = find(no_assign);
            tblappend = RNACoords.genAppendableTable(call_table, fn_count);
            
            x_round = double(round(ref_set(ridxs,1)));
            y_round = double(round(ref_set(ridxs,2)));
            z_round = double(round(ref_set(ridxs,3)));
            tblappend(:,'isnap_x') = array2table(uint16(x_round));
            tblappend(:,'isnap_y') = array2table(uint16(y_round));
            tblappend(:,'isnap_z') = array2table(uint16(z_round));
            tblappend{:,'dropout_thresh'} = 0;
            tblappend{:,'is_true'} = true;

            %Add distances
            xdist = abs(x_round - double(ref_set(ridxs,1)));
            ydist = abs(y_round - double(ref_set(ridxs,2)));
            zdist = abs(z_round - double(ref_set(ridxs,3)));
            xydist = sqrt((xdist .^ 2) + (ydist .^ 2));
            xyzdist = sqrt((xdist .^ 2) + (ydist .^ 2) + (zdist .^ 2));

            tblappend(:, 'xdist_ref') = array2table(single(xdist));
            tblappend(:, 'ydist_ref') = array2table(single(ydist));
            tblappend(:, 'zdist_ref') = array2table(single(zdist));
            tblappend(:, 'xydist_ref') = array2table(single(xydist));
            tblappend(:, 'xyzdist_ref') = array2table(single(xyzdist));

            call_table = [call_table;tblappend];
        end

    end

    function [call_table_merged, callmap] = mergeSlicedSetTo3D(call_table, snaprad_3, snap_minth)
        call_spots = size(call_table,1);
        callmtx = zeros(call_spots,3);
        callmtx(:,1) = double(table2array(call_table(:,'isnap_x')));
        callmtx(:,2) = double(table2array(call_table(:,'isnap_y')));
        callmtx(:,3) = double(table2array(call_table(:,'isnap_z')));
        callmap = zeros(1,call_spots); %Which index each spot is remapped to

        above_spots = zeros(1,call_spots);
        below_spots = zeros(1,call_spots);

        %NaN out any rows where dropout th is below minimum
        dropth = table2array(call_table(:,'dropout_thresh'));
        callmtx(dropth < snap_minth, :) = NaN;

        %Make a matrix of the distance of all spots from each other.
        tempmtx_r = double(repmat(callmtx(:,1).', [call_spots 1])); %Columns repeated to rows
        tempmtx_c = double(repmat(callmtx(:,1), [1 call_spots])); %Rows repeated to columns
        xdist = tempmtx_r - tempmtx_c;

        tempmtx_r = double(repmat(callmtx(:,2).', [call_spots 1]));
        tempmtx_c = double(repmat(callmtx(:,2), [1 call_spots]));
        ydist = tempmtx_r - tempmtx_c;

        tempmtx_r = double(repmat(callmtx(:,3).', [call_spots 1]));
        tempmtx_c = double(repmat(callmtx(:,3), [1 call_spots]));
        zdist = tempmtx_r - tempmtx_c;

        dist3 = sqrt((xdist.^2) + (ydist.^2) + (zdist.^2));
        clear xdist;
        clear ydist;

        %NaN out the diagonal
        dist3(find(eye(call_spots))) = NaN;

        %NaN out all but combos that are sufficiently close in 3D AND on
        %adjacent (not the same) z planes
        mask = (dist3 <= snaprad_3);
        mask = and(mask, or((zdist == 1), (zdist == -1)));
        dist3(~mask) = NaN;

        %Go through remaining combos and mark best candidates above and
        %below
        while nnz(~isnan(dist3)) > 0
            [~, minidx] = min(dist3, [], 'all', 'omitnan');
            [r, c] = ind2sub([call_spots call_spots], minidx);
            
            if zdist(r,c) > 0
                %Col is above row
                if (below_spots(c) == 0) & (above_spots(r) == 0)
                    %Neither is already paired elsewhere
                    below_spots(c) = r;
                    above_spots(r) = c;
                end
            else
                %Row is above col
                if (below_spots(r) == 0) & (above_spots(c) == 0)
                    below_spots(r) = c;
                    above_spots(c) = r;
                end
            end

            %NaN out this combo.
            dist3(r,c) = NaN;
            dist3(c,r) = NaN;

        end
        clear dist3;
        clear zdist;

        %Chain together and make assignments
        realloc = 0;
        for i = 1:call_spots
            if callmap(i) ~= 0
                %Already assigned
                continue;
            end

            if (below_spots(i) == 0) & (above_spots(i) == 0)
                %No matches
                callmap(i) = i;
                realloc = realloc + 1;
            else
                j = i;
                above_count = 0;
                while above_spots(j) ~= 0
                    above_count = above_count + 1;
                    j = above_spots(j);
                end

                j = i;
                below_count = 0;
                while below_spots(j) ~= 0
                    below_count = below_count + 1;
                    j = below_count(j);
                end

                %Find middle to merge other spots in chain into.
                nuci = i; %If equal, use this spot
                ztotal = above_count + below_count + 1;
                mid_slice = ceil((ztotal + 1) / 2);
                pos_slice = below_count + 1;
                if pos_slice > mid_slice
                    %This spot is on the high side
                    %Walk down
                    dist = pos_slice - mid_slice;
                    for k = 1:dist
                        nuci = below_spots(nuci);
                    end
                elseif pos_slice < mid_slice
                    %This spot on the low side
                    %Walk up
                    dist = mid_slice - pos_slice;
                    for k = 1:dist
                        nuci = above_spots(nuci);
                    end
                end

                %Assign nuci to all spots in chain
                j = i;
                callmap(i) = nuci;
                while above_spots(j) ~= 0
                    j = above_spots(j);
                    callmap(j) = nuci;
                end

                j = i;
                while below_spots(j) ~= 0
                    j = below_spots(j);
                    callmap(j) = nuci;
                end

                realloc = realloc + 1;
            end
        end

        %Create merged table
        call_table_merged = RNACoords.genAppendableTable(call_table, realloc);
        tblpos = 1;
        for i = 1:call_spots
            if callmap(i) ~= i
                continue;
            end

            call_table_merged(tblpos, :) = call_table(i, :); %Not sure if I can do that.
            tblpos = tblpos + 1;
        end

    end

    function call_table = addFitData(call_table, fit_table, snap_minth)
        %Fit table should be formatted with cols:
        %   x,y,z,int(optional)

        %Allow minimum th (for bigfish)
        if nargin < 3; snap_minth = 0; end

        snaprad_3 = 2;
        snaprad_z = 1;
        ftable_cols = size(fit_table,2);

        %TODO if fit table is too big, the memory required is absurd.
        %Maybe try assigning directly in order?
        call_count = size(call_table,1);
        fit_count = size(fit_table, 1);
        if call_count == fit_count
            fit_assign = [1:call_count];
        else
            fit_assign = RNACoords.matchSpots(call_table, fit_table, snaprad_3, snaprad_z, snap_minth);
        end

        if nnz(fit_assign) > 0
            fidxs = find(fit_assign > 0);
            cidxs = fit_assign(fidxs);
            call_table(cidxs, 'fit_x') = array2table(single(fit_table(fidxs,1)));
            call_table(cidxs, 'fit_y') = array2table(single(fit_table(fidxs,2)));
            call_table(cidxs, 'fit_z') = array2table(single(fit_table(fidxs,3)));
            if ftable_cols >= 4
                 call_table(cidxs, 'fit_intensity') = array2table(single(fit_table(fidxs,4)));
            end
        end

    end

    function call_table = addFitDataFromQuant(call_table, cell_rna_data)
        %TODO
    end

    function call_table = updateRefDistancesToUseFits(call_table, ref_table, ref_call_map)
        if nnz(ref_call_map) < 1; return; end
        fidxs = find(ref_call_map > 0);
        cidxs = ref_call_map(fidxs);

        call_x = double(table2array(call_table(cidxs,'fit_x')));
        call_y = double(table2array(call_table(cidxs,'fit_y')));
        call_z = double(table2array(call_table(cidxs,'fit_z')));
        nofit = find(isnan(call_x));
        
        if nnz(nofit) > 0
            i_x = double(table2array(call_table(cidxs,'isnap_x')));
            i_y = double(table2array(call_table(cidxs,'isnap_y')));
            i_z = double(table2array(call_table(cidxs,'isnap_z')));
            call_x(nofit) = i_x(nofit);
            call_y(nofit) = i_y(nofit);
            call_z(nofit) = i_z(nofit);
            clear i_x i_y i_z
        end

        ref_x = double(ref_table(fidxs,1));
        ref_y = double(ref_table(fidxs,2));
        ref_z = double(ref_table(fidxs,3));

        xdist = abs(ref_x - call_x);
        ydist = abs(ref_y - call_y);
        zdist = abs(ref_z - call_z);

        xydist = sqrt((xdist .^ 2) + (ydist .^ 2));
        xyzdist = sqrt((xdist .^ 2) + (ydist .^ 2) + (zdist .^ 2));

        call_table(cidxs, 'xdist_ref') = array2table(single(xdist));
        call_table(cidxs, 'ydist_ref') = array2table(single(ydist));
        call_table(cidxs, 'zdist_ref') = array2table(single(zdist));
        call_table(cidxs, 'xydist_ref') = array2table(single(xydist));
        call_table(cidxs, 'xyzdist_ref') = array2table(single(xyzdist));
    end

    function call_table = applyCellSegMask(call_table, cell_mask)
        x = table2array(call_table(:,'isnap_x'));
        y = table2array(call_table(:,'isnap_y'));

        if ndims(cell_mask) >= 3
            z = table2array(call_table(:,'isnap_z'));
            cells = cell_mask(x,y,z);
        else
            cells = cell_mask(x,y);
        end

        call_table(:,'cell') = array2table(uint16(cells));
    end

end

end