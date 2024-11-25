%
%%

classdef RNACoords

methods (Static)

    %%
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
        xx = min(tcoords(:,1), X);
        yy = min(tcoords(:,2), Y);
        zz = max(min(tcoords(:,3), Z),1);
        icoords = sub2ind(idim, yy, xx, zz);


        call_table(:,'coord_1d') = array2table(uint32(icoords));
        call_table(:,'isnap_x') = array2table(uint16(xx));
        call_table(:,'isnap_y') = array2table(uint16(yy));
        call_table(:,'isnap_z') = array2table(uint16(zz));
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
            xx = min(tcoords(:,1), X);
            yy = min(tcoords(:,2), Y);
            zz = max(min(tcoords(:,3), Z),1);

            %Find in full call table.
            ticoords = sub2ind(idim, yy, xx, zz);
            rowidxs = find(ismember(icoords, ticoords));

            call_table(rowidxs,'dropout_thresh') = array2table(single(th));
        end
    end

    %%
    function app_table = genAppendableTable(og_call_table, alloc)
        if isempty(og_call_table)
            app_table = table.empty();
            return;
        end
        varNames = og_call_table.Properties.VariableNames;

        %https://www.mathworks.com/matlabcentral/answers/304728-getting-data-types-of-table
        varTypes = varfun(@class,og_call_table,'OutputFormat','cell');
        table_size = [alloc size(varNames,2)];
        app_table = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);

        fieldCount = size(varTypes, 2);
        for i = 1:fieldCount
            fieldType = varTypes{i};
            if strcmp(fieldType, 'single') | strcmp(fieldType, 'double')
                app_table{:, i} = NaN;
            end
            if strcmp(fieldType, 'logical')
                app_table{:, i} = false;
            end
        end
    end
    
    %%
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
            varNames = cat(2, varNames, {'xFWHM', 'yFWHM', 'fit_total_intensity', 'zfitq'});
            varTypes = cat(2, varTypes, {'double', 'double', 'double', 'double'});
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
            call_table{:,'zfitq'} = NaN;
        end

        call_table{:,'is_true'} = false;
        call_table{:,'is_trimmed_out'} = false;
        call_table{:,'in_truth_region'} = false;
    end

    %%
    function call_table = extendHBOutputCallTable(og_call_table)
        if isempty(og_call_table)
            call_table = table.empty();
            return;
        end

        rowcount = size(og_call_table, 1);
        call_table = RNACoords.genNewCoordTable(rowcount, false);
        
        call_table{:,'coord_1d'} = og_call_table{:,'coord_1d'};
        call_table{:,'isnap_x'} = og_call_table{:,'isnap_x'};
        call_table{:,'isnap_y'} = og_call_table{:,'isnap_y'};
        call_table{:,'isnap_z'} = og_call_table{:,'isnap_z'};
        call_table{:,'intensity_f'} = og_call_table{:,'intensity_f'};
        call_table{:,'intensity'} = og_call_table{:,'intensity'};
        call_table{:,'dropout_thresh'} = og_call_table{:,'dropout_thresh'};
    end

    %%
    function cand_table = findMatchCandidates(callmtx, ref_set, snaprad_3, snaprad_z)
        BLOCK_SIZE = 1024;
        rcount = size(callmtx, 1);
        ccount = size(ref_set, 1);
        R = ceil(rcount/BLOCK_SIZE);
        C = ceil(ccount/BLOCK_SIZE);

        alloc = ceil(rcount * ccount * 0.001);
        combo_table = single(NaN(alloc,5)); %dist3 row col dist2 distz
        table_pos = 1;

        %Go through blocks to get combo candidates
        rr = 1;
        for rb = 1:R
            cc = 1;
            r_end = rr + BLOCK_SIZE - 1;
            if r_end > rcount
                r_end = rcount;
            end
        
            tempmtx_rx = single(repmat(callmtx(rr:r_end,1), [1 BLOCK_SIZE]));
            tempmtx_ry = single(repmat(callmtx(rr:r_end,2), [1 BLOCK_SIZE]));
            tempmtx_rz = single(repmat(callmtx(rr:r_end,3), [1 BLOCK_SIZE]));
            rb_size = r_end - rr + 1;
            for cb = 1:C
                c_end = cc + BLOCK_SIZE - 1;
                %cb_size = BLOCK_SIZE;
                if c_end > ccount
                    c_end = ccount;
                    cb_size = c_end - cc + 1;
                    tempmtx_rx = single(repmat(callmtx(rr:r_end,1), [1 cb_size]));
                    tempmtx_ry = single(repmat(callmtx(rr:r_end,2), [1 cb_size]));
                    tempmtx_rz = single(repmat(callmtx(rr:r_end,3), [1 cb_size]));
                end
                tempmtx_cols = single(repmat(ref_set(cc:c_end,1).', [rb_size 1]));
                xdist = tempmtx_rx - tempmtx_cols;

                tempmtx_cols = single(repmat(ref_set(cc:c_end,2).', [rb_size 1]));
                ydist = tempmtx_ry - tempmtx_cols;

                tempmtx_cols = single(repmat(ref_set(cc:c_end,3).', [rb_size 1]));
                zdist = tempmtx_rz - tempmtx_cols;

                dist2 = sqrt((xdist.^2) + (ydist.^2));
                dist3 = sqrt((xdist.^2) + (ydist.^2) + (zdist.^2));
                clear xdist;
                clear ydist;
                zdist = abs(zdist);

                %Eliminate any combos above radius thresholds
                dist3(zdist > snaprad_z) = NaN;
                dist3(dist3 > snaprad_3) = NaN;

                %Save remaining combos to the table.
                okayidx = find(~isnan(dist3));
                if ~isempty(okayidx)
                    if size(okayidx,2) > size(okayidx,1)
                        okayidx = transpose(okayidx);
                    end

                    matchcount = size(okayidx,1);
                    tbl_start = table_pos;
                    tbl_end = table_pos + matchcount - 1;
                    combo_table(tbl_start:tbl_end, 1) = dist3(okayidx);
                    combo_table(tbl_start:tbl_end, 4) = dist2(okayidx);
                    combo_table(tbl_start:tbl_end, 5) = zdist(okayidx);

                    [cand_rows, cand_cols] = ind2sub(size(dist3), okayidx);
                    combo_table(tbl_start:tbl_end, 2) = cand_rows + rr - 1;
                    combo_table(tbl_start:tbl_end, 3) = cand_cols + cc - 1;
                    table_pos = tbl_end + 1;
                end
                clear zdist dist2

                cc = c_end + 1;
            end
            rr = r_end + 1;
        end
        clear tempmtx_rx tempmtx_ry tempmtx_rz dist3

        cand_table = combo_table(1:table_pos-1,:);
        cand_table = sortrows(cand_table);
    end

    %%
    function cand_table = findMatchCandidates2D(callmtx, ref_set, snaprad)
        BLOCK_SIZE = 1024;
        rcount = size(callmtx, 1);
        ccount = size(ref_set, 1);
        R = ceil(rcount/BLOCK_SIZE);
        C = ceil(ccount/BLOCK_SIZE);

        alloc = ceil(rcount * ccount * 0.001);
        combo_table = single(NaN(alloc,3)); %dist row col
        table_pos = 1;

        %Go through blocks to get combo candidates
        rr = 1;
        for rb = 1:R
            cc = 1;
            r_end = rr + BLOCK_SIZE - 1;
            if r_end > rcount
                r_end = rcount;
            end
        
            tempmtx_rx = single(repmat(callmtx(rr:r_end,1), [1 BLOCK_SIZE]));
            tempmtx_ry = single(repmat(callmtx(rr:r_end,2), [1 BLOCK_SIZE]));
            rb_size = r_end - rr + 1;
            for cb = 1:C
                c_end = cc + BLOCK_SIZE - 1;
                %cb_size = BLOCK_SIZE;
                if c_end > ccount
                    c_end = ccount;
                    cb_size = c_end - cc + 1;
                    tempmtx_rx = single(repmat(callmtx(rr:r_end,1), [1 cb_size]));
                    tempmtx_ry = single(repmat(callmtx(rr:r_end,2), [1 cb_size]));
                end
                tempmtx_cols = single(repmat(ref_set(cc:c_end,1).', [rb_size 1]));
                xdist = tempmtx_rx - tempmtx_cols;

                tempmtx_cols = single(repmat(ref_set(cc:c_end,2).', [rb_size 1]));
                ydist = tempmtx_ry - tempmtx_cols;

                dist2 = sqrt((xdist.^2) + (ydist.^2));
                clear xdist;
                clear ydist;

                %Eliminate any combos above radius thresholds
                dist2(dist2 > snaprad) = NaN;

                %Save remaining combos to the table.
                okayidx = find(~isnan(dist2));
                if ~isempty(okayidx)
                    if size(okayidx,2) > size(okayidx,1)
                        okayidx = transpose(okayidx);
                    end

                    matchcount = size(okayidx,1);
                    tbl_start = table_pos;
                    tbl_end = table_pos + matchcount - 1;
                    combo_table(tbl_start:tbl_end, 1) = dist2(okayidx);
                    [cand_rows, cand_cols] = ind2sub(size(dist2), okayidx);
                    combo_table(tbl_start:tbl_end, 2) = cand_rows + rr - 1;
                    combo_table(tbl_start:tbl_end, 3) = cand_cols + cc - 1;
                    table_pos = tbl_end + 1;
                end

                cc = c_end + 1;
            end
            rr = r_end + 1;
        end
        clear tempmtx_rx tempmtx_ry dist2

        cand_table = combo_table(1:table_pos-1,:);
        cand_table = sortrows(cand_table);
    end

    %%
    function ref_assign = matchALotOfSpots(callmtx, ref_set, snaprad_3, snaprad_z)

        ccount = size(ref_set, 1);
        combo_table = RNACoords.findMatchCandidates(callmtx, ref_set, snaprad_3, snaprad_z);

        %Go through combos
        ref_assign = uint32(zeros(1,ccount));
        table_end = size(combo_table,1);
        table_pos = 1;
        while(table_pos <= table_end)
            if isnan(combo_table(table_pos,1))
                table_pos = table_pos + 1;
                continue;
            end

            my_row = combo_table(table_pos,2);
            my_col = combo_table(table_pos,3);
            ref_assign(my_col) = my_row;
            table_pos = table_pos + 1;

            %Nan out other candidates from this row and col
            outidxs = find(combo_table(table_pos:table_end, 2) == my_row) + table_pos - 1; 
            combo_table(outidxs, :) = NaN;

            outidxs = find(combo_table(table_pos:table_end, 3) == my_col) + table_pos - 1; 
            combo_table(outidxs, :) = NaN;

            %Cull nans from end
            table_end = find(~isnan(combo_table(table_pos:table_end, 1)), 1, 'last') + table_pos - 1;
            if isempty(table_end)
                break;
            end
        end
    end

    %%
    function ref_assign = matchALotOfSpots2D(callmtx, ref_set, snaprad)

        ccount = size(ref_set, 1);
        combo_table = RNACoords.findMatchCandidates2D(callmtx, ref_set, snaprad);

        %Go through combos
        ref_assign = uint32(zeros(1,ccount));
        table_end = size(combo_table,1);
        table_pos = 1;
        while(table_pos <= table_end)
            if isnan(combo_table(table_pos,1))
                table_pos = table_pos + 1;
                continue;
            end

            my_row = combo_table(table_pos,2);
            my_col = combo_table(table_pos,3);
            ref_assign(my_col) = my_row;
            table_pos = table_pos + 1;

            %Nan out other candidates from this row and col
            outidxs = find(combo_table(table_pos:table_end, 2) == my_row) + table_pos - 1; 
            combo_table(outidxs, :) = NaN;

            outidxs = find(combo_table(table_pos:table_end, 3) == my_col) + table_pos - 1; 
            combo_table(outidxs, :) = NaN;

            %Cull nans from end
            table_end = find(~isnan(combo_table(table_pos:table_end, 1)), 1, 'last') + table_pos - 1;
            if isempty(table_end)
                break;
            end
        end
    end

    %%
    function ref_assign =  matchSpots(call_table, ref_set, snaprad_3, snaprad_z, snap_minth)

        %Convert callxyz to a matrix
        ref_spots = size(ref_set,1);
        call_spots = size(call_table,1);
        callmtx = single(zeros(call_spots,3));
        callmtx(:,1) = single(table2array(call_table(:,'isnap_x')));
        callmtx(:,2) = single(table2array(call_table(:,'isnap_y')));
        callmtx(:,3) = single(table2array(call_table(:,'isnap_z')));

        %NaN out any rows where dropout th is below minimum
        dropth = table2array(call_table(:,'dropout_thresh'));
        callmtx(dropth < snap_minth, :) = NaN;

        %If the tables are huge, shunt to matchALotOfSpots
        if (call_spots > 1024) | (ref_spots > 1024)
            ref_assign = RNACoords.matchALotOfSpots(callmtx, ref_set, snaprad_3, snaprad_z);
            clear callmtx;
            return;
        end

        tempmtx_r = single(repmat(ref_set(:,1).', [call_spots 1]));
        tempmtx_c = single(repmat(callmtx(:,1), [1 ref_spots]));
        xdist = tempmtx_r - tempmtx_c;

        tempmtx_r = single(repmat(ref_set(:,2).', [call_spots 1]));
        tempmtx_c = single(repmat(callmtx(:,2), [1 ref_spots]));
        ydist = tempmtx_r - tempmtx_c;

        tempmtx_r = single(repmat(ref_set(:,3).', [call_spots 1]));
        tempmtx_c = single(repmat(callmtx(:,3), [1 ref_spots]));
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
        ref_assign = uint32(zeros(1,ref_spots)); %Index of call spot assigned to ref spot
        while nnz(~isnan(dist3)) > 0
            [~, minidx] = min(dist3, [], 'all', 'omitnan');
            [r, c] = ind2sub([call_spots ref_spots], minidx);
            dist3(:,c) = NaN;
            dist3(r,:) = NaN;
            ref_assign(c) = r;
        end
        clear dist3;

    end

    %%
    function ref_assign =  matchSpots2D(call_table, ref_set, snaprad, snap_minth)

        %Convert callxyz to a matrix
        ref_spots = size(ref_set,1);
        call_spots = size(call_table,1);
        callmtx = single(zeros(call_spots,3));
        callmtx(:,1) = single(call_table{:,'isnap_x'});
        callmtx(:,2) = single(call_table{:,'isnap_y'});

        %NaN out any rows where dropout th is below minimum
        dropth = call_table{:,'dropout_thresh'};
        callmtx(dropth < snap_minth, :) = NaN;

        %If the tables are huge, shunt to matchALotOfSpots
        if (call_spots > 1024) | (ref_spots > 1024)
            ref_assign = RNACoords.matchALotOfSpots2D(callmtx, ref_set, snaprad);
            clear callmtx;
            return;
        end

        tempmtx_r = single(repmat(ref_set(:,1).', [call_spots 1]));
        tempmtx_c = single(repmat(callmtx(:,1), [1 ref_spots]));
        xdist = tempmtx_r - tempmtx_c;

        tempmtx_r = single(repmat(ref_set(:,2).', [call_spots 1]));
        tempmtx_c = single(repmat(callmtx(:,2), [1 ref_spots]));
        ydist = tempmtx_r - tempmtx_c;
        clear tempmtx_r tempmtx_c

        dist2 = sqrt((xdist.^2) + (ydist.^2));
        clear xdist;
        clear ydist;

        %Eliminate any combos above radius thresholds
        dist2(dist2 > snaprad) = NaN;

        %Try to make matches
        ref_assign = uint32(zeros(1,ref_spots)); %Index of call spot assigned to ref spot
        while nnz(~isnan(dist2)) > 0
            [~, minidx] = min(dist2, [], 'all', 'omitnan');
            [r, c] = ind2sub([call_spots ref_spots], minidx);
            dist2(:,c) = NaN;
            dist2(r,:) = NaN;
            ref_assign(c) = r;
        end
        clear dist2;

    end

    %%
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
            temparr = double(call_table{cidxs, 'isnap_x'});
            xdist = abs(double(ref_set(fidxs,1)) - temparr);
            temparr = double(call_table{cidxs, 'isnap_y'});
            ydist = abs(double(ref_set(fidxs,2)) - temparr);
            temparr = double(call_table{cidxs, 'isnap_z'});
            zdist = abs(double(ref_set(fidxs,3)) - temparr);

            xydist = sqrt((xdist .^ 2) + (ydist .^ 2));
            xyzdist = sqrt((xdist .^ 2) + (ydist .^ 2) + (zdist .^ 2));

            call_table{cidxs, 'xdist_ref'} = single(xdist);
            call_table{cidxs, 'ydist_ref'} = single(ydist);
            call_table{cidxs, 'zdist_ref'} = single(zdist);
            call_table{cidxs, 'xydist_ref'} = single(xydist);
            call_table{cidxs, 'xyzdist_ref'} = single(xyzdist);

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
            tblappend{:,'isnap_x'} = uint16(x_round);
            tblappend{:,'isnap_y'} = uint16(y_round);
            tblappend{:,'isnap_z'} = uint16(z_round);
            tblappend{:,'dropout_thresh'} = 0;
            tblappend{:,'is_true'} = true;

            %Add distances
            xdist = abs(x_round - double(ref_set(ridxs,1)));
            ydist = abs(y_round - double(ref_set(ridxs,2)));
            zdist = abs(z_round - double(ref_set(ridxs,3)));
            xydist = sqrt((xdist .^ 2) + (ydist .^ 2));
            xyzdist = sqrt((xdist .^ 2) + (ydist .^ 2) + (zdist .^ 2));

            tblappend{:, 'xdist_ref'} = single(xdist);
            tblappend{:, 'ydist_ref'} = single(ydist);
            tblappend{:, 'zdist_ref'} = single(zdist);
            tblappend{:, 'xydist_ref'} = single(xydist);
            tblappend{:, 'xyzdist_ref'} = single(xyzdist);

            call_table = [call_table;tblappend];
        end

    end

    %%
    function [call_table, ref_assign] = updateTFCalls2D(call_table, ref_set, snaprad, snap_minth)
        if nargin < 3
            snaprad = 4;
        end
        if nargin < 4
            snap_minth = 25;
        end

        call_table{:,'is_true'} = false;
        ref_assign = RNACoords.matchSpots2D(call_table, ref_set, snaprad, snap_minth);

        %Mark trues in table
        if nnz(ref_assign) > 0
            fidxs = find(ref_assign > 0);
            cidxs = ref_assign(fidxs);
            call_table{cidxs,'is_true'} = true;

            %Record distances from match
            temparr = double(call_table{cidxs, 'isnap_x'});
            xdist = abs(double(ref_set(fidxs,1)) - temparr);
            temparr = double(call_table{cidxs, 'isnap_y'});
            ydist = abs(double(ref_set(fidxs,2)) - temparr);

            xydist = sqrt((xdist .^ 2) + (ydist .^ 2));

            call_table{cidxs, 'xdist_ref'} = single(xdist);
            call_table{cidxs, 'ydist_ref'} = single(ydist);
            call_table{cidxs, 'xydist_ref'} = single(xydist);

            clear xdist ydist xydist temparr
        end

        %Import fnegs
        no_assign = (ref_assign == 0);
        fn_count = nnz(no_assign);
        if fn_count > 0
            ridxs = find(no_assign);
            tblappend = RNACoords.genAppendableTable(call_table, fn_count);
            
            x_round = double(round(ref_set(ridxs,1)));
            y_round = double(round(ref_set(ridxs,2)));
            tblappend{:,'isnap_x'} = uint16(x_round);
            tblappend{:,'isnap_y'} = uint16(y_round);
            tblappend{:,'isnap_z'} = 1;
            tblappend{:,'dropout_thresh'} = 0;
            tblappend{:,'is_true'} = true;

            %Add distances
            xdist = abs(x_round - double(ref_set(ridxs,1)));
            ydist = abs(y_round - double(ref_set(ridxs,2)));
            xydist = sqrt((xdist .^ 2) + (ydist .^ 2));

            tblappend{:, 'xdist_ref'} = single(xdist);
            tblappend{:, 'ydist_ref'} = single(ydist);
            tblappend{:, 'xydist_ref'} = single(xydist);

            call_table = [call_table;tblappend];
        end

    end

    %%
    function [call_table_merged, callmap] = mergeSlicedSetTo3D(call_table, snaprad_3, snap_minth)
        call_table_merged = call_table;

        call_spots = size(call_table,1);
        callmtx = zeros(call_spots,3);
        callmtx(:,1) = single(table2array(call_table(:,'isnap_x')));
        callmtx(:,2) = single(table2array(call_table(:,'isnap_y')));
        callmtx(:,3) = single(table2array(call_table(:,'isnap_z')));
        callmap = uint32(zeros(1,call_spots)); %Which index each spot is remapped to

        above_spots = zeros(1,call_spots);
        below_spots = zeros(1,call_spots);

        %NaN out any rows where dropout th is below minimum
        dropth = table2array(call_table(:,'dropout_thresh'));
        callmtx(dropth < snap_minth, :) = NaN;

        cand_table = RNACoords.findMatchCandidates(callmtx, callmtx, snaprad_3, 1);
        if isempty(cand_table)
            return; 
        end

        %Recalculate zdist for these candidates
        z_rows = callmtx(cand_table(:,2), 3);
        z_cols = callmtx(cand_table(:,3), 3);
        z_dist = z_rows - z_cols;
        clear z_rows z_cols;
        z0 = find(z_dist(:,1) == 0);
        z_dist(z0) = [];
        cand_table(z0, :) = [];

        %Go through remaining combos and mark best candidates above and
        %below
        table_pos = 1;
        table_end = size(cand_table,1);
        while table_pos <= table_end
            if isnan(cand_table(table_pos, 1))
                table_pos = table_pos + 1;
                continue;
            end

            zdiff = z_dist(table_pos, 1);
            if zdiff == 0
                table_pos = table_pos + 1;
                continue;
            end

            r = cand_table(table_pos, 2);
            c = cand_table(table_pos, 3);
            if zdiff > 0
                %Row above col
                if (below_spots(r) == 0) & (above_spots(c) == 0)
                    below_spots(r) = c;
                    above_spots(c) = r;
                end
            else
                %Col above row
                if (below_spots(c) == 0) & (above_spots(r) == 0)
                    %Neither is already paired elsewhere
                    below_spots(c) = r;
                    above_spots(r) = c;
                end
            end

            table_pos = table_pos + 1;

            %NaN out
%             
%             outidxs = find(combo_table(table_pos:table_end, 2) == my_row) + table_pos - 1; 
%             cand_table(outidxs, :) = NaN;
% 
%             outidxs = find(combo_table(table_pos:table_end, 3) == my_col) + table_pos - 1; 
%             cand_table(outidxs, :) = NaN;
% 
%             %Cull nans from end
%             table_end = find(~isnan(cand_table(table_pos:table_end, 1)), 1, 'last') + table_pos - 1;
%             if isempty(table_end)
%                 break;
%             end
        end
        clear cand_table zdiff z_dist table_pos table_end

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

    %%
    function xyz_table = getThresholdCalls(call_table, thval, include_thvals)
        if nargin < 3; include_thvals = false; end

        xyz_table = [];
        if isempty(call_table); return; end

        colcount = 3;
        if include_thvals; colcount = 4; end
        tokay = find(call_table{:, 'dropout_thresh'} >= thval);
        if ~isempty(tokay)
            count = size(tokay, 1);
            xyz_table = int32(zeros(count, colcount));
            xyz_table(:,1) = call_table{tokay, 'isnap_x'};
            xyz_table(:,2) = call_table{tokay, 'isnap_y'};
            xyz_table(:,3) = call_table{tokay, 'isnap_z'};
            
            if include_thvals
                xyz_table(:,4) = call_table{tokay, 'dropout_thresh'};
            end
        end
    end

    %%
    function call_table = addFitData(call_table, fit_table, snap_minth)
        %Fit table should be formatted with cols:
        %   x,y,z,int(optional)

        %Allow minimum th (for bigfish)
        if nargin < 3; snap_minth = 0; end

        snaprad_3 = 4;
        snaprad_z = 2;
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

    %%
    function call_table = addFitDataFromQuant(call_table, cell_rna_data, snap_minth, img_size)
        %Count fits in cell_rna_data
        if isempty(call_table); return; end
        if isempty(cell_rna_data); return; end

        %Add zfitq column, if not present
        if ~ismember('zfitq', call_table.Properties.VariableNames)
            call_table.zfitq = NaN(size(call_table,1), 1);
        end

        snaprad_3 = 4;
        snaprad_z = 2;

        cell_count = size(cell_rna_data,2);
        spot_count = 0;

        for c = 1:cell_count
            my_cell = cell_rna_data(c);
            if ~isempty(my_cell.spots)
                spot_count = spot_count + size(my_cell.spots,2);
            end
        end
        if spot_count < 1; return; end

        %x y z intmax inttot xfwhm yfwhm zfitq 1dcoord
        fit_table = NaN(spot_count, 9);
        table_pos = 1;
        for c = 1:cell_count
            my_cell = cell_rna_data(c);
            if ~isempty(my_cell.spots)
                cell_spot_count = size(my_cell.spots,2);
                for s = 1:cell_spot_count
                    myspot = my_cell.spots(s);

                    spot_fit = myspot.gauss_fit;
                    if isempty(spot_fit); continue; end
                    
                    fit_table(table_pos,1) = spot_fit.xfit + my_cell.cell_loc.left - 1;
                    fit_table(table_pos,2) = spot_fit.yfit + my_cell.cell_loc.top - 1;
                    fit_table(table_pos,3) = spot_fit.zabs + my_cell.cell_loc.z_bottom - 1;
                    fit_table(table_pos,4) = spot_fit.fitMInt;
                    fit_table(table_pos,5) = spot_fit.TotFitInt;
                    fit_table(table_pos,6) = spot_fit.xFWHM;
                    fit_table(table_pos,7) = spot_fit.yFWHM;
                    fit_table(table_pos,8) = spot_fit.zqFit; %Quality of zfit.  If <0 pass if >0 curvature is positive at some point and fit is fail

                    xinit = spot_fit.xinit + my_cell.cell_loc.left - 1;
                    yinit = spot_fit.yinit + my_cell.cell_loc.top - 1;
                    zinit = spot_fit.zinit + my_cell.cell_loc.z_bottom - 1;
                    fit_table(table_pos,9) = sub2ind(img_size, yinit, xinit, zinit);

                    table_pos = table_pos + 1;
                end
            end
        end

        %Match 1 - Look for exact matches
        [exact_match, ref_assign] = ismember(fit_table(:,9), call_table{:,'coord_1d'});
        if nnz(exact_match) > 0
            rcount = size(ref_assign,1);
            for i = 1:rcount
                if ref_assign(i) < 1; continue; end
                cidx = ref_assign(i);
                call_table{cidx, 'fit_x'} = fit_table(i,1);
                call_table{cidx, 'fit_y'} = fit_table(i,2);
                call_table{cidx, 'fit_z'} = fit_table(i,3);
                call_table{cidx, 'fit_intensity'} = fit_table(i,4);
                call_table{cidx, 'fit_total_intensity'} = fit_table(i,5);
                call_table{cidx, 'xFWHM'} = fit_table(i,6);
                call_table{cidx, 'yFWHM'} = fit_table(i,7);
                call_table{cidx, 'zfitq'} = fit_table(i,8);
            end
        end

        %Match remaining
        if nnz(~exact_match) > 0
            ref_assign =  RNACoords.matchSpots(call_table, fit_table, snaprad_3, snaprad_z, snap_minth);
            rcount = size(ref_assign,2);

            %Copy to call table
            for i = 1:rcount
                if ref_assign(i) < 1; continue; end
                cidx = ref_assign(i);
                call_table{cidx, 'fit_x'} = fit_table(i,1);
                call_table{cidx, 'fit_y'} = fit_table(i,2);
                call_table{cidx, 'fit_z'} = fit_table(i,3);
                call_table{cidx, 'fit_intensity'} = fit_table(i,4);
                call_table{cidx, 'fit_total_intensity'} = fit_table(i,5);
                call_table{cidx, 'xFWHM'} = fit_table(i,6);
                call_table{cidx, 'yFWHM'} = fit_table(i,7);
                call_table{cidx, 'zfitq'} = fit_table(i,8);
            end
        end
        
    end

    %%
    function call_table = updateRefDistances(call_table, ref_table, ref_call_map)
        x_true = double(call_table{ref_call_map, 'isnap_x'});
        y_true = double(call_table{ref_call_map, 'isnap_y'});
        z_true = double(call_table{ref_call_map, 'isnap_z'});

        x_dist = abs(x_true - ref_table(:,1));
        y_dist = abs(y_true - ref_table(:,2));
        z_dist = abs(z_true - ref_table(:,3));

        call_table{ref_call_map, 'xdist_ref'} = x_dist;
        call_table{ref_call_map, 'ydist_ref'} = y_dist;
        call_table{ref_call_map, 'zdist_ref'} = z_dist;

        xsq = x_dist .^ 2;
        ysq = y_dist .^ 2;
        zsq = z_dist .^ 2;
        call_table{ref_call_map, 'xydist_ref'} = sqrt(xsq + ysq);
        call_table{ref_call_map, 'xyzdist_ref'} = sqrt(xsq + ysq + zsq);
    end

    %%
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

    %%
    function call_table = applyCellSegMask(call_table, cell_mask)
        x = table2array(call_table(:,'isnap_x'));
        y = table2array(call_table(:,'isnap_y'));

        X = size(cell_mask, 2);
        Y = size(cell_mask, 1);

        if ndims(cell_mask) >= 3
            z = table2array(call_table(:,'isnap_z'));
            Z = size(cell_mask, 3);

            ind = sub2ind([Y X Z], y, x, z);
            cells = cell_mask(ind);
        else
            ind = sub2ind([Y X], y, x);
            cells = cell_mask(ind);
        end

        call_table(:,'cell') = array2table(uint16(cells));
    end

    %%
    function gen_coords = genSpotsRandom(X, Y, Z, count)
        gen_coords = zeros(count, 3);
        
        rr = rand(count, 1);
        rr = (rr * X);
        rr = rr + 1;
        gen_coords(:,1) = floor(rr);

        rr = rand(count, 1);
        rr = (rr * Y);
        rr = rr + 1;
        gen_coords(:,2) = floor(rr);

        rr = rand(count, 1);
        rr = (rr * Z);
        rr = rr + 1;
        gen_coords(:,3) = floor(rr);
    end

    %%
    function gen_coords = genSpotsRandomWeightedZ(X, Y, Z_weights, count)
        gen_coords = zeros(count, 3);
        
        rr = rand(count, 1);
        rr = (rr * X);
        rr = rr + 1;
        gen_coords(:,1) = floor(rr);

        rr = rand(count, 1);
        rr = (rr * Y);
        rr = rr + 1;
        gen_coords(:,2) = floor(rr);

        zsum = sum(Z_weights, 'all');
        Z_weights = Z_weights ./ zsum;
        zcsum = cumsum(Z_weights);
        rr = rand(count, 1);
        cmtx = rr > zcsum;
        gen_coords(:,3) = sum(cmtx, 2);

    end

    %%
    function nearby_counts = countNearbySpotsOneChannel(call_mtx, radXY, radZ)
        nearby_counts = RNACoords.countNearbySpotsOtherChannel(call_mtx, call_mtx, radXY, radZ);
    end

    %%
    function nearby_counts = countNearbySpotsTwoChannel(call_mtx, other_mtx, rad3, radZ)
        spot_count = size(call_mtx, 1);

        cand_table = RNACoords.findMatchCandidates(call_mtx, other_mtx, rad3, radZ);
        bins = [1:1:spot_count+1];
        [nearby_counts, ~] = histcounts(cand_table(:,2), bins);
    end

    %%
    function [rad3, radZ, colocRes] = colocAutoRadiusTwoChannel(callmtx_ch_A, callmtx_ch_B, idims, voxdims_nm, maxR_nm, maxZR_nm)
        if nargin < 5; maxR_nm = 1000; end
        if nargin < 6; maxZR_nm = 2000; end

        %If maxZR_nm is set to 0, then collapse to 2D.
        if maxZR_nm <= 0
            maxZR_nm = 0;
            if size(callmtx_ch_A, 2) < 3
                aold = callmtx_ch_A;
                callmtx_ch_A = zeros(size(callmtx_ch_A,1), 3);
                callmtx_ch_A(:,1:2) = aold(:,1:2);
                clear aold;
            end
            callmtx_ch_A(:,3) = 1;

            if size(callmtx_ch_B, 2) < 3
                bold = callmtx_ch_B;
                callmtx_ch_B = zeros(size(callmtx_ch_B,1), 3);
                callmtx_ch_B(:,1:2) = bold(:,1:2);
                clear bold;
            end
            callmtx_ch_B(:,3) = 1;
        end

        if isstruct(idims)
            X = idims.x;
            Y = idims.y;
            Z = idims.z;
        else
            X = 1; Y = 1; Z = 1;
            sz = size(idims, 2);
            if sz >= 1; X = idims(1); end
            if sz >= 2; Y = idims(2); end
            if sz >= 3; Z = idims(3); end
        end

        if isstruct(voxdims_nm)
            vx = voxdims_nm.x;
            vy = voxdims_nm.y;
            vz = voxdims_nm.z;
        else
            vx = 65; vy = 65; vz = 300;
            sz = size(voxdims_nm, 2);
            if sz >= 1; vx = voxdims_nm(1); end
            if sz >= 2; vy = voxdims_nm(2); end
            if sz >= 3; vz = voxdims_nm(3); end
        end

        if maxZR_nm <= 0
            Z = 1;
        end

        rad3 = 0;
        radZ = 0;

        inCount = size(callmtx_ch_A,1);
        colocRes = cell(1,inCount);
        randCount = size(callmtx_ch_B,1);

        %Weight Z for random?
        if maxZR_nm > 0
            [Zw, ~] = histcounts(callmtx_ch_B(:,3), [1:Z+1]);
            Zw = Zw ./ randCount;
            randCoords = RNACoords.genSpotsRandomWeightedZ(X, Y, Zw, randCount);
            %randCoords = RNACoords.genSpotsRandom(X, Y, Z, randCount);
        else
            randCoords = RNACoords.genSpotsRandom(X, Y, Z, randCount);
        end

        %Convert all coords to nm (so can better weigh Z)
        callmtx_ch_A_nm = callmtx_ch_A - 1;
        callmtx_ch_A_nm(:,1) = callmtx_ch_A_nm(:,1) .* vx;
        callmtx_ch_A_nm(:,2) = callmtx_ch_A_nm(:,2) .* vy;
        callmtx_ch_A_nm(:,3) = callmtx_ch_A_nm(:,3) .* vz;

        callmtx_ch_B_nm = callmtx_ch_B - 1;
        callmtx_ch_B_nm(:,1) = callmtx_ch_B_nm(:,1) .* vx;
        callmtx_ch_B_nm(:,2) = callmtx_ch_B_nm(:,2) .* vy;
        callmtx_ch_B_nm(:,3) = callmtx_ch_B_nm(:,3) .* vz;

        randCoords_nm = randCoords - 1;
        randCoords_nm(:,1) = randCoords_nm(:,1) .* vx;
        randCoords_nm(:,2) = randCoords_nm(:,2) .* vy;
        randCoords_nm(:,3) = randCoords_nm(:,3) .* vz;

        randMatches = RNACoords.findMatchCandidates(callmtx_ch_A_nm, randCoords_nm, maxR_nm, maxZR_nm);
        callMatches = RNACoords.findMatchCandidates(callmtx_ch_A_nm, callmtx_ch_B_nm, maxR_nm, maxZR_nm);
        alldist = [randMatches(:,1)' callMatches(:,1)'];
        alldist = unique(sort(alldist));
        %alldist2 = [randMatches(:,4)' callMatches(:,4)'];
        %alldist2 = unique(sort(alldist2));
        %alldistz = [0:1:maxZR];
        alldistz = [randMatches(:,5)' callMatches(:,5)'];
        alldistz = unique(sort(alldistz));

        hbinsA = [1:1:(inCount+1)];
        %hbinsB = [1:1:(randCount+1)];

        d3Count = size(alldist, 2);
        %d2Count = size(alldist2, 2);
        dzCount = size(alldistz, 2);
        callAvgCounts = zeros(d3Count, dzCount);
        randAvgCounts = zeros(d3Count, dzCount);
        for i = 1:d3Count
            d3Val = alldist(i);
            for j = 1:dzCount
                dzVal = alldistz(j);

                fbool = and(callMatches(:,1) <= d3Val, callMatches(:,5) <= dzVal);
                cm_filt = callMatches(fbool, :);
                fbool = and(randMatches(:,1) <= d3Val, randMatches(:,5) <= dzVal);
                rn_filt = randMatches(fbool, :);

                if ~isempty(cm_filt)
                    [hh, ~] = histcounts(cm_filt(:,2), hbinsA);
                    callAvgCounts(i,j) = mean(hh, 'all', 'omitnan');
                    %callMedCounts(i) = median(hh, 'all', 'omitnan');
                    clear hh
                end

                if ~isempty(rn_filt)
                    [hh, ~] = histcounts(rn_filt(:,2), hbinsA);
                    randAvgCounts(i,j) = mean(hh, 'all', 'omitnan');
                    %randMedCounts(i) = median(hh, 'all', 'omitnan');
                    clear hh
                end
            end
        end

        lrAvgCounts = log10(randAvgCounts);
        lcAvgCounts = log10(callAvgCounts);



        %TODO

    end

end

end