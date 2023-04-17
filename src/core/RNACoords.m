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

    function call_table = genNewCoordTable(alloc, fit_incl)
        if nargin < 2
            fit_incl = 0;
        end

        varNames = {'coord_1d' 'isnap_x' 'isnap_y' 'isnap_z' 'intensity_f' 'intensity'...
            'dropout_thresh', 'is_true', 'is_trimmed_out', 'in_truth_region', 'cell'};
        varTypes = {'uint32' 'uint16' 'uint16' 'uint16' 'single' 'single'...
            'single' 'logical', 'logical', 'logical', 'uint16'};

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

        dummy16 = array2table(uint16(zeros(alloc,1)));
        call_table(:,'isnap_x') = dummy16;
        call_table(:,'isnap_y') = dummy16;
        call_table(:,'isnap_z') = dummy16;
        call_table(:,'cell') = dummy16;
        clear dummy16;

        dummy32 = array2table(uint32(zeros(alloc,1)));
        call_table(:,'coord_1d') = dummy32;
        clear dummy32;

        dummyf = array2table(single(NaN(alloc,1)));
        call_table(:,'intensity_f') = dummyf;
        call_table(:,'intensity') = dummyf;
        call_table(:,'dropout_thresh') = dummyf;
        if fit_incl >= 1
            call_table(:,'fit_x') = dummyf;
            call_table(:,'fit_y') = dummyf;
            call_table(:,'fit_z') = dummyf;
            call_table(:,'fit_intensity') = dummyf;
        end
        clear dummyf;

        if fit_incl >= 2
            dummyd = array2table(double(NaN(alloc,1)));
            call_table(:,'xFWHM') = dummyd;
            call_table(:,'yFWHM') = dummyd;
            call_table(:,'fit_total_intensity') = dummyd;
            clear dummyd;
        end

        dummyb = array2table(false(alloc,1));
        call_table(:,'is_true') = dummyb;
        call_table(:,'is_trimmed_out') = dummyb;
        call_table(:,'in_truth_region') = dummyb;
        clear dummyb;
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

        dist3 = sqrt((xdist.^2) + (ydist.^2) + (zdist.^2));
        clear xdist;
        clear ydist;
        zdist = abs(zdist);

        %Eliminate any combos above radius thresholds
        dist3(zdist > snaprad_z) = NaN;
        dist3(dist3 > snaprad_3) = NaN;

        %Try to make matches
        ref_assign = zeros(1,ref_spots); %Index of call spot assigned to ref spot
        while nnz(~isnan(dist3)) > 0
            [~, minidx] = min(dist3, [], 'all', 'omitnan');
            [r, c] = ind2sub([call_spots ref_spots], minidx);
            dist3(:,c) = NaN;
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
            aidxs = ref_assign(ref_assign > 0);
            call_table{aidxs,'is_true'} = true;
        end

        %Import fnegs
        no_assign = (ref_assign == 0);
        fn_count = nnz(no_assign);
        if fn_count > 0
            ridxs = find(no_assign);
            fitlvl = 0;
            colcount = size(call_table,2);
            if colcount > 11; fitlvl = 1; end
            if colcount > 15; fitlvl = 2; end
            tblappend = RNACoords.genNewCoordTable(fn_count, fitlvl);
            
            tblappend(:,'isnap_x') = array2table(uint16(ref_set(ridxs, 1)));
            tblappend(:,'isnap_y') = array2table(uint16(ref_set(ridxs, 2)));
            tblappend(:,'isnap_z') = array2table(uint16(ref_set(ridxs, 3)));
            tblappend{:,'dropout_thresh'} = 0;
            tblappend{:,'is_true'} = true;

            call_table = [call_table;tblappend];
        end

    end

    function call_table = addFitData(call_table, cell_rna_data)
        %TODO
    end

end

end