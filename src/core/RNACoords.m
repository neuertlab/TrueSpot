%
%%

classdef RNACoords

methods (Static)

    function call_table = convertOldCoordTable(coord_table_old, image_filtered, image_raw)

    end

    function call_table = genNewCoordTable(alloc)
        varNames = {'coord_1d' 'isnap_x' 'isnap_y' 'isnap_z' 'intensity_f' 'intensity'...
            'dropout_thresh', 'is_true'};
        varTypes = {'uint32' 'uint16' 'uint16' 'uint16' 'single' 'single'...
            'single' 'logical'};

        table_size = [alloc size(varNames,2)];
        call_table = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);

        dummy16 = array2table(uint16(zeros(alloc,1)));
        call_table(:,'isnap_x') = dummy16;
        call_table(:,'isnap_y') = dummy16;
        call_table(:,'isnap_z') = dummy16;
        clear dummy16;

        dummy32 = array2table(uint32(zeros(alloc,1)));
        call_table(:,'coord_1d') = dummy32;
        clear dummy32;

        dummyf = single(NaN(alloc,1));
        call_table(:,'intensity_f') = dummyf;
        call_table(:,'intensity') = dummyf;
        call_table(:,'dropout_thresh') = dummyf;
        clear dummyf;

        dummyb = false(alloc,1);
        call_table(:,'is_true') = dummyb;
        clear dummyb;
    end


end

end