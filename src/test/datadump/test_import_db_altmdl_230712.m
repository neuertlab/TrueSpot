%
%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%BaseDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
addpath('./test');
addpath('./test/datadump');

% ========================== Constants ==========================

START_INDEX = 67;
END_INDEX = 68;

OutputDir = [BaseDir filesep 'data' filesep 'results'];

SCRIPT_VER = 'v23.07.12.02';
COMPUTER_NAME = 'VU_NEUERTLAB_HOSPELB';

EXPTS_INITIALS = 'BH';

% ========================== Load csv Table ==========================

%InputTablePath = [BaseDir filesep 'test_images_simytc.csv'];
%InputTablePath = [BaseDir filesep 'test_images_simvarmass.csv'];
InputTablePath = [BaseDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

% ========================== Iterate through table entries ==========================
if ~isfolder(OutputDir)
    mkdir(OutputDir);
end

entry_count = size(image_table,1);

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

for r = START_INDEX:END_INDEX
    is_sim = false;
    myname = getTableValue(image_table, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

    hb_stem_base = getTableValue(image_table, r, 'OUTSTEM');
    hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
    srcpath_raw = getTableValue(image_table, r, 'IMAGEPATH');
    srcpath = [ImgDir replace(srcpath_raw, '/', filesep)];

    %Find output file
    GroupOutDir = [OutputDir filesep getSetOutputDirName(myname)];
    if ~isfolder(GroupOutDir)
        mkdir(GroupOutDir);
    end
    OutFilePath = [GroupOutDir filesep myname '_summary.mat'];
    if isfile(OutFilePath)
        load(OutFilePath, 'analysis');
    else
        analysis = struct('imgname', myname);
    end
    analysis.script_version_dbaldmdl = SCRIPT_VER;

    %Find truthset...
    ts_region = [];
    ref_coords = [];
    ref_coords = loadSimTruthset(image_table, r, ImgDir);
    if isempty(ref_coords)
        [ref_coords, ts_region] = loadExpTruthset(image_table, r, BaseDir);
    end

    if startsWith(myname, 'sim_') | startsWith(myname, 'simvar')
        is_sim = true;
    elseif startsWith(myname, 'rsfish_sim')
        is_sim = true;
    end

    %Load original image.
    if endsWith(srcpath, '.mat')
        [idir, ifname, ~] = fileparts(srcpath);
        srcpath = [idir filesep 'tif' filesep ifname '.tif'];
    end
    chcount = getTableValue(image_table, r, 'CH_TOTAL');
    trgch = getTableValue(image_table, r, 'CHANNEL');
    [channels, ~] = LoadTif(srcpath, chcount, [trgch], 1);
    my_image = channels{trgch,1};
    my_image = uint16(my_image);
    clear channels;
    X = size(my_image,2);
    Y = size(my_image,1);
    Z = size(my_image,3);

    %----------------- Import cellseg mask TODO
    cellseg_dir = getTableValue(image_table, r, 'CELLSEG_DIR');
    cellseg_sfx = getTableValue(image_table, r, 'CELLSEG_SFX');

    cellmask = [];
    if ~strcmp(cellseg_dir, '.')
        cellseg_dir_actual = [BaseDir replace(cellseg_dir, '/', filesep)];
        cellseg_path = [cellseg_dir_actual filesep 'Lab_' cellseg_sfx '.mat'];
        if isfile(cellseg_path)
            load(cellseg_path, 'cells');
            cellmask = cells;
            clear cells;
        end
    end

    rsdb_dir_ext = getRSDBGroupOutputDir(myname);
    if ~isempty(rsdb_dir_ext)
        fprintf('> Importing DeepBlink results...\n');
        db_stem = [BaseDir filesep 'data' filesep 'deepblink_simmdl' rsdb_dir_ext myname filesep 'DeepBlink_' myname];
        [db_dir, ~, ~] = fileparts(db_stem);
        coord_table_path = [db_stem '_coordTable.mat'];
        spot_table_path = [db_stem '_spotTable.mat'];

        if isfile(coord_table_path)
            load(coord_table_path, 'coord_table');
            load(spot_table_path, 'spot_table');

            call_table_raw = RNACoords.convertOldCoordTable(spot_table, coord_table, [], my_image, 1);
            if ~isempty(cellmask)
                call_table_raw = RNACoords.applyCellSegMask(call_table_raw, cellmask);
            end

            %Import Fit (BEFORE merge!)
            fit_table_path = findDBFitTable(db_dir, myname);
            if ~isempty(fit_table_path) & isfile(fit_table_path)
                import_table = readtable(fit_table_path,'Delimiter',',','ReadVariableNames',true,'Format',...
                    '%f%f%f%f');
                import_mtx = table2array(import_table);
                import_count = size(import_mtx,1);
                fit_table = NaN(import_count,3);
                fit_table(:,1:3) = import_mtx(:,[1 2 4]) + 1;
                call_table_raw = RNACoords.addFitData(call_table_raw, fit_table);

                clear import_table import_mtx fit_table import_count
            end

            %Merge slice calls for better truthset comparison
            [call_table, callmap] = RNACoords.mergeSlicedSetTo3D(call_table_raw, 4, 0.5);
            init_call_count = size(call_table,1);

            if ~isempty(ref_coords)
                [call_table, ref_call_map] = RNACoords.updateTFCalls(call_table, ref_coords, 4, 2, 0);
                full_call_count = size(call_table, 1);

                if full_call_count > init_call_count
                    %Fnegs added.
                    addst = init_call_count + 1;
                    added = full_call_count;
                    x = table2array(call_table(addst:added,'isnap_x'));
                    y = table2array(call_table(addst:added,'isnap_y'));
                    z = table2array(call_table(addst:added,'isnap_z'));
                    c1d = sub2ind([Y X Z], y, x, z);
                    call_table(addst:added,'coord_1d') = array2table(uint32(c1d));
                    call_table(addst:added,'intensity') = array2table(single(my_image(c1d)));

                    clear addst added x y z c1d
                end

                %Mask ts region (if applicable)
                if ~isempty(ts_region)
                    inside_ts_mask = true(full_call_count, 1);
                    x = table2array(call_table(:,'isnap_x'));
                    y = table2array(call_table(:,'isnap_y'));
                    z = table2array(call_table(:,'isnap_z'));

                    inside_ts_mask = and(inside_ts_mask, ~(x < ts_region.x0));
                    inside_ts_mask = and(inside_ts_mask, ~(x > ts_region.x1));
                    inside_ts_mask = and(inside_ts_mask, ~(y < ts_region.y0));
                    inside_ts_mask = and(inside_ts_mask, ~(y > ts_region.y1));
                    inside_ts_mask = and(inside_ts_mask, ~(z < ts_region.z0));
                    inside_ts_mask = and(inside_ts_mask, ~(z > ts_region.z1));
                    call_table(:,'in_truth_region') = array2table(inside_ts_mask);

                    clear x y z inside_ts_mask
                else
                    %All inside.
                    call_table{:,'in_truth_region'} = true;
                end

                if isfile(fit_table_path)
                    call_table = RNACoords.updateRefDistancesToUseFits(call_table, ref_coords, ref_call_map);
                end
            end

            %Save
            if ~isfield(analysis, 'results_db_simmdl')
                analysis.results_db_simmdl = struct('callset', table.empty());
            end
            analysis.results_db_simmdl.timestamp = datetime;
            analysis.results_db_simmdl.import_computer = COMPUTER_NAME;
            analysis.results_db_simmdl.callset = call_table;
            analysis.results_db_simmdl.callset_sliced = call_table_raw;
            analysis.results_db_simmdl.callmap_slice_merge = callmap;

            if ~isempty(ref_coords)
                %Calculate performance metrics
                analysis.results_db_simmdl = runstats(analysis.results_db_simmdl, spot_table, 0.001);
                analysis.results_db_simmdl.ref_call_map = ref_call_map;
            end

            clear call_table call_table_raw callmap ref_call_map

            %If experimental TS, mark ts fields
            if ~is_sim & ~isempty(ref_coords) & ~isempty(EXPTS_INITIALS)
                analysis.results_db_simmdl = markTSStats(analysis.results_db_simmdl, EXPTS_INITIALS);
            end
        else
            fprintf('ERROR: Could not find DeepBlink run for %s!\n', myname);
        end
        clear coord_table_path spot_table_path coord_table spot_table
    end

    %--------------------------------------- Tag exp ts data outside result
    %structs

    %--------------------------------------- Save
    fprintf('> Saving updated analysis...\n');
    save(OutFilePath, 'analysis', '-v7.3');
    clear analysis;
end

% ========================== Helper functions ==========================

function outdir = getRSDBGroupOutputDir(imgname)
    outdir = [];
    if isempty(imgname); return; end

    if startsWith(imgname, 'mESC4d_')
        outdir = [filesep 'mESC4d' filesep];
    elseif startsWith(imgname, 'scrna_')
        outdir = [filesep 'scrna' filesep];
    elseif startsWith(imgname, 'mESC_loday_')
        outdir = [filesep 'mESC_loday' filesep];
    elseif startsWith(imgname, 'scprotein_')
        outdir = [filesep 'scprotein' filesep];
    elseif startsWith(imgname, 'sim_')
        outdir = [filesep 'sim' filesep];
    elseif startsWith(imgname, 'histonesc_')
        outdir = [filesep 'histonesc' filesep];
    elseif startsWith(imgname, 'ROI0')
        outdir = [filesep 'munsky_lab' filesep];
    elseif startsWith(imgname, 'sctc_')
        inparts = split(imgname, '_');
        ch = 'CH1';
        if endsWith(imgname, 'CTT1')
            ch = 'CH2';
        end
        outdir = [filesep 'yeast_tc' filesep inparts{2,1} filesep ch filesep];
    elseif startsWith(imgname, 'rsfish_')
        outdir = [filesep 'rsfish' filesep];
    elseif startsWith(imgname, 'simvar_')
        outdir = [filesep 'simvar' filesep];
    elseif startsWith(imgname, 'simvarmass_')
        if contains(imgname, 'TMRL') | contains(imgname, 'CY5L')
            outdir = [filesep 'simytc' filesep];
        else
            outdir = [filesep 'simvarmass' filesep];
        end
    end

end

function dirname = getSetOutputDirName(imgname)
    inparts = split(imgname, '_');
    groupname = inparts{1,1};
    if strcmp(groupname, 'sctc')
        dirname = [groupname filesep inparts{2,1}];
    elseif strcmp(groupname, 'simvarmass')
        if contains(imgname, 'TMRL') | contains(imgname, 'CY5L')
            dirname = 'simytc';
        else
            dirname = groupname;
        end
    else
        if startsWith(groupname, 'ROI')
            dirname = 'munsky_lab';
        else
            dirname = groupname;
        end
    end
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

function fit_table_path = findDBFitTable(db_dir, image_name)
    fit_table_path = [db_dir filesep image_name];
    if ~isfile(fit_table_path)
        %Find another csv in that dir
        dir_contents = dir(db_dir);

        content_count = size(dir_contents,1);
        for i = 1:content_count
            if endsWith(dir_contents(i,1).name, '.csv')
                fit_table_path = [db_dir filesep dir_contents(i,1).name];
                break;
            end
        end
    end
end

function key_mtx = keyStructs2Mtx(my_key)
    ptcount = size(my_key,2);
    key_mtx = uint16(zeros(ptcount,3));
    key_mtx(:,1) = [my_key.x];
    key_mtx(:,2) = [my_key.y];
    key_mtx(:,3) = [my_key.z];
end

function ref_coords = loadSimTruthsetRS(image_table, row_index, ImgDir)
    myname = getTableValue(image_table, row_index, 'IMGNAME');
    srcpath_raw = getTableValue(image_table, row_index, 'IMAGEPATH');
    
    %I had to convert the locs to csv because MATLAB is a fussbudget
    srcpath = [ImgDir replace(srcpath_raw, '/', filesep)];
    srcpath = replace(srcpath, '.tif', '.csv');
    
    ref_coords = [];
    if ~isfile(srcpath)
        fprintf('ERROR: Could not find sim truthset for %s!\n', myname);
        return;
    end
    
    import_table = table2array(readtable(srcpath,'ReadVariableNames',false));
    import_table = import_table + 1;
    temp = import_table(:,2);
    import_table(:,2) = import_table(:,1);
    import_table(:,1) = temp;
    
    %Swap out refset and save
    ref_coords = import_table(:,1:3);
end

function ref_coords = loadSimTruthsetSF(image_table, row_index, ImgDir)
    %For simfish sims
    myname = getTableValue(image_table, row_index, 'IMGNAME');
    srcpath_raw = getTableValue(image_table, row_index, 'IMAGEPATH');
    
    srcpath = [ImgDir replace(srcpath_raw, '/', filesep)];
    
    key = [];
    if endsWith(srcpath, '.mat')
        %Use this one directly.
        load(srcpath, 'key');
    elseif endsWith(srcpath, '.tif')
        %Find mat file and load from that.
        srcpath = replace(srcpath, [filesep 'tif' filesep], filesep);
        srcpath = replace(srcpath, '.tif', '.mat');
        load(srcpath, 'key');
    end
    
    if isempty(key)
        fprintf('ERROR: Could not load sim truthset for %s!\n', myname);
        return; 
    end

    ref_coords = keyStructs2Mtx(key); %Importer already adjusts to 1 based coords.

end

function ref_coords = loadSimTruthset(image_table, row_index, ImgDir)
    ref_coords = [];
    myname = getTableValue(image_table, row_index, 'IMGNAME');
    if startsWith(myname, 'sim_')
        ref_coords = loadSimTruthsetSF(image_table, row_index, ImgDir);
    elseif startsWith(myname, 'simvar')
        ref_coords = loadSimTruthsetSF(image_table, row_index, ImgDir);
    elseif startsWith(myname, 'rsfish_sim')
        ref_coords = loadSimTruthsetRS(image_table, row_index, ImgDir);
    end
end

function [ref_coords, valid_range] = loadExpTruthset(image_table, row_index, BaseDir)
    ref_coords = [];
    valid_range = struct('x0', 0, 'x1', 0, 'y0', 0, 'y1', 0, 'z0', 0, 'z1', 0);
    
    hb_stem_base = getTableValue(image_table, row_index, 'OUTSTEM');
    hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];

    if RNA_Threshold_SpotSelector.refsetExists(hb_stem)
        spotanno = RNA_Threshold_SpotSelector.openSelector(hb_stem, true);
        ref_coords = spotanno.ref_coords;
        valid_range.z0 = spotanno.z_min;
        valid_range.z1 = spotanno.z_max;
        if ~isempty(spotanno.selmcoords)
            valid_range.x0 = spotanno.selmcoords(1,1);
            valid_range.x1 = spotanno.selmcoords(2,1);
            valid_range.y0 = spotanno.selmcoords(3,1);
            valid_range.y1 = spotanno.selmcoords(4,1);
        else
            valid_range.x0 = 1;
            valid_range.y0 = 1;
            valid_range.x1 = getTableValue(image_table, row_index, 'IDIM_X');
            valid_range.y1 = getTableValue(image_table, row_index, 'IDIM_Y');
        end
        clear spotanno;
    end

end

function rstruct = markTSStats(rstruct, tag)
    if isempty(tag); return; end

    substruct_name = ['truthset_' tag];
    if ~isfield(rstruct, substruct_name)
        rstruct.(substruct_name) = struct('timestamp', datetime);
    else
        rstruct.(substruct_name).timestamp = datetime;
    end

    if isfield(rstruct, 'callset')
        %Need to move over 'is_true' and 'in_truth_region'
        targcol_name = ['is_true_' tag];
        newcol_A = renamevars(rstruct.callset(:,'is_true'), 'is_true', targcol_name);
        targcol_name = ['in_truth_region_' tag];
        newcol_B = renamevars(rstruct.callset(:,'in_truth_region'), 'in_truth_region', targcol_name);
        rstruct.callset = [rstruct.callset newcol_A newcol_B];
        rstruct.callset{:,'is_true'} = false;
        rstruct.callset{:,'in_truth_region'} = false;
    end
    if isfield(rstruct, 'callset_sliced')
        %DB only
        targcol_name = ['is_true_' tag];
        newcol_A = renamevars(rstruct.callset_sliced(:,'is_true'), 'is_true', targcol_name);
        targcol_name = ['in_truth_region_' tag];
        newcol_B = renamevars(rstruct.callset_sliced(:,'in_truth_region'), 'in_truth_region', targcol_name);
        rstruct.callset_sliced = [rstruct.callset_sliced newcol_A newcol_B];
        rstruct.callset_sliced{:,'is_true'} = false;
        rstruct.callset_sliced{:,'in_truth_region'} = false;
    end
    if isfield(rstruct, 'ref_call_map')
        rstruct.(substruct_name).ref_call_map = rstruct.ref_call_map;
    end
    if isfield(rstruct, 'performance')
        rstruct.(substruct_name).performance = rstruct.performance;
        rstruct.(substruct_name).pr_auc = rstruct.pr_auc;
        rstruct.(substruct_name).fscore_peak = rstruct.fscore_peak;
        if isfield(rstruct, 'fscore_autoth')
            rstruct.(substruct_name).fscore_autoth = rstruct.fscore_autoth;
        end
    end
    if isfield(rstruct, 'performance_trimmed')
        rstruct.(substruct_name).performance_trimmed = rstruct.performance_trimmed;
        rstruct.(substruct_name).pr_auc_trimmed = rstruct.pr_auc_trimmed;
        rstruct.(substruct_name).fscore_peak_trimmed = rstruct.fscore_peak_trimmed;
        if isfield(rstruct, 'fscore_autoth_trimmed')
            rstruct.(substruct_name).fscore_autoth_trimmed = rstruct.fscore_autoth_trimmed;
        end
    end
    rstruct.last_ts = tag;

end

function rstruct = runstats(rstruct, spot_table, th_val)

    if nargin < 3; th_val = 0; end

    call_table = rstruct.callset;
    vec_istrimmed = table2array(call_table(:,'is_trimmed_out'));
    vec_intsreg = table2array(call_table(:,'in_truth_region'));
    vec_isreal = table2array(call_table(:,'is_true'));
    vec_dropth = table2array(call_table(:,'dropout_thresh'));

    any_trimmed = nnz(vec_istrimmed) > 0;
    
    th_count = size(spot_table,1);
    res_untrimmed = ImageResults.initializeResTable(th_count);
    res_trimmed = table.empty();

    thval_tbl = array2table(double(spot_table(:,1)));
    res_untrimmed(:,'thresholdValue') = thval_tbl;

    if any_trimmed
        res_trimmed = ImageResults.initializeResTable(th_count);
        res_trimmed(:,'thresholdValue') = thval_tbl;
    end

    sc_all = NaN(th_count,2);
    tp_all = NaN(th_count,2);
    fp_all = NaN(th_count,2);
    fn_all = NaN(th_count,2);
    for t = 1:th_count
        th = spot_table(t,1);
        pos_vec = (vec_dropth >= th) & vec_intsreg;
        tp_vec = pos_vec & vec_isreal;
        fp_vec = pos_vec & ~vec_isreal;
        fn_vec = (vec_dropth < th) & vec_intsreg & vec_isreal;

        tp_all(t,1) = nnz(tp_vec);
        fp_all(t,1) = nnz(fp_vec);
        fn_all(t,1) = nnz(fn_vec);
        sc_all(t,1) = nnz(pos_vec);

        %Repeat for trimmed, if applicable
        if any_trimmed
            pos_vec = (vec_dropth >= th) & vec_intsreg & ~vec_istrimmed;
            neg_vec = (vec_dropth < th) & vec_intsreg & ~vec_istrimmed;
            tp_vec = pos_vec & vec_isreal;
            fp_vec = pos_vec & ~vec_isreal;
            fn_vec = neg_vec & vec_isreal;

            tp_all(t,2) = nnz(tp_vec);
            fp_all(t,2) = nnz(fp_vec);
            fn_all(t,2) = nnz(fn_vec);
            sc_all(t,2) = nnz(pos_vec);
        end
    end

    %Let's speed up the easy calculations...
    res_untrimmed(:, 'spotCount') = array2table(uint32(sc_all(:,1)));
    res_untrimmed(:, 'true_pos') = array2table(uint32(tp_all(:,1)));
    res_untrimmed(:, 'false_pos') = array2table(uint32(fp_all(:,1)));
    res_untrimmed(:, 'false_neg') = array2table(uint32(fn_all(:,1)));

    recall = tp_all(:,1) ./ (tp_all(:,1) + fn_all(:,1));
    precision = tp_all(:,1) ./ (tp_all(:,1) + fp_all(:,1));
    fscores = (2 .* precision .* recall) ./ (precision + recall);
    pr_auc = RNAUtils.calculateAUC(recall, precision);
    peak_fscore = max(fscores, [], 'all');
    res_untrimmed(:, 'sensitivity') = array2table(recall);
    res_untrimmed(:, 'precision') = array2table(precision);
    res_untrimmed(:, 'fScore') = array2table(fscores);
    if any_trimmed
        res_trimmed(:, 'spotCount') = array2table(uint32(sc_all(:,2)));
        res_trimmed(:, 'true_pos') = array2table(uint32(tp_all(:,2)));
        res_trimmed(:, 'false_pos') = array2table(uint32(fp_all(:,2)));
        res_trimmed(:, 'false_neg') = array2table(uint32(fn_all(:,2)));

        recall = tp_all(:,2) ./ (tp_all(:,2) + fn_all(:,2));
        precision = tp_all(:,2) ./ (tp_all(:,2) + fp_all(:,2));
        fscores = (2 .* precision .* recall) ./ (precision + recall);
        pr_auc_trim = RNAUtils.calculateAUC(recall, precision);
        peak_fscore_trim = max(fscores, [], 'all');
        res_trimmed(:, 'sensitivity') = array2table(recall);
        res_trimmed(:, 'precision') = array2table(precision);
        res_trimmed(:, 'fScore') = array2table(fscores);
    end

    th_idx = 0;
    if th_val > 0
        th_idx = RNAUtils.findThresholdIndex(th_val, spot_table(:,1).');
    end

    %Save to output struct
    rstruct.performance = res_untrimmed;
    rstruct.pr_auc = pr_auc;
    rstruct.fscore_peak = peak_fscore;
    if th_idx > 0
        rstruct.fscore_autoth = res_untrimmed{th_idx, 'fScore'};
    end
    if any_trimmed
        rstruct.performance_trimmed = res_trimmed;
        rstruct.pr_auc_trimmed = pr_auc_trim;
        rstruct.fscore_peak_trimmed = peak_fscore_trim;
        if th_idx > 0
            rstruct.fscore_autoth_trimmed = res_trimmed{th_idx, 'fScore'};
        end
    end

end
