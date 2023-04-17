%
%%  !! UPDATE TO YOUR BASE DIR
%BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
BaseDir = 'D:\usr\bghos\labdat\imgproc';

%ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
addpath('./test');
addpath('./test/datadump');

% ========================== Constants ==========================

START_INDEX = 1014;
END_INDEX = 1014;

DO_HOMEBREW = true;
DO_BIGFISH = true;
DO_RSFISH = true;

DO_TRUTHSET = true;

OutputDir = [BaseDir filesep 'data' filesep 'results'];

RS_TH_IVAL = 0.1/250;
SCRIPT_VER = 'v23.04.17.00';

% ========================== Load csv Table ==========================

%InputTablePath = [DataDir filesep 'test_images_simytc.csv'];
%InputTablePath = [DataDir filesep 'test_images_simvarmass.csv'];
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
    analysis.script_version = SCRIPT_VER;

    %Find truthset...
    ts_region = [];
    ref_coords = [];
    ref_coords = loadSimTruthset(image_table, r, ImgDir);
    if isempty(ref_coords)
        [ref_coords, ts_region] = loadExpTruthset(image_table, r, BaseDir);
    end

    if ~isempty(ts_region)
        analysis.truthset_region = ts_region;
    end

    %--------------------------------------- Truthset
    if DO_TRUTHSET
        if ~isempty(ref_coords)
            %If it's a sim image, go back and find the whole key and load
            %that in.
            if startsWith(myname, 'sim_') | startsWith(myname, 'simvar')
                key = [];
                if endsWith(srcpath, '.mat')
                    %Use this one directly.
                    load(srcpath, 'key');
                elseif endsWith(srcpath, '.tif')
                    %Find mat file and load from that.
                    matpath = replace(srcpath, [filesep 'tif' filesep], filesep);
                    matpath = replace(matpath, '.tif', '.mat');
                    load(matpath, 'key');
                end
                if ~isempty(key)
                    analysis.simkey = key;
                end
                clear key;
            elseif startsWith(myname, 'rsfish_sim')
                rsrefpath = replace(srcpath, '.tif', '.csv');

                if isfile(srcpath)
                    import_table = table2array(readtable(rsrefpath,'ReadVariableNames',false));
                    import_table = import_table + 1;
                    temp = import_table(:,2);
                    import_table(:,2) = import_table(:,1);
                    import_table(:,1) = temp;
                    analysis.simkey = import_table;
                end
            end
        end
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

    %--------------------------------------- Homebrew
    if DO_HOMEBREW
        %Look for run and load in coord and spot tables.
        spotsrun = RNASpotsRun.loadFrom(hb_stem);
        if ~isempty(spotsrun)
            spotsrun.out_stem = hb_stem;
            [~, spot_table] = spotsrun.loadSpotsTable();
            [~, coord_table] = spotsrun.loadCoordinateTable();

            %Apply filter to image.
            [IMG_filtered] = RNA_Threshold_SpotDetector.run_spot_detection_pre(my_image, '.', true, 7, false);

            %Do table transfer
            call_table = RNACoords.convertOldCoordTable(spot_table, coord_table, IMG_filtered, my_image, 2);
            init_call_count = size(call_table,1);

            gaussrad = spotsrun.dtune_gaussrad;
            if gaussrad < 1
                gaussrad = 7;
            end

            %TF call, if applicable
            if ~isempty(ref_coords)
                %Determine minimum threhsold for snapping...
                snapminth = 25;
                if spotsrun.intensity_threshold > 0
                    if spotsrun.intensity_threshold < snapminth
                        snapminth = spotsrun.intensity_threshold;
                    end
                end
                [call_table, ~] = RNACoords.updateTFCalls(call_table, ref_coords, 4, 2, snapminth);
                full_call_count = size(call_table, 1);

                %if fnegs were added, get the intensity values for those
                %Also need to update the 1D coords
                if full_call_count > init_call_count
                    addst = init_call_count + 1;
                    added = full_call_count;
                    x = table2array(call_table(addst:added,'isnap_x'));
                    y = table2array(call_table(addst:added,'isnap_y'));
                    z = table2array(call_table(addst:added,'isnap_z'));
                    c1d = sub2ind([Y X Z], y, x, z);
                    call_table(addst:added,'coord_1d') = array2table(uint32(c1d));
                    call_table(addst:added,'intensity_f') = array2table(single(IMG_filtered(c1d)));
                    call_table(addst:added,'intensity') = array2table(single(my_image(c1d)));

                    %Flag any fnegs that are outside detection area (trimmed
                    %out)
                    fneg_count = added - addst + 1;
                    fneg_trimmed = false(fneg_count,1);
                    fneg_trimmed = or(fneg_trimmed, x <= gaussrad);
                    fneg_trimmed = or(fneg_trimmed, x > (X - gaussrad));
                    fneg_trimmed = or(fneg_trimmed, y <= gaussrad);
                    fneg_trimmed = or(fneg_trimmed, y > (Y - gaussrad));
                    fneg_trimmed = or(fneg_trimmed, z < spotsrun.z_min_apply);
                    fneg_trimmed = or(fneg_trimmed, z > spotsrun.z_max_apply);
                    call_table(addst:added,'is_trimmed_out') = array2table(fneg_trimmed);
                    clear addst added x y z fneg_trimmed c1d
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

            end
            

            %TODO Fit matching, if applicable



            %Save updated results
            if ~isfield(analysis, 'results_hb')
                analysis.results_hb = struct('callset', table.empty());
            end
            analysis.results_hb.callset = call_table;
            analysis.results_hb.threshold = spotsrun.intensity_threshold;
            analysis.results_hb.threshold_details = spotsrun.threshold_results;
            analysis.results_hb.gaussrad_xy = gaussrad;
            analysis.results_hb.x_min = gaussrad;
            analysis.results_hb.x_max = X - gaussrad;
            analysis.results_hb.y_min = gaussrad;
            analysis.results_hb.y_max = Y - gaussrad;
            analysis.results_hb.z_min = spotsrun.z_min_apply;
            analysis.results_hb.z_max = spotsrun.z_max_apply;
            analysis.results_hb.th_scan_min = spotsrun.t_min;
            analysis.results_hb.th_scan_max = spotsrun.t_max;

            if ~isempty(ref_coords)
                %Calculate performance metrics
                analysis.results_hb = runstats(analysis.results_hb, spot_table, spotsrun.intensity_threshold);
            end

            clear call_table coord_table spot_table IMG_filtered gaussrad;
        else
            fprintf('ERROR: Could not find HB run for %s!\n', myname);
        end
        clear spotsrun;
    end

    %--------------------------------------- BigFISH
    if DO_BIGFISH
        bf_stem_base = replace(getTableValue(image_table, r, 'BIGFISH_OUTSTEM'), '/bigfish/', '/bigfish/_rescaled/');
        bf_stem = [BaseDir replace(bf_stem_base, '/', filesep)];

        bf_ct_path = [bf_stem '_coordTable.mat'];
        if isfile(bf_ct_path)
            [bf_dir, ~, ~] = fileparts(bf_stem);
            summary_path = [bf_dir filesep 'summary.txt'];
            bf_st_path = [bf_stem '_spotTable.mat'];

            load(bf_ct_path, 'coord_table');
            load(bf_st_path, 'spot_table');
            [zmin, zmax, bfthresh] = BigfishCompare.readSummaryTxt(summary_path);

            call_table = RNACoords.convertOldCoordTable(spot_table, coord_table, [], my_image, 1);
            init_call_count = size(call_table,1);

            if ~isempty(ref_coords)
                if bfthresh > 0
                    if bfthresh < snapminth
                        snapminth = bfthresh;
                    end
                end
                [call_table, ~] = RNACoords.updateTFCalls(call_table, ref_coords, 4, 2, snapminth);
                full_call_count = size(call_table, 1);

                if full_call_count > init_call_count
                    addst = init_call_count + 1;
                    added = full_call_count;
                    x = table2array(call_table(addst:added,'isnap_x'));
                    y = table2array(call_table(addst:added,'isnap_y'));
                    z = table2array(call_table(addst:added,'isnap_z'));
                    c1d = sub2ind([Y X Z], y, x, z);
                    call_table(addst:added,'coord_1d') = array2table(uint32(c1d));
                    call_table(addst:added,'intensity') = array2table(single(my_image(c1d)));

                    %Flag any fnegs that are outside detection area (trimmed
                    %out)
                    fneg_count = added - addst + 1;
                    fneg_trimmed = false(fneg_count,1);
                    fneg_trimmed = or(fneg_trimmed, z < zmin);
                    fneg_trimmed = or(fneg_trimmed, z > zmax);
                    call_table(addst:added,'is_trimmed_out') = array2table(fneg_trimmed);
                    clear addst added x y z fneg_trimmed c1d
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
            end

            %Fit matching TODO
            
            if ~isfield(analysis, 'results_bf')
                analysis.results_bf = struct('callset', table.empty());
            end
            analysis.results_bf.callset = call_table;
            analysis.results_bf.threshold = bfthresh;
            analysis.results_bf.z_min = zmin;
            analysis.results_bf.z_max = zmax;

            if ~isempty(ref_coords)
                %Calculate performance metrics
                analysis.results_bf = runstats(analysis.results_bf, spot_table, bfthresh);
            end

            clear call_table coord_table spot_table zmin zmax bfthresh;
        end
    end

    %--------------------------------------- RSFISH
    if DO_RSFISH
        rsdb_dir_ext = getRSDBGroupOutputDir(myname);
        if ~isempty(rsdb_dir_ext)
            rs_stem = [BaseDir filesep 'data' filesep 'rsfish' rsdb_dir_ext myname filesep 'RSFISH_' myname];
            [rs_dir, ~, ~] = fileparts(rs_stem);
            coord_table_path = [rs_stem '_coordTable.mat'];
            spot_table_path = [rs_stem '_spotTable.mat'];

            if isfile(coord_table_path)
                load(coord_table_path, 'coord_table');
                load(spot_table_path, 'spot_table');

                spot_table(:,1) = spot_table(:,1) .* RS_TH_IVAL;

                call_table = RNACoords.convertOldCoordTable(spot_table, coord_table, [], my_image, 1);
                init_call_count = size(call_table,1);

                if ~isempty(ref_coords)
                    [call_table, ~] = RNACoords.updateTFCalls(call_table, ref_coords, 4, 2, RS_TH_IVAL);
                    full_call_count = size(call_table, 1);

                    if full_call_count > init_call_count
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
                end

                %Fit matching TODO
            
                if ~isfield(analysis, 'results_rs')
                    analysis.results_rs = struct('callset', table.empty());
                end
                analysis.results_rs.callset = call_table;

                if ~isempty(ref_coords)
                    %Calculate performance metrics
                    analysis.results_rs = runstats(analysis.results_rs, spot_table, 0);
                end

            else
                fprintf('ERROR: Could not find RSFISH run for %s!\n', myname);
            end
            clear coord_table_path spot_table_path coord_table spot_table
        end
    end

    %--------------------------------------- Save
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
        dirname = groupname;
    end
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
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
%TODO
ref_coords = [];
valid_range = struct('x0', 0, 'x1', 0, 'y0', 0, 'y1', 0, 'z0', 0, 'z1', 0);
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
