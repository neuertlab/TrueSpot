%
%%  !! UPDATE TO YOUR BASE DIR
%BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
BaseDir = 'D:\usr\bghos\labdat\imgproc';

%ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

ResultsDir = [BaseDir filesep 'data' filesep 'results'];

addpath('./core');
addpath('./test');
addpath('./test/datadump');

% ========================== Constants ==========================

TablePath_Main = [BaseDir filesep 'test_images.csv'];
TablePath_Mass = [BaseDir filesep 'test_images_simvarmass.csv'];
TablePath_YTC = [BaseDir filesep 'test_images_simytc.csv'];

AllTablePaths = {TablePath_Main, TablePath_Mass, TablePath_YTC};
GroupPrefixes = {'sim_', 'simvar_', 'rsfish_sim_', 'simvarmass_'};
ImgTableCount = size(AllTablePaths, 2);
SimGroupCount = size(GroupPrefixes, 2);

delete_fields = {'cell_type', 'probe_target', 'probe', 'probe_target_type',...
    'species', 'last_ts', 'truthset_BH'};

% ========================== Do Things ==========================

for t = 1:ImgTableCount
    fprintf('Trying Table %s...\n', AllTablePaths{t});
    image_table = testutil_opentable(AllTablePaths{t});

    entry_count = size(image_table, 1);
    for r = 1:entry_count

        %Check if sim
        myname = getTableValue(image_table, r, 'IMGNAME');
        fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

        is_sim = false;
        for j = 1:SimGroupCount
            if startsWith(myname, GroupPrefixes{j})
                is_sim = true;
                break;
            end
        end

        if ~is_sim
            fprintf('> %s is not a sim image. Skipping...\n', myname);
            continue;
        end

        %Get res file path
        set_group_dir = getSetOutputDirName(myname);
        ResFilePath = [ResultsDir filesep set_group_dir filesep myname '_summary.mat'];

        if isfile(ResFilePath)
            load(ResFilePath, 'analysis');
        else
            fprintf('> Result file not found... Skipping...\n');
            continue;
        end

        %Remove from main struct
        rmfcount = size(delete_fields, 2);
        for i = 1:rmfcount
            if isfield(analysis, delete_fields{i})
                analysis = rmfield(analysis, delete_fields{i});
            end
        end

        %Go through results structs and move datasets back
        if isfield(analysis, 'results_hb')
            if isfield(analysis.results_hb, 'truthset_BH')
                analysis.results_hb = rmfield(analysis.results_hb, 'truthset_BH');
            end

            if isfield(analysis.results_hb, 'last_ts')
                analysis.results_hb = rmfield(analysis.results_hb, 'last_ts');
            end

            %Fix fields in callset
            if tableHasColumn(analysis.results_hb.callset, 'is_true_BH')
                analysis.results_hb.callset{:,'is_true'} = analysis.results_hb.callset{:,'is_true_BH'};
                analysis.results_hb.callset = removevars(analysis.results_hb.callset, 'is_true_BH');
            end
            if tableHasColumn(analysis.results_hb.callset, 'in_truth_region_BH')
                analysis.results_hb.callset{:,'in_truth_region'} = analysis.results_hb.callset{:,'in_truth_region_BH'};
                analysis.results_hb.callset = removevars(analysis.results_hb.callset, 'in_truth_region_BH');
            end

            %Update performance
            T = size(analysis.results_hb.performance, 1);
            spot_table = NaN(T,2);
            spot_table(:,1) = analysis.results_hb.performance{:,'thresholdValue'};
            spot_table(:,2) = analysis.results_hb.performance{:,'spotCount'};
            analysis.results_hb = runstats(analysis.results_hb, spot_table, analysis.results_hb.threshold);
        end

        if isfield(analysis, 'results_bf')
            if isfield(analysis.results_bf, 'truthset_BH')
                analysis.results_bf = rmfield(analysis.results_bf, 'truthset_BH');
            end

            if isfield(analysis.results_bf, 'last_ts')
                analysis.results_bf = rmfield(analysis.results_bf, 'last_ts');
            end

            %Fix fields in callset
            if isfield(analysis.results_bf, 'callset')
                if tableHasColumn(analysis.results_bf.callset, 'is_true_BH')
                    analysis.results_bf.callset{:,'is_true'} = analysis.results_bf.callset{:,'is_true_BH'};
                    analysis.results_bf.callset = removevars(analysis.results_bf.callset, 'is_true_BH');
                end
                if tableHasColumn(analysis.results_bf.callset, 'in_truth_region_BH')
                    analysis.results_bf.callset{:,'in_truth_region'} = analysis.results_bf.callset{:,'in_truth_region_BH'};
                    analysis.results_bf.callset = removevars(analysis.results_bf.callset, 'in_truth_region_BH');
                end

                %Update performance
                T = size(analysis.results_bf.performance, 1);
                spot_table = NaN(T,2);
                spot_table(:,1) = analysis.results_bf.performance{:,'thresholdValue'};
                spot_table(:,2) = analysis.results_bf.performance{:,'spotCount'};
                analysis.results_bf = runstats(analysis.results_bf, spot_table, analysis.results_bf.threshold);
            end
        end

        if isfield(analysis, 'results_rs')
            if isfield(analysis.results_rs, 'truthset_BH')
                analysis.results_rs = rmfield(analysis.results_rs, 'truthset_BH');
            end

            if isfield(analysis.results_rs, 'last_ts')
                analysis.results_rs = rmfield(analysis.results_rs, 'last_ts');
            end

            %Fix fields in callset
            if tableHasColumn(analysis.results_rs.callset, 'is_true_BH')
                analysis.results_rs.callset{:,'is_true'} = analysis.results_rs.callset{:,'is_true_BH'};
                analysis.results_rs.callset = removevars(analysis.results_rs.callset, 'is_true_BH');
            end
            if tableHasColumn(analysis.results_rs.callset, 'in_truth_region_BH')
                analysis.results_rs.callset{:,'in_truth_region'} = analysis.results_rs.callset{:,'in_truth_region_BH'};
                analysis.results_rs.callset = removevars(analysis.results_rs.callset, 'in_truth_region_BH');
            end

            %Update performance
            T = size(analysis.results_rs.performance, 1);
            spot_table = NaN(T,2);
            spot_table(:,1) = analysis.results_rs.performance{:,'thresholdValue'};
            spot_table(:,2) = analysis.results_rs.performance{:,'spotCount'};
            analysis.results_rs = runstats(analysis.results_rs, spot_table, 0);
        end

        if isfield(analysis, 'results_db')
            if isfield(analysis.results_db, 'truthset_BH')
                analysis.results_db = rmfield(analysis.results_db, 'truthset_BH');
            end

            if isfield(analysis.results_db, 'last_ts')
                analysis.results_db = rmfield(analysis.results_db, 'last_ts');
            end

            %Fix fields in callset
            if tableHasColumn(analysis.results_db.callset, 'is_true_BH')
                analysis.results_db.callset{:,'is_true'} = analysis.results_db.callset{:,'is_true_BH'};
                analysis.results_db.callset = removevars(analysis.results_db.callset, 'is_true_BH');
            end
            if tableHasColumn(analysis.results_db.callset, 'in_truth_region_BH')
                analysis.results_db.callset{:,'in_truth_region'} = analysis.results_db.callset{:,'in_truth_region_BH'};
                analysis.results_db.callset = removevars(analysis.results_db.callset, 'in_truth_region_BH');
            end

            %Update performance
            T = size(analysis.results_db.performance, 1);
            spot_table = NaN(T,2);
            spot_table(:,1) = analysis.results_db.performance{:,'thresholdValue'};
            spot_table(:,2) = analysis.results_db.performance{:,'spotCount'};
            analysis.results_db = runstats(analysis.results_db, spot_table, 0);
        end

        save(ResFilePath, 'analysis');
        clear analysis;
    end
end

% ========================== Helper Functions ==========================

function boolresult = tableHasColumn(table, colname)
    memcheck = ismember(table.Properties.VariableNames, colname);
    boolresult = (nnz(memcheck) > 0);
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