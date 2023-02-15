%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
% ========================== Constants ==========================

THR_NAMES = ["HB" "BF" "RS" "HBBF" "HBRS"];
THR_COUNT = size(THR_NAMES,2);

SKIP_FSCORES = false;

%GroupPrefix = [];
GroupPrefix = 'sctc_E2R3_';

% ========================== Load csv Table ==========================

InputTablePath = [ImgDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Prepare output table ==========================

report_path = [ImgDir filesep 'th_report.tsv'];
report_file = fopen(report_path, 'w');

fprintf(report_file, 'IMGNAME\t');
fprintf(report_file, 'CURVE_CLASS\t');
fprintf(report_file, 'HB_TH_PRESET\t');
fprintf(report_file, 'PEAK_FSCORE\t'); %HB only
fprintf(report_file, 'AVG_SPOTS\t');
fprintf(report_file, 'STD_SPOTS\t');
for i = 1:THR_COUNT
    fprintf(report_file, [char(THR_NAMES(i)) '_COUNT\t']);
end
for i = 1:THR_COUNT
    fprintf(report_file, [char(THR_NAMES(i)) '_STDEVS\t']);
end
for i = 1:THR_COUNT
    fprintf(report_file, [char(THR_NAMES(i)) '_TH\t']);
end
for i = 1:THR_COUNT
    fprintf(report_file, [char(THR_NAMES(i)) '_FSCORE\t']);
end
fprintf(report_file, 'PROBE\n');


% ========================== Iterate through table entries ==========================
entry_count = size(imgtbl,1);

for i = 1:entry_count
    
    iname = getTableValue(imgtbl, i, 'IMGNAME');
    if ~isempty(GroupPrefix)
        if ~startsWith(iname, GroupPrefix) 
            fprintf("Image %d of %d - %s not in requested group! Skipping...\n", i, entry_count, iname);
            continue; 
        end
    end

    mystem = replace(getTableValue(imgtbl, i, 'OUTSTEM'), '/', filesep);
    mystem = [ImgDir mystem];
    
%     skipme = getTableValue(imgtbl, i, 'SKIP_RERUN');
%     if skipme ~= 0
%         fprintf("Image %d of %d - Skip requested for %s! Skipping...\n", i, entry_count, mystem);
%         continue;
%     end
    
    spotsrun = RNASpotsRun.loadFrom(mystem);
    if isempty(spotsrun)
        fprintf("Image %d of %d - Run could not be found for %s! Skipping...\n", i, entry_count, mystem);
        continue;
    end
    
    if isempty(spotsrun.threshold_results)
        fprintf("Image %d of %d - Threshold results not found for %s! Skipping...\n", i, entry_count, mystem);
        continue;
    end
    
    if spotsrun.threshold_results.threshold < 1
        fprintf("Image %d of %d - Threshold results not valid for %s! Skipping...\n", i, entry_count, mystem);
        continue;
    end
    
    fprintf("Image %d of %d - Run found for %s. Loading... \n", i, entry_count, mystem);
    
    spotsrun.out_stem = mystem;
    %spotsrun = spotsrun.fixSpotcountTables();
    th_preset = getTableValue(imgtbl, i, 'THRESH_SETTING');
    probename = getTableValue(imgtbl, i, 'PROBE');
    
    all_spotcounts = NaN(1,5); %For putting spotcounts at each picked th
    
    %Get HB Info
    peak_fscore = NaN;
    thresh_hb = makeThStruct();
    thresh_hb.threshold = spotsrun.threshold_results.threshold;
    [~, spot_table] = spotsrun.loadSpotsTable();
    tidx = find(spot_table(:,1) == thresh_hb.threshold, 1);
    thresh_hb.spot_count = spot_table(tidx, 2);
    all_spotcounts(1,1) = thresh_hb.spot_count;
    
    if ~SKIP_FSCORES & RNA_Threshold_SpotSelector.refsetExists(mystem)
        fscores_vals = RNA_Threshold_SpotSelector.loadFScores(mystem, spotsrun.z_min_apply, spotsrun.z_max_apply);
        thresh_hb.fscore = fscores_vals(tidx, 1);
        peak_fscore = max(fscores_vals(:,1), [], 'all');
    else
        thresh_hb.fscore = NaN;
    end
    
    %Get BF Info
    bfstem_raw = getTableValue(imgtbl, i, 'BIGFISH_OUTSTEM');
    bfstem = [ImgDir replace(bfstem_raw, '/', filesep)];
    [bf_dir, ~, ~] = fileparts(bfstem);
    bf_st_path = [bfstem '_spotTable.mat'];
    thresh_bf = makeThStruct();
    thresh_hbbf = makeThStruct();
    if isfile(bf_st_path)
        [~, ~, bfthresh] = BigfishCompare.readSummaryTxt([bf_dir filesep 'summary.txt']);
        load(bf_st_path, 'spot_table');
        thresh_bf.threshold = bfthresh;
        tidx = find(spot_table(:,1) == thresh_bf.threshold, 1);
        if ~isempty(tidx)
            thresh_bf.spot_count = spot_table(tidx, 2);
            all_spotcounts(1,2) = thresh_bf.spot_count;
        end
        
        tidx2 = 0;
        try
            t_res = RNAThreshold.runSavedParameters(spotsrun, 0, spot_table, []);
            thresh_hbbf.threshold = t_res.threshold;
            tidx2 = find(spot_table(:,1) == thresh_hbbf.threshold, 1);
            thresh_hbbf.spot_count = spot_table(tidx2, 2);
            all_spotcounts(1,3) = thresh_hbbf.spot_count;
        catch EXCP
            fprintf("Error thresholding BF run... Skipping HBBF...\n");
        end
        
        %Check if fscores present...
        if ~SKIP_FSCORES
            fscores = RNA_Threshold_SpotSelector.loadFScores(bfstem);
            if ~isempty(fscores)
                if ~isempty(tidx)
                    thresh_bf.fscore = fscores(tidx, 1);
                end
                if tidx2 > 0
                    thresh_hbbf.fscore = fscores(tidx2, 1);
                end
            end
        end
    end
    
    %Get RS Info
    rsstem = replace(bfstem_raw, '/bigfish/', '/bigfish/_rescaled/');
    rsstem = [ImgDir replace(rsstem, '/', filesep)];
    [rs_dir, ~, ~] = fileparts(rsstem);
    rs_st_path = [rsstem '_spotTable.mat'];
    thresh_rs = makeThStruct();
    thresh_hbrs = makeThStruct();
    if isfile(rs_st_path)
        [~, ~, bfthresh] = BigfishCompare.readSummaryTxt([rs_dir filesep 'summary.txt']);
        load(rs_st_path, 'spot_table');
        thresh_rs.threshold = bfthresh;
        tidx = find(spot_table(:,1) == thresh_rs.threshold, 1);
        if ~isempty(tidx)
            thresh_rs.spot_count = spot_table(tidx, 2);
            all_spotcounts(1,4) = thresh_rs.spot_count;
        end
        
        tidx2 = 0;
        try
            t_res = RNAThreshold.runSavedParameters(spotsrun, 0, spot_table, []);
            thresh_hbrs.threshold = t_res.threshold;
            tidx2 = find(spot_table(:,1) == thresh_hbrs.threshold, 1);
            thresh_hbrs.spot_count = spot_table(tidx2, 2);
            all_spotcounts(1,5) = thresh_hbrs.spot_count;
        catch EXCP
            fprintf("Error thresholding Rescaled BF run... Skipping HBRS...\n");
        end
        
        %Check if fscores present...
        if ~SKIP_FSCORES
            fscores = RNA_Threshold_SpotSelector.loadFScores(rsstem);
            if ~isempty(fscores)
                if ~isempty(tidx)
                    thresh_rs.fscore = fscores(tidx, 1);
                end
                if tidx2 > 0
                    thresh_hbrs.fscore = fscores(tidx2, 1);
                end
            end
        end
    end
    
    %Put th info in array
    thstructs(5) = thresh_hbrs;
    thstructs(4) = thresh_hbbf;
    thstructs(3) = thresh_rs;
    thstructs(2) = thresh_bf;
    thstructs(1) = thresh_hb;
    
    %Calculate things
    spotcount_avg = nanmean(all_spotcounts, 'all');
    spotcount_std = nanstd(all_spotcounts, 0, 'all');
    for j = 1:5
        thstructs(j).stdevs = (thstructs(j).spot_count - spotcount_avg) / spotcount_std;
    end

    %Print
    fprintf(report_file, "%s\t.\t%d\t", iname, th_preset);
    fprintf(report_file, "%f\t%f\t%f\t", peak_fscore, spotcount_avg, spotcount_std);
    for j = 1:5
        fprintf(report_file, "%d\t", thstructs(j).spot_count);
    end
    for j = 1:5
        fprintf(report_file, "%f\t", thstructs(j).stdevs);
    end
    for j = 1:5
        fprintf(report_file, "%d\t", thstructs(j).threshold);
    end
    for j = 1:5
        fprintf(report_file, "%f\t", thstructs(j).fscore);
    end
    fprintf(report_file, "%s\n", probename);
end

fclose(report_file);

%%
function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

%%
function th_struct = makeThStruct()
    th_struct = struct('threshold', NaN);
    th_struct.spot_count = NaN;
    th_struct.fscore = NaN;
    th_struct.stdevs = NaN;
end
