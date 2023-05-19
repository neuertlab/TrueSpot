%
%%  !! UPDATE TO YOUR BASE DIR
%BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
BaseDir = 'D:\usr\bghos\labdat\imgproc';

%ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

TblOutDir = [BaseDir filesep 'tables'];
ResultsDir = [BaseDir filesep 'data' filesep 'results'];

addpath('./core');
addpath('./test');

% ========================== Constants ==========================

TablePath_Main = [BaseDir filesep 'test_images.csv'];
TablePath_Mass = [BaseDir filesep 'test_images_simvarmass.csv'];
TablePath_YTC = [BaseDir filesep 'test_images_simytc.csv'];

AllTablePaths = {TablePath_Main, TablePath_Mass, TablePath_YTC};
ImgTableCount = size(AllTablePaths, 2);

OutTablePath = [TblOutDir filesep 'thcheck.tsv'];

ImageTableCols = {'IMGNAME', 'PASS_TH_CHECK', 'HB_AUC', 'HB_FSCORE', 'LOG10_SCMAX', 'LOG10_SELTHSPOTS'};
ImageTableColCount = size(ImageTableCols,2);

% ========================== Prep ==========================

if ~isfolder(TblOutDir)
    mkdir(TblOutDir);
end

OutTableFile = fopen(OutTablePath, 'w');

for i = 1:ImageTableColCount
    if i ~= 1; fprintf(OutTableFile, '\t'); end
    fprintf(OutTableFile, ImageTableCols{i});
end
fprintf(OutTableFile, '\n');

% ========================== Do Things ==========================

for t = 1:ImgTableCount
    fprintf('Trying Table %s...\n', AllTablePaths{t});
    image_table = testutil_opentable(AllTablePaths{t});

    entry_count = size(image_table, 1);
    for r = 1:entry_count
        myname = getTableValue(image_table, r, 'IMGNAME');
        fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

        %Get res file path
        set_group_dir = getSetOutputDirName(myname);
        ResFilePath = [ResultsDir filesep set_group_dir filesep myname '_summary.mat'];

        if isfile(ResFilePath)
            load(ResFilePath, 'analysis');
            fprintf(OutTableFile, '%s\t', myname);
            
            if isfield(analysis, 'results_hb')
                if isfield(analysis.results_hb, 'performance')
                    T = size(analysis.results_hb.performance, 1);
                    spot_table = NaN(T,2);

                    spot_table(:,1) = analysis.results_hb.performance{:,'thresholdValue'};
                    spot_table(:,2) = analysis.results_hb.performance{:,'spotCount'};
                    passth = RNA_Threshold_Common.curveIsTestable(spot_table);
                    if passth
                        fprintf(OutTableFile, '1\t');
                    else
                        fprintf(OutTableFile, '0\t');
                    end
                    fprintf(OutTableFile, '%f\t%f\t', analysis.results_hb.pr_auc, analysis.results_hb.fscore_autoth);

                    scmax = max(spot_table(:,2), [], 'all', 'omitnan');
                    logscmax = log10(double(scmax));
                    fprintf(OutTableFile, '%f\t', logscmax);

                    thidx = RNAUtils.findThresholdIndex(analysis.results_hb.threshold, spot_table(:,1).');
                    scval = spot_table(thidx, 2);
                    logscval = log10(double(scval));
                    fprintf(OutTableFile, '%f\n', logscval);
                else
                    thres = analysis.results_hb.threshold_details;
                    T = size(analysis.results_hb.threshold_details.x, 1);
                    spot_table = NaN(T,2);
                    spot_table(:,1) = analysis.results_hb.threshold_details.x(:,1);
                    dropout_th = analysis.results_hb.callset{:,'dropout_thresh'};
                    for j = 1:T
                        thval = spot_table(j,1);
                        spot_table(j,2) = nnz(dropout_th >= thval);
                    end
                    passth = RNA_Threshold_Common.curveIsTestable(spot_table);
                    if passth
                        fprintf(OutTableFile, '1\t');
                    else
                        fprintf(OutTableFile, '0\t');
                    end
                    fprintf(OutTableFile, 'NaN\tNaN\t');

                    scmax = max(spot_table(:,2), [], 'all', 'omitnan');
                    logscmax = log10(double(scmax));
                    fprintf(OutTableFile, '%f\t', logscmax);

                    thidx = RNAUtils.findThresholdIndex(analysis.results_hb.threshold, spot_table(:,1).');
                    scval = spot_table(thidx, 2);
                    logscval = log10(double(scval));
                    fprintf(OutTableFile, '%f\n', logscval);
                end
            else
                fprintf(OutTableFile, '-1\tNaN\tNaN\n');
            end

        else
            fprintf('Could not find summary file. Skipping...\n');
        end

    end
end

% ========================== Close Output Files ==========================

fclose(OutTableFile);

% ========================== Helper Functions ==========================

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

