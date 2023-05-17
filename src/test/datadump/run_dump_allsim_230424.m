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
GroupPrefixes = {'sim_', 'simvar_', 'rsfish_sim_', 'simvarmass_'};
ImgTableCount = size(AllTablePaths, 2);
SimGroupCount = size(GroupPrefixes, 2);

OutTablePath = [TblOutDir filesep 'stats_sim.tsv'];
OutSpotsPathHB = [TblOutDir filesep 'spots_sim_hb.tsv'];
OutSpotsPathBF = [TblOutDir filesep 'spots_sim_bf.tsv'];
OutSpotsPathRS = [TblOutDir filesep 'spots_sim_rs.tsv'];
OutSpotsPathDB = [TblOutDir filesep 'spots_sim_db.tsv'];

ImageTableCols = {'IMGNAME', 'SPOTS_ACTUAL', 'BKG_LVL', 'AMP_LVL', 'BKG_VAR', 'AMP_VAR',...
    'CLUSTERING', 'HB_SPOTS', 'PRAUC_HB', 'FSPEAK_HB', 'HB_FSCORE', ...
    'HBTr_SPOTS', 'PRAUC_HBTr', 'FSPEAK_HBTr', 'HBTr_FSCORE',...
    'BF_SPOTS', 'PRAUC_BF', 'FSPEAK_BF', 'BF_FSCORE',...
    'PRAUC_RS', 'FSPEAK_RS', 'PRAUC_DB', 'FSPEAK_DB', 'FILT_PROP_ZERO'};
SpotTableCols = {'SPOT_IDX', 'IMGNAME', 'PASS_TH', 'XY_DIST', 'Z_DIST', 'XYZ_DIST',...
    'X_FIT', 'Y_FIT', 'Z_FIT', 'X_REF', 'Y_REF', 'Z_REF','ZFIT_Q'};

ImageTableColCount = size(ImageTableCols,2);
SpotTableColCount = size(SpotTableCols,2);

% ========================== Prep ==========================

if ~isfolder(TblOutDir)
    mkdir(TblOutDir);
end

OutTableFile = fopen(OutTablePath, 'w');
OutSpotsFileHB = fopen(OutSpotsPathHB, 'w');
OutSpotsFileBF = fopen(OutSpotsPathBF, 'w');
OutSpotsFileRS = fopen(OutSpotsPathRS, 'w');
OutSpotsFileDB = fopen(OutSpotsPathDB, 'w');

for i = 1:ImageTableColCount
    if i ~= 1; fprintf(OutTableFile, '\t'); end
    fprintf(OutTableFile, ImageTableCols{i});
end
fprintf(OutTableFile, '\n');

for i = 1:SpotTableColCount
    if i ~= 1; fprintf(OutSpotsFileHB, '\t'); end
    fprintf(OutSpotsFileHB, SpotTableCols{i});

    if i ~= 1; fprintf(OutSpotsFileBF, '\t'); end
    fprintf(OutSpotsFileBF, SpotTableCols{i});

    if i ~= 1; fprintf(OutSpotsFileRS, '\t'); end
    fprintf(OutSpotsFileRS, SpotTableCols{i});

    if i ~= 1; fprintf(OutSpotsFileDB, '\t'); end
    fprintf(OutSpotsFileDB, SpotTableCols{i});
end
fprintf(OutSpotsFileHB, '\n');
fprintf(OutSpotsFileBF, '\n');
fprintf(OutSpotsFileRS, '\n');
fprintf(OutSpotsFileDB, '\n');

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
            fprintf(OutTableFile, '%s\t', myname);

            if isfield(analysis, 'simkey')
                if isstruct(analysis.simkey)
                    fprintf(OutTableFile, '%d\t', size(analysis.simkey,2));
                else
                    fprintf(OutTableFile, '%d\t', size(analysis.simkey,1));
                end
            else
                fprintf(OutTableFile, 'NaN\t');
            end

            if isfield(analysis, 'simparam')
                if isfield(analysis.simparam, 'bg_level')
                    fprintf(OutTableFile, '%f\t', analysis.simparam.bg_level);
                else
                    fprintf(OutTableFile, 'NaN\t');
                end

                if isfield(analysis.simparam, 'amplitude_mean')
                    fprintf(OutTableFile, '%f\t', analysis.simparam.amplitude_mean);
                else
                    fprintf(OutTableFile, 'NaN\t');
                end

                if isfield(analysis.simparam, 'bg_var')
                    fprintf(OutTableFile, '%f\t', analysis.simparam.bg_var);
                else
                    fprintf(OutTableFile, 'NaN\t');
                end

                if isfield(analysis.simparam, 'amplitude_var')
                    fprintf(OutTableFile, '%f\t', analysis.simparam.amplitude_var);
                else
                    fprintf(OutTableFile, 'NaN\t');
                end

                if isfield(analysis.simparam, 'cluster_count')
                    fprintf(OutTableFile, '%d;%d\t', analysis.simparam.cluster_count, analysis.simparam.cluster_size);
                else
                    fprintf(OutTableFile, '<UNK>\t');
                end
            else
                fprintf(OutTableFile, 'NaN\tNaN\tNaN\tNaN\t<UNK>\t');
            end

            if isfield(analysis, 'results_hb')
                if isfield(analysis.results_hb, 'performance')
                    th_table = table2array(analysis.results_hb.performance(:,'thresholdValue'));
                    th_idx = RNAUtils.findThresholdIndex(analysis.results_hb.threshold, transpose(th_table));
                    
                    fprintf(OutTableFile, '%d\t', analysis.results_hb.performance{th_idx, 'spotCount'});
                    fprintf(OutTableFile, '%f\t', analysis.results_hb.pr_auc);
                    fprintf(OutTableFile, '%f\t', analysis.results_hb.fscore_peak);
                    
                    if isfield(analysis.results_hb, 'fscore_autoth')
                        fprintf(OutTableFile, '%f\t', analysis.results_hb.fscore_autoth);
                    else
                        fprintf(OutTableFile, 'NaN\t');
                    end

                    if isfield(analysis.results_hb, 'performance_trimmed')
                        fprintf(OutTableFile, '%d\t', analysis.results_hb.performance_trimmed{th_idx, 'spotCount'});
                        fprintf(OutTableFile, '%f\t', analysis.results_hb.pr_auc_trimmed);
                        fprintf(OutTableFile, '%f\t', analysis.results_hb.fscore_peak_trimmed);
                        if isfield(analysis.results_hb, 'fscore_autoth_trimmed')
                            fprintf(OutTableFile, '%f\t', analysis.results_hb.fscore_autoth_trimmed);
                        else
                            fprintf(OutTableFile, 'NaN\t');
                        end
                    else
                        fprintf(OutTableFile, '%d\t', analysis.results_hb.performance{th_idx, 'spotCount'});
                        fprintf(OutTableFile, '%f\t', analysis.results_hb.pr_auc);
                        fprintf(OutTableFile, '%f\t', analysis.results_hb.fscore_peak);
                        if isfield(analysis.results_hb, 'fscore_autoth')
                            fprintf(OutTableFile, '%f\t', analysis.results_hb.fscore_autoth);
                        else
                            fprintf(OutTableFile, 'NaN\t');
                        end
                    end
                else
                    fprintf(OutTableFile, 'NaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t');
                end

                writeSpotTable(OutSpotsFileHB, analysis.results_hb, analysis.simkey, myname);

            else
                fprintf(OutTableFile, 'NaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t');
            end

            if isfield(analysis, 'results_bf')
                if isfield(analysis.results_bf, 'performance')
                    th_table = table2array(analysis.results_bf.performance(:,'thresholdValue'));
                    th_idx = RNAUtils.findThresholdIndex(analysis.results_bf.threshold, transpose(th_table));
                    
                    fprintf(OutTableFile, '%d\t', analysis.results_bf.performance{th_idx, 'spotCount'});
                    fprintf(OutTableFile, '%f\t', analysis.results_bf.pr_auc);
                    fprintf(OutTableFile, '%f\t', analysis.results_bf.fscore_peak);
                    fprintf(OutTableFile, '%f\t', analysis.results_bf.fscore_autoth);
                else
                    fprintf(OutTableFile, 'NaN\tNaN\tNaN\tNaN\t');
                end

                writeSpotTable(OutSpotsFileBF, analysis.results_bf, analysis.simkey, myname);
            else
                fprintf(OutTableFile, 'NaN\tNaN\tNaN\tNaN\t');
            end

            if isfield(analysis, 'results_rs')
                if isfield(analysis.results_rs, 'performance')
                    fprintf(OutTableFile, '%f\t', analysis.results_rs.pr_auc);
                    fprintf(OutTableFile, '%f\t', analysis.results_rs.fscore_peak);
                else
                    fprintf(OutTableFile, 'NaN\tNaN\t');
                end

                writeSpotTable(OutSpotsFileRS, analysis.results_rs, analysis.simkey, myname);
            else
                fprintf(OutTableFile, 'NaN\tNaN\t');
            end

            if isfield(analysis, 'results_db')
                if isfield(analysis.results_db, 'performance')
                    fprintf(OutTableFile, '%f\t', analysis.results_db.pr_auc);
                    fprintf(OutTableFile, '%f\t', analysis.results_db.fscore_peak);
                else
                    fprintf(OutTableFile, 'NaN\tNaN\t');
                end
                writeSpotTable(OutSpotsFileDB, analysis.results_db, analysis.simkey, myname);
            else
                fprintf(OutTableFile, 'NaN\tNaN\t');
            end

            if isfield(analysis, 'results_hb')
                if isfield(analysis.results_hb, 'fprop_nz')
                    fprop_z = 1.0 - analysis.results_hb.fprop_nz;
                    fprintf(OutTableFile, '%f\n', fprop_z);
                else
                    fprintf(OutTableFile, 'NaN\n');
                end
            else
                fprintf(OutTableFile, 'NaN\n');
            end

        else
            fprintf('> Results file for %s could not be found. Skipping...\n', myname);
        end
    end

end

% ========================== Close Output Files ==========================

fclose(OutTableFile);
fclose(OutSpotsFileHB);
fclose(OutSpotsFileBF);
fclose(OutSpotsFileRS);
fclose(OutSpotsFileDB);

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
        dirname = groupname;
    end
end

function writeSpotTable(tablefile, res_struct, simkey, myname)
if isempty(res_struct); return; end
if isempty(tablefile); return; end

if isfield(res_struct, 'callset')
    rcount = size(res_struct.ref_call_map,2);
    for s = 1:rcount
        callidx = res_struct.ref_call_map(s);
        if callidx > 0
            fprintf(tablefile, '%d\t%s\t', callidx, myname);

            if isfield(res_struct, 'threshold')
                pass_th = res_struct.callset{callidx, 'dropout_thresh'};
                if pass_th >= res_struct.threshold
                    fprintf(tablefile, '1\t');
                else
                    fprintf(tablefile, '0\t');
                end
            else
                fprintf(tablefile, '1\t');
            end
            

            fprintf(tablefile, '%f\t', res_struct.callset{callidx, 'xydist_ref'});
            fprintf(tablefile, '%f\t', res_struct.callset{callidx, 'zdist_ref'});
            fprintf(tablefile, '%f\t', res_struct.callset{callidx, 'xyzdist_ref'});
            fprintf(tablefile, '%f\t', res_struct.callset{callidx, 'fit_x'});
            fprintf(tablefile, '%f\t', res_struct.callset{callidx, 'fit_y'});
            fprintf(tablefile, '%f\t', res_struct.callset{callidx, 'fit_z'});

            if ~isempty(simkey)
                if isstruct(simkey)
                    fprintf(tablefile, '%f\t', simkey(s).x);
                    fprintf(tablefile, '%f\t', simkey(s).y);
                    fprintf(tablefile, '%f\t', simkey(s).z);
                else
                    fprintf(tablefile, '%f\t', simkey(s,1));
                    fprintf(tablefile, '%f\t', simkey(s,2));
                    fprintf(tablefile, '%f\t', simkey(s,3));
                end
            else
                fprintf(tablefile, 'NaN\tNaN\tNaN\t');
            end
            if nnz(ismember(res_struct.callset.Properties.VariableNames, 'zfitq'))
                fprintf(tablefile, '%f\n', res_struct.callset{callidx, 'zfitq'});
            else
                fprintf(tablefile, 'NaN\n');
            end
        end
    end
end
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end