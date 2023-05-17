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

FiltOutTablePath = [TblOutDir filesep 'stats_histo_f.tsv'];
RawOutTablePath = [TblOutDir filesep 'stats_histo_r.tsv'];

ImageTableCols = {'IMGNAME', 'MIN', 'MAX', 'PROPZERO', 'P_50', 'P_75', 'P_80',...
    'P_85', 'P_90', 'P_95', 'P_99', 'P_999'};

ImageTableColCount = size(ImageTableCols,2);

% ========================== Prep ==========================

if ~isfolder(TblOutDir)
    mkdir(TblOutDir);
end

FiltOutTableFile = fopen(FiltOutTablePath, 'w');
RawOutTableFile = fopen(RawOutTablePath, 'w');

for i = 1:ImageTableColCount
    if i ~= 1; fprintf(FiltOutTableFile, '\t'); end
    fprintf(FiltOutTableFile, ImageTableCols{i});

    if i ~= 1; fprintf(RawOutTableFile, '\t'); end
    fprintf(RawOutTableFile, ImageTableCols{i});
end
fprintf(FiltOutTableFile, '\n');
fprintf(RawOutTableFile, '\n');

% ========================== Do Things ==========================

for t = 1:ImgTableCount
    fprintf('Trying Table %s...\n', AllTablePaths{t});
    image_table = testutil_opentable(AllTablePaths{t});

    entry_count = size(image_table, 1);
    for r = 1:entry_count

        %Check if sim
        myname = getTableValue(image_table, r, 'IMGNAME');
        fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

        %Get res file path
        set_group_dir = getSetOutputDirName(myname);
        ResFilePath = [ResultsDir filesep set_group_dir filesep myname '_summary.mat'];

        if isfile(ResFilePath)
            load(ResFilePath, 'analysis');
            fprintf(FiltOutTableFile, '%s\t', myname);
            fprintf(RawOutTableFile, '%s\t', myname);

            if isfield(analysis, 'iprctile')
                fprintf(RawOutTableFile, '%d\t', analysis.imin);
                fprintf(RawOutTableFile, '%d\t', analysis.imax);
                fprintf(RawOutTableFile, 'NaN\t');
                for i = 1:7
                    fprintf(RawOutTableFile, '%d\t', analysis.iprctile(i,2));
                end
                fprintf(RawOutTableFile, '%d\n', analysis.i999);
            else
                for i = 1:10
                    fprintf(RawOutTableFile, 'NaN\t');
                end
                fprintf(RawOutTableFile, 'NaN\n');
            end

            if isfield(analysis, 'results_hb')
                if isfield(analysis.results_hb, 'iprctile_f')
                    fprintf(FiltOutTableFile, '%d\t', analysis.results_hb.imin_f);
                    fprintf(FiltOutTableFile, '%d\t', analysis.results_hb.imax_f);
                    fprintf(FiltOutTableFile, '%f\t', (1.0 - analysis.results_hb.fprop_nz));
                    for i = 1:7
                        fprintf(FiltOutTableFile, '%d\t', analysis.results_hb.iprctile_f(i,2));
                    end
                    fprintf(FiltOutTableFile, '%d\n', analysis.results_hb.i999_f);
                else
                    for i = 1:10
                        fprintf(FiltOutTableFile, 'NaN\t');
                    end
                    fprintf(FiltOutTableFile, 'NaN\n');
                end
            else
                for i = 1:10
                    fprintf(FiltOutTableFile, 'NaN\t');
                end
                fprintf(FiltOutTableFile, 'NaN\n');
            end
            

        else
            fprintf('Could not find summary file. Skipping...\n');
        end

    end
end

% ========================== Close Output Files ==========================

fclose(FiltOutTableFile);
fclose(RawOutTableFile);

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
