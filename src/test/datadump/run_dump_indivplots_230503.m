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

% ========================== Constants ==========================

ROW_INDEX = 65;

DUMP_SPOTCOUNTS = true;
DUMP_FSCORES = true;
DUMP_PRC = true;

USE_TRIMMED = true;

COLOR_HB = [0.667 0.220 0.220];
COLOR_BF = [0.231 0.231 0.702]; %#3b3bb3
COLOR_RS = [0.318 0.541 0.318]; %#518a51
COLOR_DB = [0.700 0.700 0.000];

% ========================== Other Paths ==========================

ImageDumpDir = [ImgProcDir filesep 'figures' filesep 'all_curves'];

ResultsDir = [BaseDir filesep 'data' filesep 'results'];

% ========================== Read Table ==========================

%InputTablePath = [BaseDir filesep 'test_images_simytc.csv'];
InputTablePath = [BaseDir filesep 'test_images_simvarmass.csv'];
%InputTablePath = [BaseDir filesep 'test_images.csv'];

image_table = testutil_opentable(InputTablePath);

% ========================== Go through ==========================

%Find summary file.
myname = getTableValue(image_table, ROW_INDEX, 'IMGNAME');
set_group_dir = getSetOutputDirName(myname);
ResFilePath = [ResultsDir filesep set_group_dir filesep myname '_summary.mat'];
if isfile(ResFilePath)
    load(ResFilePath, 'analysis');
else
    return;
end

if DUMP_SPOTCOUNTS
    if isfield(analysis, 'results_hb')
        if isfield(analysis.results_hb, 'performance')

            if USE_TRIMMED
                if isfield(analysis.results_hb, 'performance_trimmed')
                    x = analysis.results_hb.performance_trimmed{:,'thresholdValue'};
                    y = double(analysis.results_hb.performance_trimmed{:,'spotCount'});
                else
                    x = analysis.results_hb.performance{:,'thresholdValue'};
                    y = double(analysis.results_hb.performance{:,'spotCount'});
                end
            else
                x = analysis.results_hb.performance{:,'thresholdValue'};
                y = double(analysis.results_hb.performance{:,'spotCount'});
            end
            y = log10(y);

            SpotPlots.renderLogSpotCountPlot(x, y,...
                COLOR_HB, analysis.results_hb.threshold,...
                analysis.results_hb.threshold_details, 1, []);
        end
    end

    if isfield(analysis, 'results_bf')
        if isfield(analysis.results_bf, 'performance')
            if USE_TRIMMED
                if isfield(analysis.results_bf, 'performance_trimmed')
                    x = analysis.results_bf.performance_trimmed{:,'thresholdValue'};
                    y = double(analysis.results_bf.performance_trimmed{:,'spotCount'});
                else
                    x = analysis.results_bf.performance{:,'thresholdValue'};
                    y = double(analysis.results_bf.performance{:,'spotCount'});
                end
            else
                x = analysis.results_bf.performance{:,'thresholdValue'};
                y = double(analysis.results_bf.performance{:,'spotCount'});
            end
            y = log10(y);

            SpotPlots.renderLogSpotCountPlot(x, y,...
                COLOR_BF, analysis.results_bf.threshold,...
                [], 2, []);
        end
    end

    if isfield(analysis, 'results_rs')
        if isfield(analysis.results_rs, 'performance')
            x = analysis.results_rs.performance{:,'thresholdValue'};
            y = double(analysis.results_rs.performance{:,'spotCount'});
            y = log10(y);

            SpotPlots.renderLogSpotCountPlot(x, y,...
                COLOR_RS, 0,...
                [], 3, []);
        end
    end

    if isfield(analysis, 'results_db')
        if isfield(analysis.results_db, 'performance')
            x = analysis.results_db.performance{:,'thresholdValue'};
            y = double(analysis.results_db.performance{:,'spotCount'});
            y = log10(y);

            SpotPlots.renderLogSpotCountPlot(x, y,...
                COLOR_DB, 0,...
                [], 4, []);
        end
    end
end

if DUMP_FSCORES
    if isfield(analysis, 'results_hb')
        if isfield(analysis.results_hb, 'performance')
            if USE_TRIMMED
                if isfield(analysis.results_hb, 'performance_trimmed')
                    x = analysis.results_hb.performance_trimmed{:,'thresholdValue'};
                    y = double(analysis.results_hb.performance_trimmed{:,'fScore'});
                else
                    x = analysis.results_hb.performance{:,'thresholdValue'};
                    y = double(analysis.results_hb.performance{:,'fScore'});
                end
            else
                x = analysis.results_hb.performance{:,'thresholdValue'};
                y = double(analysis.results_hb.performance{:,'fScore'});
            end

            SpotPlots.renderFScorePlot(x, y,...
                COLOR_HB, analysis.results_hb.threshold,...
                analysis.results_hb.threshold_details, 5, []);
        end
    end

    if isfield(analysis, 'results_bf')
        if isfield(analysis.results_bf, 'performance')
            if USE_TRIMMED
                if isfield(analysis.results_bf, 'performance_trimmed')
                    x = analysis.results_bf.performance_trimmed{:,'thresholdValue'};
                    y = double(analysis.results_bf.performance_trimmed{:,'fScore'});
                else
                    x = analysis.results_bf.performance{:,'thresholdValue'};
                    y = double(analysis.results_bf.performance{:,'fScore'});
                end
            else
                x = analysis.results_bf.performance{:,'thresholdValue'};
                y = double(analysis.results_bf.performance{:,'fScore'});
            end

            SpotPlots.renderFScorePlot(x, y,...
                COLOR_BF, analysis.results_bf.threshold,...
                [], 6, []);
        end
    end

    if isfield(analysis, 'results_rs')
        if isfield(analysis.results_rs, 'performance')
            x = analysis.results_rs.performance{:,'thresholdValue'};
            y = double(analysis.results_rs.performance{:,'fScore'});

            SpotPlots.renderFScorePlot(x, y,...
                COLOR_RS, 0,...
                [], 7, []);
        end
    end

    if isfield(analysis, 'results_db')
        if isfield(analysis.results_db, 'performance')
            x = analysis.results_db.performance{:,'thresholdValue'};
            y = double(analysis.results_db.performance{:,'fScore'});

            SpotPlots.renderFScorePlot(x, y,...
                COLOR_DB, 0,...
                [], 8, []);
        end
    end
end

if DUMP_PRC
    if isfield(analysis, 'results_hb')
        if isfield(analysis.results_hb, 'performance')
            if USE_TRIMMED
                if isfield(analysis.results_hb, 'performance_trimmed')
                    x = analysis.results_hb.performance_trimmed{:,'sensitivity'};
                    y = double(analysis.results_hb.performance_trimmed{:,'precision'});
                else
                    x = analysis.results_hb.performance{:,'sensitivity'};
                    y = double(analysis.results_hb.performance{:,'precision'});
                end
            else
                x = analysis.results_hb.performance{:,'sensitivity'};
                y = double(analysis.results_hb.performance{:,'precision'});
            end

            SpotPlots.renderPRPlot(x, y, COLOR_HB, 9, []);
        end
    end

    if isfield(analysis, 'results_bf')
        if isfield(analysis.results_bf, 'performance')
            if USE_TRIMMED
                if isfield(analysis.results_bf, 'performance_trimmed')
                    x = analysis.results_bf.performance_trimmed{:,'sensitivity'};
                    y = double(analysis.results_bf.performance_trimmed{:,'precision'});
                else
                    x = analysis.results_bf.performance{:,'sensitivity'};
                    y = double(analysis.results_bf.performance{:,'precision'});
                end
            else
                x = analysis.results_bf.performance{:,'sensitivity'};
                y = double(analysis.results_bf.performance{:,'precision'});
            end

            SpotPlots.renderPRPlot(x, y, COLOR_BF, 10, []);
        end
    end

    if isfield(analysis, 'results_rs')
        if isfield(analysis.results_rs, 'performance')
            x = analysis.results_rs.performance{:,'sensitivity'};
            y = double(analysis.results_rs.performance{:,'precision'});

            SpotPlots.renderPRPlot(x, y, COLOR_RS, 11, []);
        end
    end

    if isfield(analysis, 'results_db')
        if isfield(analysis.results_db, 'performance')
            x = analysis.results_db.performance{:,'sensitivity'};
            y = double(analysis.results_db.performance{:,'precision'});

            SpotPlots.renderPRPlot(x, y, COLOR_DB, 12, []);
        end
    end
end

% ========================== Helper Functions ==========================

function cleanupFormatting()
    set(gca,'FontSize',12);
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
