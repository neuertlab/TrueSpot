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

START_INDEX = 45;
END_INDEX = 45;

DUMP_SPOTCOUNTS = true;
DUMP_FSCORES = true;
DUMP_PRC = true;

% ========================== Other Paths ==========================

ImageDumpDir = [ImgProcDir filesep 'figures' filesep 'all_curves'];
%DataFilePath = [BaseDir filesep 'data' filesep 'DataComposite.mat'];

spc_figdir = [ImageDumpDir filesep 'spot_count'];
fsc_figdir = [ImageDumpDir filesep 'f_score'];
prc_figdir = [ImageDumpDir filesep 'prc'];

if ~isfolder(spc_figdir); mkdir(spc_figdir); end
if ~isfolder(fsc_figdir); mkdir(fsc_figdir); end
if ~isfolder(prc_figdir); mkdir(prc_figdir); end

ResultsDir = [BaseDir filesep 'data' filesep 'results'];

% ========================== Read Table ==========================

%InputTablePath = [BaseDir filesep 'test_images_simytc.csv'];
%InputTablePath = [BaseDir filesep 'test_images_simvarmass.csv'];
InputTablePath = [BaseDir filesep 'test_images.csv'];

image_table = testutil_opentable(InputTablePath);

% ========================== Go through ==========================

%TODO use indiv summar files instead!
%load(DataFilePath, 'image_analyses');

for i = START_INDEX:END_INDEX

    %Find summary file.
    myname = getTableValue(image_table, i, 'IMGNAME');
    set_group_dir = getSetOutputDirName(myname);
    ResFilePath = [ResultsDir filesep set_group_dir filesep myname '_summary.mat'];
    if isfile(ResFilePath)
        load(ResFilePath, 'analysis');
    else
        continue;
    end

    if DUMP_SPOTCOUNTS
        figpath = [spc_figdir filesep 'spc__' myname '.png'];
        %fig_handle = GenMultiTool_SpotPlot(image_analyses(i).analysis);
        fig_handle = GenMultiTool_SpotPlot(analysis);
        saveas(fig_handle, figpath);
        close(fig_handle);
    end

    if DUMP_FSCORES
        figpath = [fsc_figdir filesep 'fsc__' myname '.png'];
        %fig_handle = GenMultiTool_FScorePlot(image_analyses(i).analysis);
        fig_handle = GenMultiTool_FScorePlot(analysis);
        saveas(fig_handle, figpath);
        close(fig_handle);
    end

    if DUMP_PRC
        figpath = [prc_figdir filesep 'prc__' myname '.png'];
        %fig_handle = GenMultiTool_ROCPlot(image_analyses(i).analysis);
        fig_handle = GenMultiTool_ROCPlot(analysis);
        saveas(fig_handle, figpath);
        close(fig_handle);
    end

end

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

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

