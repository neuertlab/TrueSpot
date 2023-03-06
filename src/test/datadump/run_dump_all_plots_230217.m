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

START_INDEX = 1;
END_INDEX = 38;

DUMP_SPOTCOUNTS = true;
DUMP_FSCORES = true;
DUMP_ROC = true;

% ========================== Other Paths ==========================

ImageDumpDir = [ImgProcDir filesep 'figures' filesep 'all_curves'];
DataFilePath = [BaseDir filesep 'data' filesep 'DataComposite.mat'];

spc_figdir = [ImageDumpDir filesep 'spot_count'];
fsc_figdir = [ImageDumpDir filesep 'f_score'];
roc_figdir = [ImageDumpDir filesep 'roc'];

if ~isfolder(spc_figdir); mkdir(spc_figdir); end
if ~isfolder(fsc_figdir); mkdir(fsc_figdir); end
if ~isfolder(roc_figdir); mkdir(roc_figdir); end

% ========================== Go through ==========================

load(DataFilePath, 'image_analyses');

for i = START_INDEX:END_INDEX

if DUMP_SPOTCOUNTS
    figpath = [spc_figdir filesep 'spc__' image_analyses(i).imgname '.png'];
    fig_handle = GenMultiTool_SpotPlot(image_analyses(i).analysis);
    saveas(fig_handle, figpath);
    close(fig_handle);
end

if DUMP_FSCORES
    figpath = [fsc_figdir filesep 'fsc__' image_analyses(i).imgname '.png'];
    fig_handle = GenMultiTool_FScorePlot(image_analyses(i).analysis);
    saveas(fig_handle, figpath);
    close(fig_handle);
end

if DUMP_ROC
    figpath = [roc_figdir filesep 'roc__' image_analyses(i).imgname '.png'];
    fig_handle = GenMultiTool_ROCPlot(image_analyses(i).analysis);
    saveas(fig_handle, figpath);
    close(fig_handle);
end

end

