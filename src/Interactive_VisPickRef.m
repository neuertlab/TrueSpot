%
%%

addpath('./core');
addpath('./thirdparty');

% ------------------------ INPUT PATHS ------------------------

ResultsBasePath = 'D:\Users\hospelb\labdata\RNAFISH\Analysis';
ImagesBasePath = 'D:\Users\hospelb\labdata\RNAFISH\Images';
%ImagesBasePath = 'C:\Users\hospelb\labdata\imgproc\img';

%ResultsBasePath = 'D:\usr\bghos\labdat\RNAFISH\Analysis';
%ImagesBasePath = 'D:\usr\bghos\labdat\RNAFISH\Images';

AnalysisDir = [ResultsBasePath filesep 'SimerlyLab\40x\MEO\FemaleX_MEO_BSTam_40xstack_I1'];
RunPath = [AnalysisDir filesep 'CH2\FemaleX_MEO_BSTam_40xstack_I1_CH2_spotCall_rnaspotsrun.mat'];

%Overrides

useOverrides = true;
saveOverrides = false;
ImagePath = [ImagesBasePath filesep 'SimerlyLab\40x\MEO\FemaleX_MEO_BSTam_40xstack_I1.tif'];
CellSegPath = [AnalysisDir filesep 'CellSeg_FemaleX_MEO_BSTam_40xstack_I1.mat'];

% ------------------------ LOAD/UPDATE RUN ------------------------

spotsrun = RNASpotsRun.loadFrom(RunPath, useOverrides);

if useOverrides
    spotsrun.paths.img_path = ImagePath;
    spotsrun.paths.cellseg_path = CellSegPath;

    if saveOverrides
        spotsrun.saveMe();
    end
end

quantPath = [spotsrun.getFullOutStem() '_quantData.mat'];
if ~isfile(quantPath)
    quantPath = replace(quantPath, '_spotCall', '');
    if ~isfile(quantPath)
        fprintf('WARNING: Quant file not found at "%s"\n', quantPath);
        quantPath = [];
    end
end

% ------------------------ GUI LOAD ------------------------

mygui = VisInteractive;
mygui = mygui.initRefSelectView(spotsrun, quantPath);

mygui.cellLayerOn = false;
mygui.nucLayerOn = false;
mygui.quantFitLayerOn = false;
mygui.spotCircleLayerOn = true;

mygui = mygui.launchFigureGUI();
