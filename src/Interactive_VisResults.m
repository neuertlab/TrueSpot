%
%%

addpath('./core');
addpath('./thirdparty');

% ------------------------ INPUT PATHS ------------------------

ResultsBasePath = 'D:\Users\hospelb\labdata\RNAFISH\Analysis';
ImagesBasePath = 'C:\Users\hospelb\labdata\imgproc\img';

% AnalysisDir = [ResultsBasePath filesep 'JAJH_ProbeTest_20240827\JA_20240828_Jurkat_TrueSpot_hARF4_Cy5_1-2000_I1'];
% RunPath = [AnalysisDir filesep 'CH3\JA_20240828_Jurkat_TrueSpot_hARF4_Cy5_1-2000_I1_CH3_spotCall_rnaspotsrun.mat'];
AnalysisDir = [ResultsBasePath filesep 'sctc2017\E2R2\Exp2_rep2_10min_im3'];
RunPath = [AnalysisDir filesep 'CH2\Exp2_rep2_10min_im3_CH2_spotCall_rnaspotsrun.mat'];

%Overrides

useOverrides = true;
saveOverrides = false;
% ImagePath = [ImagesBasePath filesep 'JURKAT FISH\TrueSpots\JA_2024_08_28_Jurkat_TrueSpothARF4_Cy5_1-2000_1_MMStack_Pos0.ome.tif'];
% CellSegPath = [AnalysisDir filesep 'CellSeg_JA_20240828_Jurkat_TrueSpot_hARF4_Cy5_1-2000_I1.mat'];
ImagePath = [ImagesBasePath filesep 'yeast_full\Exp2_rep2\Exp2_rep2_10min_im3.tif'];
CellSegPath = [AnalysisDir filesep 'CellSeg_Exp2_rep2_10min_im3.mat'];

% ------------------------ INITIAL SETTINGS ------------------------

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
mygui = mygui.initCallView(spotsrun, quantPath);

mygui.cellLayerOn = true;
mygui.nucLayerOn = true;
mygui.quantFitLayerOn = false;
mygui.spotCircleLayerOn = true;

mygui = mygui.launchFigureGUI();
