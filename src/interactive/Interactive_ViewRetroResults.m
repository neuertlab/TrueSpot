%
%%

addpath('./core');
addpath('./thirdparty');

% ------------------------ INPUT PATHS ------------------------

ResultsBasePath = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImagesBasePath = 'C:\Users\hospelb\labdata\imgproc\img';

ImagePath = [ImagesBasePath filesep 'mESC_tc\mESC_loday_D1I15.tif'];
CellSegPath = [AnalysisDir filesep 'data\cell_seg\mESC_loday\*mESC_loday_D1I15'];
ResultsPath = [AnalysisDir filesep 'data\results\mescE\mESC_loday_D1I15_Xist_summary.mat'];

ChTotal = 4;
ChSample = 2;
ChNuc = 1;
ChLight = 4;

% ------------------------ GUI LOAD ------------------------

