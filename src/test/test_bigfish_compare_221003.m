%
%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%BaseDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Paths ==========================

InputDir = [BaseDir '\data\bigfish\mESC_4d\AF594\Tsix'];
RefStem = [BaseDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];
BFName = 'BIGFISH_Tsix-AF594';

% ========================== Conversion via Main ==========================
Main_Bigfish2Mat(InputDir, [InputDir filesep BFName]);

% ========================== Run Normal ==========================
% addpath('./core');
% addpath('./test');
% leave_fig = true;
% 
% success_bool = BigfishCompare.doBigfishCompare(RefStem, InputDir, BFName, leave_fig);
% if success_bool
%     fprintf("Success: Compare method returned true\n");
% else
%     fprintf("Error: Compare method returned false\n");
% end