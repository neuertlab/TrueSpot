%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

ScriptDir = 'C:\Users\hospelb.VUDS\Desktop';

ClusterWorkDir = '/scratch/hospelb/imgproc';
ClusterScriptsDir = '/scratch/hospelb/scripts';
ClusterPyenvDir = '/scratch/hospelb/pyvenv';

% ========================== Constants ==========================

TH_MIN = 10;
TH_MAX = 1000;
Z_TRIM = 1;
BF_SOBJSZ = 50;
BF_NUCSZ = 200;
BF_RESCALE = false;

% ========================== Load csv Table ==========================
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

ImagePrefix='scprotein_';

% ========================== Scan Table ==========================
addpath('./core');

%Generate a bf wrapper script and hb wrapper script command for every img
%source directory found.

