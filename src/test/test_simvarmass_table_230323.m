%
%% BASE DIR

ImgProcBaseDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Constants ==========================

addpath('./core');

UseDir = [ImgProcBaseDir filesep 'img' filesep 'simvarmass'];

VOX_XY = 65;
VOX_Z = 200;
POINT_XY = 95;
POINT_Z = 210;

% ========================== Scan Directory ==========================

dirContents = dir(UseDir);

filecount = size(dirContents,1);
for i = 1:filecount
    fileinfo = dirContents(i);
    if fileinfo.isdir; continue; end
    if ~(endsWith(fileinfo.name, '.mat')); continue; end

    imgname = replace(fileinfo.name, '.mat', '');
    fprintf('%s\t', imgname);
    fprintf('%s\t', ['/img/simvarmass/tif/' imgname '.tif']);
    fprintf('1\t0\t0\t1\t');
    fprintf('%s\t', ['/data/preprocess/simvarmass/' imgname '/' imgname '_all_3d']);
    fprintf('Sim\tSim\tSim\tSim\tSim\t');
    fprintf('%s\t', ['/data/bigfish/simvarmass/' imgname '/BIGFISH_' imgname]);
    fprintf('5\t.\t.\t.\t');
    fprintf('512\t512\t16\t');
    fprintf('%d\t%d\t%d\t', VOX_XY, VOX_XY, VOX_Z);
    fprintf('%d\t%d\t%d\t', POINT_XY, POINT_XY, POINT_Z);
    fprintf('0\t0\t0\t.');

    fprintf('\n');
end