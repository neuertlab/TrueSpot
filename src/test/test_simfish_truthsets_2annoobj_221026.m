%
%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Image Channels ==========================

% img_stems = cell(48,1);
% mtx_paths = cell(48,1);
% i = 1;
% 
% %----- Test
% img_stems{i,1} = [DataDir '\data\preprocess\sim\mESC_RNA_100x_1\RNA-AF594\mESC_RNA_100x_1_all_3d']; %Results stem
% mtx_paths{i,1} = [ImgDir '\img\sim\mESC_RNA_100x_1']; i = i+1; %Path to sim mat file

InputTablePath = [DataDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Cycle ==========================

addpath('./core');
%path_count = i-1;
path_count = size(imgtbl,1);

for j = 1:path_count
    %mystem = img_stems{j,1};
    %mypath = mtx_paths{j,1};
    
    myname = getTableValue(imgtbl, j, 'IMGNAME');
    mystem = [DataDir replace(getTableValue(imgtbl, j, 'OUTSTEM'), '/', filesep)];
    mypath = [ImgDir replace(getTableValue(imgtbl, j, 'IMAGEPATH'), '/', filesep)];

    fprintf("Now processing %d of %d (%s)...\n", j, path_count, myname);
    skipme = getTableValue(imgtbl, j, 'SKIP_RERUN');
    if skipme ~= 0
        fprintf("Image skip requested! Skipping...\n");
        continue;
    end
    
    if startsWith(myname, 'sim_')
        fprintf("Image marked as simulation! Checking for output...\n");
    else
        continue;
    end

    spotsrun = RNASpotsRun.loadFrom(mystem);
    if isempty(spotsrun)
        fprintf("Run info not found. Skipping.\n");
        continue;
    end
    
    spotsrun.out_stem = mystem;
    spotsrun.z_min = -1;
    spotsrun.z_max = -1;
    spotsrun = spotsrun.saveMe();

    spotsrun = RNA_Threshold_Common.cullRunThresholds(spotsrun, 2000000000);
    
    if isempty(spotsrun.idims_sample)
        load(mypath, 'imgdat');
        spotsrun.idims_sample = struct('x', 0, 'y', 0, 'z', 0);
        spotsrun.idims_sample.x = size(imgdat,2);
        spotsrun.idims_sample.y = size(imgdat,1);
        spotsrun.idims_sample.z = size(imgdat,3);
        clear imgdat;
    end
    spotsrun = spotsrun.updateZTrimParams();
    spotsrun = spotsrun.saveMe();
    
    selector = RNA_Threshold_SpotSelector;
    selector = selector.initializeNew(mystem, 100);

    %Get truthset
    load(mypath, 'key');
    rcount = size(key,2);
    selector.ref_coords = zeros(rcount,3);
    for i = 1:rcount
        selector.ref_coords(i,1) = key(i).x;
        selector.ref_coords(i,2) = key(i).y;
        selector.ref_coords(i,3) = key(i).z;
    end

    selector = selector.updateFTable();
    selector = selector.saveMe();
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end