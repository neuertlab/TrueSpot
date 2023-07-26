%
%%  !! UPDATE TO YOUR BASE DIR
DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Settings ==========================
addpath('./core');

ImgName = 'mESC_loday_D0I04_Xist';

% ========================== Read Table ==========================

%InputTablePath = [DataDir filesep 'test_images_simytc.csv'];
%InputTablePath = [DataDir filesep 'test_images_simvarmass.csv'];
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

rec_row = 0;
rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    if strcmp(iname, ImgName)
        rec_row = r;
        break;
    end
end

if rec_row < 1
    fprintf('Image with name %s could not be found!\n', ImgName);
    return;
end

% ========================== Load Image View Struct ==========================

outstem = [DataDir replace(getTableValue(image_table, rec_row, 'OUTSTEM'), '/', filesep)];

spotsrun = RNASpotsRun.loadFrom(outstem);
if isempty(spotsrun)
    fprintf('Run with prefix "%s" was not found!\n', outstem);
    return;
end

%Update spotsrun paths
spotsrun.out_stem = outstem;
spotsrun.saveMe();

imgstruct_path = [outstem '_imgviewstructs.mat'];

if isfile(imgstruct_path)
    load(imgstruct_path, 'my_images');

    figure(1);
    clf;
    myimg = my_images(1);
    imshow(myimg.image, [myimg.Lmin myimg.Lmax]);
    hold on;
    title('Max Projection - Filtered');
    impixelinfo;

    figure(2);
    clf;
    myimg = my_images(2);
    imshow(myimg.image, [myimg.Lmin myimg.Lmax]);
    hold on;
    title('Max Projection - Raw');
    impixelinfo;
else
    fprintf('Stored max projections not found at stem %s...', outstem);
    return;
end

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end
