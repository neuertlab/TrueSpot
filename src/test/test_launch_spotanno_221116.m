%
%%  !! UPDATE TO YOUR BASE DIR
DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Settings ==========================
addpath('./core');

ImgName = 'mESC_loday_D1I15_Xist';
RefMode = true;
NewAnno = false;
JustLoad = false;

% ========================== Read Table ==========================

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

% ========================== Load Spots Anno ==========================

outstem = [DataDir replace(getTableValue(image_table, rec_row, 'OUTSTEM'), '/', filesep)];

spotsrun = RNASpotsRun.loadFrom(outstem);
if isempty(spotsrun)
    fprintf('Run with prefix "%s" was not found!\n', outstem);
    return;
end

%Update spotsrun paths
spotsrun.out_stem = outstem;
spotsrun.saveMe();

if NewAnno | ~RNA_Threshold_SpotSelector.selectorExists(outstem)
    selector = RNA_Threshold_SpotSelector;
    th_idx = round((spotsrun.t_max - spotsrun.t_min) ./ 2);
    selector = selector.initializeNew(outstem, th_idx, []);
    selector.z_min = 1;
    selector.z_max = spotsrun.idims_sample.z;
else
    selector = RNA_Threshold_SpotSelector.openSelector(outstem, true);
end

if JustLoad
    return;
end

selector.crosshair_color = [1.0 1.0 0.0];
if RefMode
    selector.launchRefSelectGUI();
else
    selector.launchGUI();
end

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end
