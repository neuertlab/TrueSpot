%
%%  !! UPDATE TO YOUR BASE DIR
DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Settings ==========================
addpath('./core');

ImgName = 'sctc_E1R1_15m_I2_CTT1';
RefMode = false;
NewAnno = true;
JustLoad = false;

LoadRescaled = true;

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

bfstem_raw = getTableValue(image_table, rec_row, 'BIGFISH_OUTSTEM');

if LoadRescaled
    bfstem_raw = replace(bfstem_raw, '/bigfish/', '/bigfish/_rescaled/');
end

outstem = [DataDir replace(bfstem_raw, '/', filesep)];
refstem = [DataDir replace(getTableValue(image_table, rec_row, 'OUTSTEM'), '/', filesep)];

%Get initial threshold.
[bfdir, ~, ~] = fileparts(outstem);
[zmin, zmax, bfthresh] = BigfishCompare.readSummaryTxt([bfdir filesep 'summary.txt']);
load([outstem '_spotTable.mat'], 'spot_table');
th_idx = RNAUtils.findThresholdIndex(bfthresh, transpose(spot_table(:,1)));

if NewAnno | ~RNA_Threshold_SpotSelector.selectorExists(outstem)
    selector = RNA_Threshold_SpotSelector;
    selector = selector.initializeNew(outstem, th_idx, spot_table(:,1));
    selector.z_min = zmin;
    selector.z_max = zmax;
    %selector.threshold_idx = th_idx;
else
    selector = RNA_Threshold_SpotSelector.openSelector(outstem, true);
    selector.threshold_idx = th_idx;
end

if JustLoad
    return;
end

%Load in image data from reference...
refanno = RNA_Threshold_SpotSelector.openSelector(refstem, true);
selector.img_structs = refanno.img_structs;
selector.imgdat_path = refanno.imgdat_path;
if isfile(selector.imgdat_path)
    load(selector.imgdat_path, 'img_filter');
    selector.loaded_ch = double(img_filter);
end
clear refanno;

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
