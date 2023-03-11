%
%%  !! UPDATE TO YOUR BASE DIR
DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Settings ==========================
addpath('./core');

ImgName = 'sim_mESC_histone_100x_1';

USE_HB = 0;
USE_BF = 1;
USE_BFNR = 2;

use_which = USE_BF;

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

if use_which == USE_HB
    outstem = [DataDir replace(getTableValue(image_table, rec_row, 'OUTSTEM'), '/', filesep)];
elseif use_which == USE_BF
    bf_out = getTableValue(image_table, rec_row, 'BIGFISH_OUTSTEM');
    bfrs_out = replace(bf_out, '/bigfish/', '/bigfish/_rescaled/');
    outstem = [DataDir replace(bfrs_out, '/', filesep)];
elseif use_which == USE_BFNR
    outstem = [DataDir replace(getTableValue(image_table, rec_row, 'BIGFISH_OUTSTEM'), '/', filesep)];
else
    fprintf('Not sure which run to use. Exiting...\n');
    return;
end

if ~RNA_Threshold_SpotSelector.selectorExists(outstem)
    fprintf('Spotanno does not exist for this run! Exiting...\n');
    return;
end

RNA_Threshold_SpotSelector.fixPosTable(outstem);

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end