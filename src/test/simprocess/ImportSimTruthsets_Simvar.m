%Imports ref sets from the original sim .mat files, generates a spotanno
%for the HB runs on these sim images, and inserts the truthset into that
%spotanno.
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Load csv Table ==========================
addpath('./core');
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

% ========================== Go through simvar entries ==========================

rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if ~startsWith(iname, 'simvar_'); continue; end
    
    %Load the source file
    tifpath_sfx = getTableValue(image_table, r, 'IMAGEPATH');
    matpath_sfx = replace(replace(tifpath_sfx, '/tif/', '/'), '.tif', '.mat');
    matpath = [ImgDir replace(matpath_sfx, '/', filesep)];
    load(matpath, 'key');
    
    %Load or generate spotanno
    outstem = [DataDir replace(getTableValue(image_table, r, 'OUTSTEM'), '/', filesep)];
    if ~RNA_Threshold_SpotSelector.selectorExists(outstem)
        selector = RNA_Threshold_SpotSelector;
        th_idx = 50;
        selector = selector.initializeNew(outstem, th_idx, []);
    else
        selector = RNA_Threshold_SpotSelector.openSelector(outstem, true);
    end
    
    %Copy key to ref
    scount = size(key,2);
    selector.ref_coords = uint16(zeros(scount,3));
    for s = 1:scount
        selector.ref_coords(s,1) = key(s).x;
        selector.ref_coords(s,2) = key(s).y;
        selector.ref_coords(s,3) = key(s).z;
    end
    selector.f_scores_dirty = true;
    
    %Save spotanno and clear
    selector.saveMe();
    clear selector;
    
end

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end