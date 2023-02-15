%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Constants ==========================
addpath('./core');

TH_MIN = 10;
TH_MAX = 500;
Z_TRIM = 3;

% ========================== Load csv Table ==========================
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

GroupPrefix = 'mESC4d_Tsix-AF594';

% ========================== Do things ==========================

rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if ~startsWith(iname, GroupPrefix); continue; end
    tifpath = [ImgDir replace(getTableValue(image_table, r, 'IMAGEPATH'), '/', filesep)];
    outstem = [DataDir replace(getTableValue(image_table, r, 'OUTSTEM'), '/', filesep)];
    [outdir, ~, ~] = fileparts(outstem);
    
    outdir = [outdir filesep 'maxproj'];
    mkdir(outdir);
    
    total_ch = getTableValue(image_table, r, 'CH_TOTAL');
    trg_ch = getTableValue(image_table, r, 'CHANNEL');
    sens_levl = getTableValue(image_table, r, 'THRESH_SETTING');
    
    sens_keyword = '-sensitivity';
    if sens_levl > 3
        sens_levl = sens_levl - 3;
    elseif sens_levl < 3
        sens_keyword = '-specificity';
        sens_levl = 3 - sens_levl;
    else
        sens_levl = 0;
    end
    
    fprintf('Running "%s"...\n', iname);
    %Not bothering with metadata fields
    Main_RNASpots('-imgname', iname, '-tif', tifpath, '-outdir', outdir, '-chsamp', trg_ch, ...
        '-chtotal', total_ch, '-thmin', TH_MIN, '-thmax', TH_MAX, '-ztrim', Z_TRIM, ...
        sens_keyword, sens_levl,'-maxzproj','-verbose');
    
    %Delete the big files and leave just the spot and coord tables
    delfile_stem = [outdir filesep iname '_max_proj'];
    delfile = [delfile_stem '_deadpix.mat'];
    if isfile(delfile); dres = delete(delfile); end
    delfile = [delfile_stem '_imgviewstructs.mat'];
    if isfile(delfile); dres = delete(delfile); end
    delfile = [delfile_stem '_prefilteredIMG.mat'];
    if isfile(delfile); dres = delete(delfile); end
    
end

% ========================== Helper funcs ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end