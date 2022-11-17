%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';

% ========================== Constants ==========================

DETECT_THREADS = 2;
TARGET_MOL = 'lncRNA';
SPECIES = 'Mus musculus';

TH_MIN = 50;
TH_MAX = 70;
Z_TRIM = 5;

% ========================== Load csv Table ==========================

InputTablePath = [DataDir filesep 'test_images.csv'];
imgtbl = readtable(InputTablePath,'Delimiter',',','ReadVariableNames',true,'Format','%s%s%d%d%d%d%s%s%s%s%s%d%s%s%s%d%d%d%d%d%d%d%d%s');

SingleImgName = 'mESC4d_Tsix-AF594';
OverrideName = 'Tsix-AF594_IMG1';

% ========================== Find Record ==========================
addpath('./core');

rec_row = 0;
rec_count = size(imgtbl,1);
for r = 1:rec_count
    myname = getTableValue(imgtbl, r, 'IMGNAME');
    if strcmp(myname, SingleImgName); rec_row = r; break; end
end

if rec_row < 1
    fprintf("Couldn't find image: %s!\n", SingleImgName);
    return;
end

% ========================== Prepare Spotsrun ==========================

rnaspots_run = RNASpotsRun.initDefault();
rnaspots_run.img_name = OverrideName;
rnaspots_run.tif_path = [ImgDir replace(getTableValue(imgtbl, rec_row, 'IMAGEPATH'), '/', filesep)];

fieldstr = replace(getTableValue(imgtbl, rec_row, 'OUTSTEM'), '/', filesep);
fieldstr = [DataDir fieldstr];
[rnaspots_run.out_dir, ~, ~] = fileparts(fieldstr);

fieldstr = getTableValue(imgtbl, rec_row, 'CELLSEG_DIR');
if ~strcmp(fieldstr, '.')
    cellsegdir = [DataDir replace(fieldstr, '/', filesep)];
    fieldstr = getTableValue(imgtbl, rec_row, 'CELLSEG_SFX');
    rnaspots_run.cellseg_path = [cellsegdir filesep 'Lab_' fieldstr '.mat'];
else
    rnaspots_run.cellseg_path = [];
end

fieldstr = getTableValue(imgtbl, rec_row, 'CONTROL_TIF');
if ~strcmp(fieldstr, '.')
    rnaspots_run.ctrl_path = [ImgDir replace(fieldstr, '/', filesep)];
    rnaspots_run.ctrl_ch = imgtbl{rec_row, 'CHANNEL'};
    rnaspots_run.ctrl_chcount = imgtbl{rec_row, 'CH_TOTAL'};
else
    rnaspots_run.ctrl_path = [];
end

rnaspots_run.rna_ch = imgtbl{rec_row, 'CHANNEL'};
rnaspots_run.light_ch = imgtbl{rec_row, 'CH_LIGHT'};
rnaspots_run.total_ch = imgtbl{rec_row, 'CH_TOTAL'};

rnaspots_run.idims_voxel.x = imgtbl{rec_row, 'VOXEL_X'};
rnaspots_run.idims_voxel.y = imgtbl{rec_row, 'VOXEL_Y'};
rnaspots_run.idims_voxel.z = imgtbl{rec_row, 'VOXEL_Z'};
rnaspots_run.idims_expspot.x = imgtbl{rec_row, 'POINT_X'};
rnaspots_run.idims_expspot.y = imgtbl{rec_row, 'POINT_Y'};
rnaspots_run.idims_expspot.z = imgtbl{rec_row, 'POINT_Z'};

rnaspots_run.type_probe = getTableValue(imgtbl, rec_row, 'PROBE');
rnaspots_run.type_cell = getTableValue(imgtbl, rec_row, 'CELLTYPE');
rnaspots_run.type_target = getTableValue(imgtbl, rec_row, 'TARGET');
rnaspots_run.type_targetmol = TARGET_MOL;
rnaspots_run.type_species = SPECIES;

rnaspots_run.t_min = TH_MIN;
rnaspots_run.t_max = TH_MAX;
rnaspots_run.ztrim = Z_TRIM;

rnaspots_run = RNAThreshold.applyPreset(rnaspots_run, imgtbl{rec_row,'THRESH_SETTING'});
%rnaspots_run.saveMe();

RNA_Pipeline_Core(rnaspots_run, 2, [], true, DETECT_THREADS);

% ========================== Do run ==========================

% ========================== Helpers ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    val = val{1,1};
end