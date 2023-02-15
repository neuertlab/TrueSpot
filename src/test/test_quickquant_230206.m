%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Load csv Table ==========================
addpath('./core');
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

ImageName = 'mESC4d_Tsix-AF594';
FixedThresholdVal = 0;

% ========================== Try ==========================

row = 0;

rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    if strcmp(iname, ImageName)
        row = r;
        break;
    end
end

if row < 1
    fprintf('Image %s not found in table!\n', ImageName);
    return;
end

%Process outstem
hb_outstem = [DataDir replace(getTableValue(image_table,row,'OUTSTEM'), '/', filesep)];

%Get coord table path
coord_tbl_path = [hb_outstem '_coordTable.mat'];

%Generate output path
out_path = [hb_outstem '_quickCellQuant.mat'];

%Get cell seg path
cell_seg_dir = [DataDir replace(getTableValue(image_table,row,'CELLSEG_DIR'), '/', filesep)];
cell_seg_basename = getTableValue(image_table,row,'CELLSEG_SFX');
cell_seg_path = [cell_seg_dir filesep 'Lab_' cell_seg_basename '.mat'];

%Get selected threshold (From spotsrun if not specified)
if FixedThresholdVal == 0
    %Find value and index from spotsrun threshold_results
    spotsrun = RNASpotsRun.loadFrom(hb_outstem);
    if ~isempty(spotsrun)
        if ~isempty(spotsrun.threshold_results)
            %Use that call and x table to find idx
            thval = spotsrun.threshold_results.threshold;
            xtbl = transpose(spotsrun.threshold_results.x);
            th_idx = RNAUtils.findThresholdIndex(thval, xtbl);
            clear xtbl;
        else
            %Try intensity_threshold and spots table. If either of these
            %are no good, just set to 1.
            spots_tbl_path = [hb_outstem '_spotTable.mat'];
            if (isfile(spots_tbl_path)) & (spotsrun.intensity_threshold > 0)
                load(spots_tbl_path, 'spot_table');
                th_idx = RNAUtils.findThresholdIndex(spotsrun.intensity_threshold, transpose(spot_table(:,1)));
                clear spot_table;
            else
                fprintf('Threshold for %s could not be determined. Using minimum threshold...\n', iname);
                th_idx = 1;
            end
        end
    else
        fprintf('Spotsrun for %s could not be found. Using minimum threshold...\n', iname);
        th_idx = 1;
    end
else
    %It's fixed. Just figure out the index.
    spots_tbl_path = [hb_outstem '_spotTable.mat'];
    if isfile(spots_tbl_path)
        load(spots_tbl_path, 'spot_table');
        th_idx = RNAUtils.findThresholdIndex(FixedThresholdVal, transpose(spot_table(:,1)));
        clear spot_table;
    else
        %Just use direct value as index, I guess...
        th_idx = FixedThresholdVal;
    end
end

%Image dims
idim_x = getTableValue(image_table,row,'IDIM_X');
idim_y = getTableValue(image_table,row,'IDIM_Y');
idim_z = getTableValue(image_table,row,'IDIM_Z');
dimstr = sprintf('(%d,%d,%d)', idim_x, idim_y, idim_z);

%Call cell quant function.
Main_QuickCellQuant(coord_tbl_path, cell_seg_path, out_path, th_idx, dimstr);

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

