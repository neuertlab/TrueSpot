%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
% ========================== Constants ==========================

ImageName = 'scrna_E1R2I4_STL1';
%dbgcell = 19;
worker_count = 2;

% ========================== Load csv Table ==========================

InputTablePath = [ImgProcDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Find table entry ==========================

row_idx = 0;
entry_count = size(imgtbl,1);
for r = 1:entry_count
    myname = getTableValue(imgtbl, r, 'IMGNAME');
    if strcmp(myname, ImageName)
        row_idx = r;
        break;
    end
end
if row_idx < 1
    fprintf('Could not find image "%s" in table! Exiting...\n', ImageName);
    return;
end

% ========================== Run ==========================

hb_outstem = [ImgProcDir replace(getTableValue(imgtbl, row_idx, 'OUTSTEM'), '/', filesep)];
tif_path = [ImgDir replace(getTableValue(imgtbl, row_idx, 'IMAGEPATH'), '/', filesep)];
cs_dir = [ImgProcDir replace(getTableValue(imgtbl, row_idx, 'CELLSEG_DIR'), '/', filesep)];
cs_name = getTableValue(imgtbl, row_idx, 'CELLSEG_SFX');

[run_dir, ~, ~] = fileparts(hb_outstem);
run_path = [hb_outstem '_rnaspotsrun.mat'];

%View cell mask
% cellseg_path = [cs_dir filesep 'Lab_' cs_name '.mat'];
% if isfile(cellseg_path)
%     load(cellseg_path, 'cells');
%     rndr = CellsegDrawer.renderCellMask(cells, [], []);
%     fighandle = figure(1);
%     CellsegDrawer.drawCellsToFigure(rndr, 1, fighandle, [1 0 0]);
% end

%Actual run
% Main_RNAQuant('-runinfo', run_path, '-tif', tif_path, '-outdir', run_dir, ...
%     '-cellsegdir', cs_dir, '-cellsegname', cs_name, '-workers', worker_count,...
%     '-dbgcell', dbgcell);
Main_RNAQuant('-runinfo', run_path, '-tif', tif_path, '-outdir', run_dir, ...
    '-cellsegdir', cs_dir, '-cellsegname', cs_name, '-workers', worker_count,...
    '-noclouds');

% ========================== Helper functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end