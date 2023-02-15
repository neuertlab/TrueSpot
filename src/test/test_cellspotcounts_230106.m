%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
% ========================== Constants ==========================

GroupPrefix = 'sctc_E2R2_';

% ========================== Load csv Table ==========================

InputTablePath = [ImgDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Prepare output ==========================

report_path = [ImgDir filesep 'sc_cell_report_counts.tsv'];
report_file = fopen(report_path, 'w');

fprintf(report_file, 'EXPERIMENT\t');
fprintf(report_file, 'REPLICATE\t');
fprintf(report_file, 'IMAGE\t');
fprintf(report_file, 'TIMEPOINT\t');
fprintf(report_file, 'PROBE\t');
fprintf(report_file, 'CELL_IDX\t');
fprintf(report_file, 'SPOTS_HB_COORDS\t');
fprintf(report_file, 'SPOTS_BF_COORDS\t');
fprintf(report_file, 'SPOTS_HB_QUANT\t');
fprintf(report_file, 'CELL_LOC\n');

% ========================== Go through images ==========================
entry_count = size(imgtbl,1);

for i = 1:entry_count
    
    iname = getTableValue(imgtbl, i, 'IMGNAME');
    if ~isempty(GroupPrefix)
        if ~startsWith(iname, GroupPrefix) 
            fprintf("Image %d of %d - %s not in requested group! Skipping...\n", i, entry_count, iname);
            continue; 
        end
    end
    
    mystem = [ImgDir replace(getTableValue(imgtbl, i, 'OUTSTEM'), '/', filesep)];
    [datdir, ~, ~] = fileparts(mystem);
    
    qdatpath = [datdir filesep iname '_quantData.mat'];
    if ~isfile(qdatpath)
        fprintf("Image %d of %d - quant data not found for %s! Skipping...\n", i, entry_count, iname);
        continue; 
    end
    
    load(qdatpath, 'quant_results');
    mycells = quant_results.cell_rna_data;
    cell_count = size(mycells,2);
    
    %Load HB coord table
    hb_coord_table = loadHBCoordTableForSelectedThresh(mystem);
    
    %Load BF coord table
    bf_stem = replace(getTableValue(imgtbl, i, 'BIGFISH_OUTSTEM'), '/bigfish/', '/bigfish/_rescaled/');
    bf_stem = [ImgDir replace(bf_stem, '/', filesep)];
    bf_coord_table = loadBFCoordTableForSelectedThresh(bf_stem);
     
    %Determinie time point from image name
    [expno, repno, timepoint, imgno] = getImageNameInfo(iname);
    
    %Get probe name
    probename = getTableValue(imgtbl, i, 'PROBE');
    
    %Loop through cells...
    for c = 1:cell_count
        this_cell = mycells(c);
        
        fprintf(report_file, '%d\t%d\t%d\t%d\t%s\t%d\t', expno, repno, imgno, timepoint, probename, c);
        hb_cell_coords = this_cell.getCoordsSubset(hb_coord_table);
        if isempty(hb_cell_coords)
            fprintf(report_file, '0\t');
        else
            fprintf(report_file, '%d\t', size(hb_cell_coords,1));
        end
        
        bf_cell_coords = this_cell.getCoordsSubset(bf_coord_table);
        if isempty(bf_cell_coords)
            fprintf(report_file, '0\t');
        else
            fprintf(report_file, '%d\t', size(bf_cell_coords,1));
        end
        
        if ~isempty(this_cell.spots)
            fprintf(report_file, '%d\t', size(this_cell.spots,2));
        else
            fprintf(report_file, '0\t');
        end
        
        fprintf(report_file, '[%d,%d]', this_cell.cell_loc.left, this_cell.cell_loc.right);
        fprintf(report_file, '[%d,%d]', this_cell.cell_loc.top, this_cell.cell_loc.bottom);
        fprintf(report_file, '[%d,%d]\n', this_cell.cell_loc.z_bottom, this_cell.cell_loc.z_top);
    end
    
end

% ========================== Close output ==========================

fclose(report_file);

% ========================== Helper functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

function coords = loadHBCoordTableForSelectedThresh(hb_stem)
    coords = [];
    spotsrun = RNASpotsRun.loadFrom(hb_stem);
    if isempty(spotsrun); return; end
    if isempty(spotsrun.threshold_results); return; end
    hb_coord_path = [hb_stem '_coordTable.mat'];
    if ~isfile(hb_coord_path); return; end
    load(hb_coord_path, 'coord_table');
    tth = spotsrun.threshold_results.x;
    th_idx = find(tth == spotsrun.threshold_results.threshold);
    if isempty(th_idx); return; end
    coords = coord_table{th_idx,1};
end

function coords = loadBFCoordTableForSelectedThresh(bf_stem)
    coords = [];
    bf_coord_path = [bf_stem '_coordTable.mat'];
    if ~isfile(bf_coord_path); return; end
    
    [bf_dir, ~, ~] = fileparts(bf_stem);
    summ_path = [bf_dir filesep 'summary.txt'];
    [~, ~, bfthresh] = BigfishCompare.readSummaryTxt(summ_path);
    load(bf_coord_path, 'coord_table');
    
    bf_st_path = [bf_stem '_spotTable.mat'];
    load(bf_st_path, 'spot_table');
    th_idx = find(spot_table(:,1) == bfthresh);
    if isempty(th_idx); return; end
    coords = coord_table{th_idx,1};
end

function [expno, repno, timepoint, imgno] = getImageNameInfo(iname)
    expno = NaN; repno = NaN; timepoint = NaN; imgno = NaN;
    iname_parts = split(iname, '_');
    pcount = size(iname_parts,1);
    for j = 1:pcount
        if endsWith(iname_parts{j,1}, 'm')
            timepoint = str2num(replace(iname_parts{j,1}, 'm', ''));
        end
        if startsWith(iname_parts{j,1}, 'E')
            parts2 = split(iname_parts{j,1}, 'R');
            expno = str2num(replace(parts2{1,1}, 'E', ''));
            repno = str2num(parts2{2,1});
        end
        if startsWith(iname_parts{j,1}, 'I')
            imgno = str2num(replace(iname_parts{j,1}, 'I', ''));
        end
    end
end
