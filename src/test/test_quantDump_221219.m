%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
% ========================== Constants ==========================

GroupPrefix = 'sctc_E2R3_';

%Because the channels are switched (?), switch these too
% Rna_on_tmr = 2;
% Rna_on_cy5 = 8;

Rna_on_tmr = 8;
Rna_on_cy5 = 2;

% ========================== Load csv Table ==========================

InputTablePath = [ImgDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Prepare output table ==========================

report_path = [ImgDir filesep 'qu_report.tsv'];
report_file = fopen(report_path, 'w');

fprintf(report_file, 'IMAGENAME\t');
fprintf(report_file, 'CELL_COUNT\t');
fprintf(report_file, 'CELLS_ON\t'); %>2 txn
fprintf(report_file, 'CELLS_ON_PROP\t');
fprintf(report_file, 'MEAN_PER_CELL\t');
fprintf(report_file, 'MEAN_PER_CELL_NUC\t');
fprintf(report_file, 'STD_PER_CELL\t');
fprintf(report_file, 'STD_PER_CELL_NUC\t');
fprintf(report_file, 'MEAN_PER_ON_CELL\t');
fprintf(report_file, 'MEAN_PER_ON_CELL_NUC\t');
fprintf(report_file, 'STD_PER_ON_CELL\t');
fprintf(report_file, 'STD_PER_ON_CELL_NUC\t');
fprintf(report_file, 'PROBE\n');

% ========================== Iterate through table entries ==========================
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
    
    spot_counts = NaN(1,cell_count);
    nuc_counts = NaN(1,cell_count);
    for c = 1:cell_count
        spot_counts(c) = mycells(c).spotcount_total;
        nuc_counts(c) = mycells(c).spotcount_nuc;
    end
    
    rna_for_on = Rna_on_tmr;
    if contains(iname, 'CTT1')
    	rna_for_on = Rna_on_cy5;
    end
    
    on_cells = spot_counts > rna_for_on;
    on_count = nnz(on_cells);
    
    all_mean = nanmean(spot_counts,'all');
    all_mean_nuc = nanmean(nuc_counts,'all');
    all_std = nanstd(spot_counts,0,'all');
    all_std_nuc = nanstd(nuc_counts,0,'all');
    
    off_cells_index = find(~on_cells);
    spot_counts(off_cells_index) = NaN;
    nuc_counts(off_cells_index) = NaN;
    on_mean = nanmean(spot_counts,'all');
    on_mean_nuc = nanmean(nuc_counts,'all');
    on_std = nanstd(spot_counts,0,'all');
    on_std_nuc = nanstd(nuc_counts,0,'all');
    
    fprintf(report_file, "%s\t", iname);
    fprintf(report_file, "%d\t", cell_count);
    fprintf(report_file, "%d\t", on_count);
    fprintf(report_file, "%f\t", (on_count/cell_count));
    fprintf(report_file, "%f\t", all_mean);
    fprintf(report_file, "%f\t", all_mean_nuc);
    fprintf(report_file, "%f\t", all_std);
    fprintf(report_file, "%f\t", all_std_nuc);
    fprintf(report_file, "%f\t", on_mean);
    fprintf(report_file, "%f\t", on_mean_nuc);
    fprintf(report_file, "%f\t", on_std);
    fprintf(report_file, "%f\t", on_std_nuc);
    fprintf(report_file, "%s\n", getTableValue(imgtbl, i, 'PROBE'));
end

fclose(report_file);

% ========================== Helper functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end