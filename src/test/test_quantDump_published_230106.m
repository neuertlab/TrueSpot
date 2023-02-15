%
%%
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
InputDir = [ImgDir '\data\ref\yeast_tc'];

% ========================== Constants ==========================

ExpNumber = 1;
RepNumber = 2;

Rna_on_tmr = 2;
Rna_on_cy5 = 8;

input_path = [InputDir filesep 'Results_Exp' num2str(ExpNumber) '_rep' num2str(RepNumber) '.mat'];

% ========================== Ready output ==========================

report_path_tmr = [ImgDir filesep 'qupub_report_tmr.tsv'];
report_path_cy5 = [ImgDir filesep 'qupub_report_cy5.tsv'];
report_file_tmr = fopen(report_path_tmr, 'w');
report_file_cy5 = fopen(report_path_cy5, 'w');

% ========================== Read file ==========================

if isfile(input_path)
    load(input_path, 'Cells', 'RNA_CY5', 'RNA_TMR', 'timepoints');
    header_string = 'Timepoint\tCellCount\tMeanTot\tStdTot\tMeanNuc\tStdNuc\tONCountTot\tONPropTot\tONCountNuc\tONPropNuc\tMeanONTot\tStdONTot\tMeanONNuc\tStdONNuc\n';
    fprintf(report_file_tmr, header_string);
    fprintf(report_file_cy5, header_string);
    
    tp_count = size(timepoints,2);
    cells_points = size(Cells,2);
    if (tp_count > cells_points); tp_count = cells_points; end
    for i = 1:tp_count
        fprintf(report_file_tmr, '%d\t%d\t', timepoints(1,i), Cells(1,i));
        fprintf(report_file_cy5, '%d\t%d\t', timepoints(1,i), Cells(1,i));
        
        %TMR channel
        fprintf(report_file_tmr, '%f\t', nanmean(RNA_TMR.tot(:,i), 'all'));
        fprintf(report_file_tmr, '%f\t', nanstd(RNA_TMR.tot(:,i), 0, 'all'));
        fprintf(report_file_tmr, '%f\t', nanmean(RNA_TMR.nuc(:,i), 'all'));
        fprintf(report_file_tmr, '%f\t', nanstd(RNA_TMR.nuc(:,i), 0, 'all'));
        
        on_cells_idx_tot = find(RNA_TMR.tot(:,i) > Rna_on_tmr);
        on_cells_idx_nuc = find(RNA_TMR.nuc(:,i) > Rna_on_tmr);
        on_count = 0;
        if ~isempty(on_cells_idx_tot)
            on_count = size(on_cells_idx_tot,1);
        end
        fprintf(report_file_tmr, '%d\t', on_count);
        fprintf(report_file_tmr, '%f\t', (on_count/Cells(1,i)));
        on_count = 0;
        if ~isempty(on_cells_idx_nuc)
            on_count = size(on_cells_idx_nuc,1);
        end
        fprintf(report_file_tmr, '%d\t', on_count);
        fprintf(report_file_tmr, '%f\t', (on_count/Cells(1,i)));
        fprintf(report_file_tmr, '%f\t', nanmean(RNA_TMR.tot(on_cells_idx_tot,i), 'all'));
        fprintf(report_file_tmr, '%f\t', nanstd(RNA_TMR.tot(on_cells_idx_tot,i), 0, 'all'));
        fprintf(report_file_tmr, '%f\t', nanmean(RNA_TMR.nuc(on_cells_idx_nuc,i), 'all'));
        fprintf(report_file_tmr, '%f\n', nanstd(RNA_TMR.nuc(on_cells_idx_nuc,i), 0, 'all'));
        
        %CY5 channel
        fprintf(report_file_cy5, '%f\t', nanmean(RNA_CY5.tot(:,i), 'all'));
        fprintf(report_file_cy5, '%f\t', nanstd(RNA_CY5.tot(:,i), 0, 'all'));
        fprintf(report_file_cy5, '%f\t', nanmean(RNA_CY5.nuc(:,i), 'all'));
        fprintf(report_file_cy5, '%f\t', nanstd(RNA_CY5.nuc(:,i), 0, 'all'));
        
        on_cells_idx_tot = find(RNA_CY5.tot(:,i) > Rna_on_cy5);
        on_cells_idx_nuc = find(RNA_CY5.nuc(:,i) > Rna_on_cy5);
        on_count = 0;
        if ~isempty(on_cells_idx_tot)
            on_count = size(on_cells_idx_tot,1);
        end
        fprintf(report_file_cy5, '%d\t', on_count);
        fprintf(report_file_cy5, '%f\t', (on_count/Cells(1,i)));
        on_count = 0;
        if ~isempty(on_cells_idx_nuc)
            on_count = size(on_cells_idx_nuc,1);
        end
        fprintf(report_file_cy5, '%d\t', on_count);
        fprintf(report_file_cy5, '%f\t', (on_count/Cells(1,i)));
        fprintf(report_file_cy5, '%f\t', nanmean(RNA_CY5.tot(on_cells_idx_tot,i), 'all'));
        fprintf(report_file_cy5, '%f\t', nanstd(RNA_CY5.tot(on_cells_idx_tot,i), 0, 'all'));
        fprintf(report_file_cy5, '%f\t', nanmean(RNA_CY5.nuc(on_cells_idx_nuc,i), 'all'));
        fprintf(report_file_cy5, '%f\n', nanstd(RNA_CY5.nuc(on_cells_idx_nuc,i), 0, 'all'));
        
    end
else
    fprintf('File %s could not be found!\n', input_path);
end

% ========================== Close output ==========================

fclose(report_file_tmr);
fclose(report_file_cy5);