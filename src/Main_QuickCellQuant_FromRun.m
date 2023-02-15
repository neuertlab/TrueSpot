%
%%
%For just taking the coord table and counting points by cell.

function Main_QuickCellQuant_FromRun(run_path_stem, cell_seg_path, fixed_th)

if nargin < 3
    fixed_th = 0;
end

addpath('./core');
spotsrun = RNASpotsRun.loadFrom(run_path_stem);

if isempty(spotsrun)
    fprintf('Main_QuickCellQuant_FromRun || ERROR: Run at %s could not be found!\n', run_path_stem);
    return;
end
spotsrun.out_stem = run_path_stem;

%Find threshold index
if fixed_th == 0
    if ~isempty(spotsrun.threshold_results)
        th_val = spotsrun.threshold_results.threshold;
        if th_val > 0
            th_x = transpose(spotsrun.threshold_results.x);
            th_idx = RNAUtils.findThresholdIndex(th_val,th_x);
            clear th_x;
        else
            fprintf('Main_QuickCellQuant_FromRun || Threshold could not be determined. Using lowest threshold...\n');
            th_idx = 1;
        end
    else
        th_val = spotsrun.intensity_threshold;
        if th_val > 0
            [~,spot_table] = spotsrun.loadSpotsTable();
            th_x = transpose(spot_table(:,1));
            th_idx = RNAUtils.findThresholdIndex(th_val,th_x);
            clear spot_table;
            clear th_x;
        else
            fprintf('Main_QuickCellQuant_FromRun || Threshold could not be determined. Using lowest threshold...\n');
            th_idx = 1;
        end
    end
else
    [~,spot_table] = spotsrun.loadSpotsTable();
    th_idx = RNAUtils.findThresholdIndex(fixed_th,transpose(spot_table(:,1)));
    clear spot_table;
end

%Get image dims
if ~isempty(spotsrun.idims_sample)
    dimstr = sprintf('(%d,%d,%d)', spotsrun.idims_sample.x, spotsrun.idims_sample.y, spotsrun.idims_sample.z);
else
    fprintf('Main_QuickCellQuant_FromRun || Image dimensions could not be determined. Exiting...\n');
    return;
end

%Run.
clear spotsrun;
Main_QuickCellQuant([run_path_stem '_coordTable.mat'], cell_seg_path, [run_path_stem '_quickCellQuant.mat'], th_idx, dimstr);

end