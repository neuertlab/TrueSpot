%
%%  !! UPDATE TO YOUR BASE DIR
%BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
BaseDir = 'D:\usr\bghos\labdat\imgproc';

%ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
addpath('./test');

% ========================== Constants ==========================

SimResCSVPath = [BaseDir filesep 'sim_results.csv'];

% ========================== Parameters ==========================

DO_HB = true;
DO_BF = true;
DO_RS = true;
DO_DB = true;

HB_TRIMMED = true;
MAX_ZEROPROP = 0;

DO_PRAUC = true;
DO_FSCORE = true;
DO_SPOTSVS = true;

% ========================== Read Table ==========================

fmt_string = ['%s' repmat('%f', 1, 5) '%s' repmat('%f', 1, 17) '%s%s'];
simres_table = readtable(SimResCSVPath,'Delimiter',',','ReadVariableNames',true,'Format',...
    fmt_string);

% ========================== Do plot ==========================

%Filter..
if MAX_ZEROPROP > 0
    keeprows = find(~isnan(simres_table{:, 'FILT_PROP_ZERO'}));
    simres_table = simres_table(keeprows,:);
    keeprows = find(simres_table{:, 'FILT_PROP_ZERO'} <= MAX_ZEROPROP);
    simres_table = simres_table(keeprows,:);
end

bkg_lvl = simres_table{:, 'BKG_LVL'};
amp_lvl = simres_table{:, 'AMP_LVL'};
%snr = amp_lvl ./ bkg_lvl;

actual_spots = simres_table{:, 'SPOTS_ACTUAL'};
fzprop = simres_table{:, 'FILT_PROP_ZERO'};
bkg_var = simres_table{:, 'BKG_VAR'};
amp_var = simres_table{:, 'AMP_VAR'};
snr = amp_lvl ./ (bkg_lvl .* bkg_var);

if DO_HB
    if HB_TRIMMED
        pr_auc = simres_table{:, 'PRAUC_HBTr'};
        f_scores = simres_table{:, 'HBTr_FSCORE'};
        spot_det = simres_table{:, 'HBTr_SPOTS'};
    else
        pr_auc = simres_table{:, 'PRAUC_HB'};
        f_scores = simres_table{:, 'HB_FSCORE'};
        spot_det = simres_table{:, 'HB_SPOTS'};
    end

    if DO_PRAUC
        %MakeRedBlueScatterplot(snr, pr_auc, amp_var, bkg_var, 'SNR', 'PR-AUC', 1);
        MakeRedBlueScatterplot(snr, pr_auc, actual_spots, amp_var, 'SNR', 'PR-AUC', 1);
        ylim([0 1]);
        title('Simulated Image PR-AUC (Homebrew)');
        cleanupFormatting();
    end

    if DO_FSCORE
        %MakeRedBlueScatterplot(snr, f_scores, amp_var, bkg_var, 'SNR', 'F-Score', 2);
        MakeRedBlueScatterplot(snr, f_scores, actual_spots, amp_var, 'SNR', 'F-Score', 2);
        %MakeRedBlueScatterplot(snr, f_scores, actual_spots, fzprop, 'SNR', 'F-Score', 2);
        ylim([0 1]);
        title('Simulated Image F-Score (Homebrew)');
        cleanupFormatting();
    end

    if DO_SPOTSVS
        MakeRedBlueScatterplot(actual_spots, spot_det, snr, amp_var, 'Actual Spots', 'Detected Spots', 3);
        maxval = max(actual_spots, [], 'all', 'omitnan');
        maxval = max(maxval, max(spot_det, [], 'all', 'omitnan'));
        ylim([0 maxval]);
        xlim([0 maxval]);
        drawXeqYLine();
        title('Detected vs. Simulated Spots (Homebrew)');
        cleanupFormatting();
    end

end

if DO_BF
    pr_auc = simres_table{:, 'PRAUC_BF'};
    f_scores = simres_table{:, 'BF_FSCORE'};
    spot_det = simres_table{:, 'BF_SPOTS'};

    if DO_PRAUC
        %MakeRedBlueScatterplot(snr, pr_auc, amp_var, bkg_var, 'SNR', 'PR-AUC', 4);
        MakeRedBlueScatterplot(snr, pr_auc, actual_spots, amp_var, 'SNR', 'PR-AUC', 4);
        ylim([0 1]);
        title('Simulated Image PR-AUC (Big-FISH)');
        cleanupFormatting();
    end

    if DO_FSCORE
        %MakeRedBlueScatterplot(snr, f_scores, amp_var, bkg_var, 'SNR', 'F-Score', 5);
        MakeRedBlueScatterplot(snr, f_scores, actual_spots, amp_var, 'SNR', 'F-Score', 5);
        %MakeRedBlueScatterplot(snr, f_scores, actual_spots, fzprop, 'SNR', 'F-Score', 5);
        ylim([0 1]);
        title('Simulated Image F-Score (Big-FISH)');
        cleanupFormatting();
    end

    if DO_SPOTSVS
        MakeRedBlueScatterplot(actual_spots, spot_det, snr, amp_var, 'Actual Spots', 'Detected Spots', 6);
        maxval = max(actual_spots, [], 'all', 'omitnan');
        maxval = max(maxval, max(spot_det, [], 'all', 'omitnan'));
        ylim([0 maxval]);
        xlim([0 maxval]);
        drawXeqYLine();
        title('Detected vs. Simulated Spots (Big-FISH)');
        cleanupFormatting();
    end

end

if DO_RS
    pr_auc = simres_table{:, 'PRAUC_RS'};

    if DO_PRAUC
        %MakeRedBlueScatterplot(snr, pr_auc, amp_var, bkg_var, 'SNR', 'PR-AUC', 7);
        MakeRedBlueScatterplot(snr, pr_auc, actual_spots, amp_var, 'SNR', 'PR-AUC', 7);
        ylim([0 1]);
        title('Simulated Image PR-AUC (RS-FISH)');
        cleanupFormatting();
    end
end

if DO_DB
    pr_auc = simres_table{:, 'PRAUC_DB'};

    if DO_PRAUC
        %MakeRedBlueScatterplot(snr, pr_auc, amp_var, bkg_var, 'SNR', 'PR-AUC', 10);
        MakeRedBlueScatterplot(snr, pr_auc, actual_spots, amp_var, 'SNR', 'PR-AUC', 10);
        ylim([0 1]);
        title('Simulated Image PR-AUC (DeepBlink)');
        cleanupFormatting();
    end
end

if DO_HB & DO_BF & DO_FSCORE
    fs_hb = simres_table{:, 'HB_FSCORE'};
    fs_bf = simres_table{:, 'BF_FSCORE'};
    MakeRedBlueScatterplot(fs_bf, fs_hb, actual_spots, amp_var, 'F-Score (BF)', 'F-Score (HB)', 11);
    ylim([0 1]);
    xlim([0 1]);
    drawXeqYLine();
    title('F-Scores Homebrew vs. Big-FISH');
    cleanupFormatting();
end

% ========================== Helper functions ==========================

function cleanupFormatting()
    set(gca,'FontSize',12);
end

function drawXeqYLine()
    ax = gca;
    xmax = ax.XLim(2);
    ymax = ax.YLim(2);
    xymax = max(xmax, ymax);
    plot([0 xymax], [0 xymax], 'LineStyle', '-', 'LineWidth', 1, 'Color', 'black');
    hold on;
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end