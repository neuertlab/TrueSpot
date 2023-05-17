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

DateDir = '20230517';
DateSuffix = '230517';
OutDir = [ImgProcDir filesep 'figures' filesep DateDir];

% ========================== Parameters ==========================

DO_HB = true;
DO_BF = true;
DO_RS = true;
DO_DB = true;

HB_TRIMMED = true;
MAX_ZEROPROP = 0;
SNR_MAX = 240;

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
    zpstr = sprintf('%02d', MAX_ZEROPROP * 100);
else
    zpstr = 'n';
end

bkg_lvl = simres_table{:, 'BKG_LVL'};
amp_lvl = simres_table{:, 'AMP_LVL'};
%snr = amp_lvl ./ bkg_lvl;

actual_spots = simres_table{:, 'SPOTS_ACTUAL'};
fzprop = simres_table{:, 'FILT_PROP_ZERO'};
bkg_var = simres_table{:, 'BKG_VAR'};
amp_var = simres_table{:, 'AMP_VAR'};
snr = amp_lvl ./ (bkg_lvl .* bkg_var);

%snr_unit = max(snr, [], 'all', 'omitnan') / 20;
snr_unit = SNR_MAX / 20;

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
        fig = doHeatmap(1, snr, pr_auc, snr_unit, 0.05);
        title('Simulated Image PR-AUC (Homebrew)');
        ylabel('PR-AUC');
        saveas(fig, [OutDir filesep 'zp' zpstr '_prauc_hb_heatmap_' DateSuffix '.svg']);
    end

    if DO_FSCORE
        fig = doHeatmap(2, snr, f_scores, snr_unit, 0.05);
        title('Simulated Image F-Scores (Homebrew)');
        ylabel('F-Score');
        saveas(fig, [OutDir filesep 'zp' zpstr '_fscore_hb_heatmap_' DateSuffix '.svg']);
    end

end

if DO_BF
    pr_auc = simres_table{:, 'PRAUC_BF'};
    f_scores = simres_table{:, 'BF_FSCORE'};
    spot_det = simres_table{:, 'BF_SPOTS'};

    if DO_PRAUC
        fig = doHeatmap(3, snr, pr_auc, snr_unit, 0.05);
        title('Simulated Image PR-AUC (BigFISH)');
        ylabel('PR-AUC');
        saveas(fig, [OutDir filesep 'zp' zpstr '_prauc_bf_heatmap_' DateSuffix '.svg']);
    end

    if DO_FSCORE
        fig = doHeatmap(4, snr, f_scores, snr_unit, 0.05);
        title('Simulated Image F-Scores (BigFISH)');
        ylabel('F-Score');
        saveas(fig, [OutDir filesep 'zp' zpstr '_fscore_bf_heatmap_' DateSuffix '.svg']);
    end
    
end

if DO_RS
    pr_auc = simres_table{:, 'PRAUC_RS'};

    if DO_PRAUC
        fig = doHeatmap(5, snr, pr_auc, snr_unit, 0.05);
        title('Simulated Image PR-AUC (RS-FISH)');
        ylabel('PR-AUC');
        saveas(fig, [OutDir filesep 'zp' zpstr '_prauc_rs_heatmap_' DateSuffix '.svg']);
    end
end

if DO_DB
    pr_auc = simres_table{:, 'PRAUC_DB'};

    if DO_PRAUC
        fig = doHeatmap(6, snr, pr_auc, snr_unit, 0.05);
        title('Simulated Image PR-AUC (DeepBlink)');
        ylabel('PR-AUC');
        saveas(fig, [OutDir filesep 'zp' zpstr '_prauc_db_heatmap_' DateSuffix '.svg']);
    end
end

% ========================== Helper functions ==========================

function fig_handle = doHeatmap(figno, x, y, x_unit, y_unit, x_boxes, y_boxes)

    if nargin < 6; x_boxes = 20; end
    if nargin < 6; y_boxes = 20; end
    
    x_max = x_boxes * x_unit;
    y_max = y_boxes * y_unit;

    x_bounds = [0:x_unit:x_max];
    y_bounds = [0:y_unit:y_max];

    countmtx = NaN(y_boxes, x_boxes);
    for i = 1:x_boxes
        x_lo = x_bounds(i);
        x_hi = x_bounds(i+1);
        boolx_lo = (x >= x_lo);
        if i == x_boxes
            boolx_hi = (x <= x_hi);
        else
            boolx_hi = (x < x_hi);
        end
        boolx = boolx_hi & boolx_lo;

        for j = 1:y_boxes
            y_lo = y_bounds(j);
            y_hi = y_bounds(j+1);
            booly_lo = (y >= y_lo);
            if i == y_boxes
                booly_hi = (y <= y_hi);
            else
                booly_hi = (y < y_hi);
            end
            booly = booly_hi & booly_lo;
            boolxy = boolx & booly;
            countmtx(j,i) = nnz(boolxy);
        end
    end

    xlbl = x_bounds(1:x_boxes);
    ylbl = y_bounds(1:y_boxes);

    %Invert y
    ylbl = flip(ylbl);
    countmtx = flip(countmtx, 1);

    fig_handle = figure(figno);
    clf;
    hm = heatmap(xlbl, ylbl, countmtx);
    hm.Colormap = turbo;
    hm.CellLabelColor = 'none';
    xlabel('SNR');
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end