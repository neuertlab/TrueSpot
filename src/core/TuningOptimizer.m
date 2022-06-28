%
%%

function [window_size, mad_factor, max_scan_thresh] = TuningOptimizer(spotcount_table, bkg_spotcount_table, fscore_table)

    mad_min = -1.0;
    mad_max = 2.0;
    win_min = 5;
    win_max = 50;
    mad_start_ival = 0.25;
    mad_min_ival = 0.05;
    
    T = size(spotcount_table,1);
    
    %Max Scan
    frac_max_fscore = 0.5; %Call max scan at this fraction of max fscore (after peak)
    max_scan_thresh = 0;
    peak_fscore = 0.0;
    call_thresh = 0.0;
    for t = 1:T
        fscore = fscore_table(t,1);
        if fscore > peak_fscore
            peak_fscore = fscore;
            call_thresh = frac_max_fscore * peak_fscore;
            continue;
        end
        if fscore <= call_thresh
            max_scan_thresh = spotcount_table(t,1);
            break;
        end
    end
    
    %Window size & MAD factor (simultaneously)
    mad_values = [mad_min:mad_start_ival:mad_max];
    win_values = [win_min:5:win_max];
    tscan_table = zeros(mad_values, win_values);
    fscan_table = NaN(mad_values, win_values);
    
    

end