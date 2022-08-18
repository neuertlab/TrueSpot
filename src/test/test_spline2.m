%
%%
%For histones - run the older dataset - more homogeneous


%Try three part spline? (Current algorithm may not work)
%Try matlab cubic spline functions?
%Try different criterion for "best" fits? (seg r2 weighted WAY in favor of
%   initial downslope)

BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc'; %workstation

RefStem = [BaseDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];

madf = -0.5;
winsz = 20;

addpath('./core');

spotsrun = RNASpotsRun.loadFrom(RefStem);
[spotsrun, spot_counts_ref] = spotsrun.loadSpotsTable();

[thresh_hb, win_out_hb, score_thresh_hb, scanst_hb] = RNA_Threshold_Common.estimateThreshold(spot_counts_ref, [], winsz, 0.5, madf);

%Trim everything left of the max.
ptct = size(win_out_hb,1);
[~,maxidx] = max(win_out_hb);
win_out_trimmed = win_out_hb(maxidx:ptct);
trimmed_points = size(win_out_trimmed,1);

%Prepare data table
mydata = NaN(trimmed_points, 2);
mydata(:,2) = win_out_trimmed(:,1);
mydata(:,1) = spot_counts_ref(maxidx:ptct,1);

%Try spline
linfit = Seglr2.fitTo(mydata);

%Plot
b1 = linfit.left.yintr;
b2 = linfit.right.yintr;
m1 = linfit.left.slope;
m2 = linfit.right.slope;

l1x1 = 0;
l1x2 = mydata(linfit.break_index,1);
l2x1 = l1x2;
l2x2 = spotsrun.t_max;
l1y1 = b1;
l1y2 = b1 + (l1x2 * m1);
l2y1 = b2 + (l2x1 * m2);
l2y2 = b2 + (l2x2 * m2);

figure(505);
hold on;
plot(mydata(:,1), mydata(:,2), 'LineWidth', 2);
line([l1x1 l1x2], [l1y1 l1y2], 'Color', 'red', 'LineStyle','--');
line([l2x1 l2x2], [l2y1 l2y2], 'Color', 'green', 'LineStyle','--');

% splx = [0:0.5:spotsrun.t_max];
% spl = spline(mydata(:,1),mydata(:,2),splx);
% plot(splx, spl, 'LineWidth', 2);

