%
%%

BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc'; %workstation

RefStem = [BaseDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];

addpath('./core');

spotsrun = RNASpotsRun.loadFrom(RefStem);
threshold_results = RNAThreshold.runDefaultParametersSpotsRun(spotsrun,1);

winsz_count = size(threshold_results.window_sizes,2);
data_points = size(threshold_results.x,1);
plotdata = NaN(data_points-1,1);
plotdata(:,1) = threshold_results.x(1:data_points-1);
for i = 1:winsz_count
    res = threshold_results.test_winsc(1,i);
    plotdata(:,2) = log10(threshold_results.window_scores(:,i));
    figh = Seglr2.renderFit(plotdata, res.spline_fit, true, i);
end

