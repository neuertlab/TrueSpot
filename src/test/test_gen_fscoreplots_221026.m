%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
% ========================== Image Channels ==========================

save_stem_rna = [ImgDir '\data\preprocess\feb2018\Tsix_AF594\Tsix\Tsix-AF594_IMG1_all_3d'];

% ========================== Generate Plots ==========================

path_spotsrun = [save_stem_rna '_rnaspotsrun.mat'];
if isfile(path_spotsrun)
    %Load run data
    spotsrun = RNASpotsRun.loadFrom(path_spotsrun);
    spotsrun.out_stem = save_stem_rna;
    spotsrun.saveMe();

    %Check threshold
    if isempty(spotsrun.threshold_results)
        spotsrun.threshold_results = RNAThreshold.runSavedParameters(rnaspots_run);
    end

    %Generate spot count graph
    [spotsrun, spots_table, ~] = spotsrun.loadZTrimmedTables_Sample();
    logspots = log10(spots_table);
    logspots(:,1) = spots_table(:,1);
    handle1 = RNAThreshold.plotThreshRanges(spotsrun, logspots, 'log10(Spot Count)', [], 505);
    %clear logspots;

    %Load refset if available and generate fscore plot
    fscores = RNA_Threshold_SpotSelector.loadFScores(save_stem_rna);
    if ~isempty(fscores)
        tcount = size(fscores,1);
        plot_scores = NaN(tcount,2);
        plot_scores(:,1) = spots_table(:,1);
        plot_scores(:,2) = fscores(:,1);
        clear fscores;
        handle2 = RNAThreshold.plotThreshRanges(spotsrun, plot_scores, 'F-Score', [0 1], 615);
        clear plot_scores;
    end

    clear spots_table;
end