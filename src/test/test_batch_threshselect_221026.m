%
%%  !! UPDATE TO YOUR BASE DIR
%ImgDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%ImgDir = 'C:\Users\hospelb\labdata\imgproc';
ImgDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
% ========================== Image Channels ==========================

RefStem = [ImgDir filesep ''];
SearchDir = [ImgDir filesep ''];

ProbeName = '';
TargetName = ''; %Have to match so don't do all channels in image.

% ========================== Load Reference FScores ==========================

fprintf("Loading reference fscore curve!\n");
fscores = RNA_Threshold_SpotSelector.loadFScores(RefStem);
if isempty(fscores)
    fprintf("Provided reference path has no fscore curve!\n");
    return;
end

[~, sugg_min, sugg_max] = RNAThreshold.threshSuggestionFromFScores(fscores);

% ========================== Scan search directory for spotruns to rethresh ==========================

fprintf("Recursively scanning search directory...\n");
scandir_rec(SearchDir);

function scandir_rec(dirpath, probe, target, fscores)
    dircontents = dir(dirpath);
    ccount = size(dircontents,1);
    for i = 1:ccount
        if dircontents(i).isdir
            scandir_rec([dirpath filesep dircontents(i).name], probe, target, fscores);
        else
            fname = dircontents(i).name;
            if endsWith(fname, '_rnaspotsrun.mat')
                fullpath = [dirpath filesep fname];
                rnaspots_run = RNASpotsRun.loadFrom(fullpath);
                if isempty(rnaspots_run.type_probe); continue; end
                if isempty(rnaspots_run.type_target); continue; end
                if ~strcmp(rnaspots_run.type_probe, ProbeName); continue; end
                if ~strcmp(rnaspots_run.type_target, TargetName); continue; end

                fprintf("Image run found: %s\n", rnaspots_run.img_name);
                param_info = RNAThreshold.paramsFromSpotsrun(rnaspots_run);
                param_info.suggestion = [sugg_min sugg_max];
                [rnaspots_run, param_info.sample_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Sample();
                [rnaspots_run, param_info.control_spot_table, ~] = rnaspots_run.loadZTrimmedTables_Control();
                rnaspots_run.threshold_results = RNA_Threshold_Common.estimateThreshold(param_info);
                rnaspots_run.intensity_threshold = rnaspots_run.threshold;
                rnaspots_run.saveMe();
            end
        end
    end
end
