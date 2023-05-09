%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');

% ========================== Constants ==========================

START_INDEX = 1;
END_INDEX = 1000;

NEW_TH_SETTING = 5;

% ========================== Load csv Table ==========================
%InputTablePath = [DataDir filesep 'test_images_simytc.csv'];
InputTablePath = [DataDir filesep 'test_images_simvarmass.csv'];
%InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

% ========================== Go through entries ==========================

entry_count = size(image_table,1);

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

for r = START_INDEX:END_INDEX
    myname = getTableValue(image_table, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

    hb_stem_base = getTableValue(image_table, r, 'OUTSTEM');
    hb_stem = [DataDir replace(hb_stem_base, '/', filesep)];

    spotsrun = RNASpotsRun.loadFrom(hb_stem);
    if ~isempty(spotsrun)
        spotsrun.out_stem = hb_stem;
        spotsrun = RNAThreshold.applyPreset(spotsrun, NEW_TH_SETTING);
        threshold_results = RNAThreshold.runSavedParameters(spotsrun);
        if ~isempty(threshold_results)
            spotsrun.threshold_results = threshold_results;
            spotsrun.intensity_threshold = threshold_results.threshold;
        end
        spotsrun.saveMe();

    else
        fprintf('Run not found. Skipping...\n');
    end

end

% ========================== Helper functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end