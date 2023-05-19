%
%%  !! UPDATE TO YOUR BASE DIR
%DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Settings ==========================
addpath('./core');

ImgName = 'simvarmass_a41fe574';
RefMode = false;
NewAnno = false;
JustLoad = false;

GuessMask = false;
NewMaskZMin = 10;
NewMaskZMax = 40;

% ========================== Read Table ==========================

%InputTablePath = [DataDir filesep 'test_images_simytc.csv'];
InputTablePath = [DataDir filesep 'test_images_simvarmass.csv'];
%InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

rec_row = 0;
rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    if strcmp(iname, ImgName)
        rec_row = r;
        break;
    end
end

if rec_row < 1
    fprintf('Image with name %s could not be found!\n', ImgName);
    return;
end

% ========================== Load Spots Anno ==========================

outstem = [DataDir replace(getTableValue(image_table, rec_row, 'OUTSTEM'), '/', filesep)];

spotsrun = RNASpotsRun.loadFrom(outstem);
if isempty(spotsrun)
    fprintf('Run with prefix "%s" was not found!\n', outstem);
    return;
end

%Update spotsrun paths
spotsrun.out_stem = outstem;
spotsrun.saveMe();

if NewAnno | ~RNA_Threshold_SpotSelector.selectorExists(outstem)
    selector = RNA_Threshold_SpotSelector;
    th_idx = round((spotsrun.t_max - spotsrun.t_min) ./ 2);
    selector = selector.initializeNew(outstem, th_idx, []);
    selector.z_min = 1;
    selector.z_max = spotsrun.idims_sample.z;
else
    selector = RNA_Threshold_SpotSelector.openSelector(outstem, true);
end

if GuessMask
    if NewMaskZMin > 0
        selector.z_min = NewMaskZMin;
        if selector.current_slice < NewMaskZMin
            selector.current_slice = NewMaskZMin;
        end
    end
    if NewMaskZMax > 0
        selector.z_max = NewMaskZMax;
        if selector.current_slice > NewMaskZMax
            selector.current_slice = NewMaskZMax;
        end
    end
  

    if ~isempty(selector.ref_coords)
        x_good = find(selector.ref_coords(:,1));
        x_min = min(selector.ref_coords(x_good,1), [], 'all', 'omitnan');
        x_max = max(selector.ref_coords(x_good,1), [], 'all', 'omitnan');
        y_good = find(selector.ref_coords(:,2));
        y_min = min(selector.ref_coords(y_good,2), [], 'all', 'omitnan');
        y_max = max(selector.ref_coords(y_good,2), [], 'all', 'omitnan');
        selector.selmcoords = zeros(4,1);
        selector.selmcoords(1,1) = x_min;
        selector.selmcoords(2,1) = x_max;
        selector.selmcoords(3,1) = y_min;
        selector.selmcoords(4,1) = y_max;
    end
end

if JustLoad
    return;
end

selector.crosshair_color = [1.0 1.0 0.0];
if RefMode
    selector.launchRefSelectGUI();
else
    if spotsrun.intensity_threshold > 0
        th_idx = RNAUtils.findThresholdIndex(spotsrun.intensity_threshold, transpose(selector.threshold_table(:,1)));
        selector.threshold_idx = th_idx;
    end
    
    selector.launchGUI();
end

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end
