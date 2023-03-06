%
%%  !! UPDATE TO YOUR BASE DIR
DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Settings ==========================
addpath('./core');

ImgName = 'mESC4d_Tsix-AF594';
RefMode = false;
NewAnno = true;
JustLoad = false;

% ========================== Read Table ==========================

InputTablePath = [DataDir filesep 'test_images.csv'];
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

outstem = [DataDir filesep 'data' filesep 'rsfish' getRSDBGroupOutputDir(ImgName) ImgName filesep 'RSFISH_' ImgName];
refstem_raw = getTableValue(image_table, r, 'OUTSTEM');
refstem = [DataDir replace(refstem_raw, '/', filesep)];

if NewAnno | ~RNA_Threshold_SpotSelector.selectorExists(outstem)
    addpath('./test/datadump');
    RefCompare_RSFish(refstem, outstem, false);
    selector = RNA_Threshold_SpotSelector.openSelector(outstem, false);
else
    selector = RNA_Threshold_SpotSelector.openSelector(outstem, false);
end

if JustLoad
    return;
end

selector.crosshair_color = [1.0 1.0 0.0];
if RefMode
    selector.launchRefSelectGUI();
else
    selector.launchGUI();
end

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

function outdir = getRSDBGroupOutputDir(imgname)
    outdir = [];
    if isempty(imgname); return; end

    if startsWith(imgname, 'mESC4d_')
        outdir = [filesep 'mESC4d' filesep];
    elseif startsWith(imgname, 'scrna_')
        outdir = [filesep 'scrna' filesep];
    elseif startsWith(imgname, 'mESC_loday_')
        outdir = [filesep 'mESC_loday' filesep];
    elseif startsWith(imgname, 'scprotein_')
        outdir = [filesep 'scprotein' filesep];
    elseif startsWith(imgname, 'sim_')
        outdir = [filesep 'sim' filesep];
    elseif startsWith(imgname, 'histonesc_')
        outdir = [filesep 'histonesc' filesep];
    elseif startsWith(imgname, 'ROI0')
        outdir = [filesep 'munsky_lab' filesep];
    elseif startsWith(imgname, 'sctc_')
        inparts = split(imgname, '_');
        ch = 'CH1';
        if endsWith(imgname, 'CTT1')
            ch = 'CH2';
        end
        outdir = [filesep 'yeast_tc' filesep inparts{2,1} filesep ch];
    elseif startsWith(imgname, 'rsfish_')
        outdir = [filesep 'rsfish' filesep];
    elseif startsWith(imgname, 'simvar_')
        outdir = [filesep 'simvar' filesep];
    end

end