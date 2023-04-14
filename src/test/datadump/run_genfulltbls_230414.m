%
%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%BaseDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = [BaseDir '\data'];
OutputDir = [DataDir '\callsets'];

addpath('./core');
addpath('./test');

% ========================== Constants ==========================

START_INDEX = 1;
END_INDEX = 1;

DO_HOMEBREW = true;
DO_BIGFISH = true;
DO_RSFISH = false;
DO_DEEPBLINK = false;

% ========================== Load csv Table ==========================
%InputTablePath = [BaseDir filesep 'test_images_simytc.csv'];
%InputTablePath = [BaseDir filesep 'test_images_simvarmass.csv'];
InputTablePath = [BaseDir filesep 'test_images.csv'];

image_table = testutil_opentable(InputTablePath);

% ========================== Iterate through table entries ==========================
entry_count = size(image_table,1);

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

for r = START_INDEX:END_INDEX
    myname = getTableValue(image_table, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);
    
    hb_stem_base = getTableValue(image_table, r, 'OUTSTEM');
    hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
    
    callfiledir = [OutputDir filesep getSetOutputDirName(myname)];
    if ~isfolder(callfiledir)
        mkdir(callfiledir);
    end
    callfilepath = [callfiledir filesep myname '_allcalls.mat'];
    if isfile(callfilepath)
        load(callfilepath, 'allcalls_hb', 'allcalls_bf', 'allcalls_rs', 'allcalls_db', 'roc_hb', 'roc_bf', 'roc_rs', 'roc_db');
    else
        allcalls_hb = table.empty();
        allcalls_bf = table.empty();
        allcalls_rs = table.empty();
        allcalls_db = table.empty();
    end
    
    if DO_HOMEBREW
        if RNA_Threshold_SpotSelector.refsetExists(hb_stem)
            spotanno = RNA_Threshold_SpotSelector.openSelector(hb_stem);
            rawtbl = spotanno.genScoreResponseTable();
            allcalls_hb = tablify(rawtbl);
            
            
        end
    end
    
    
    %Save
    
end

% ========================== Helper functions ==========================

function respTable = tablify(raw_table)

varNames = {'isnap_x' 'isnap_y' 'isnap_z' 'dropout_thresh' 'is_true'};
varTypes = {'uint16' 'uint16' 'uint16' 'double' 'logical'};

spot_count = size(raw_table,1);
table_size = [spot_count size(varNames,2)];
respTable = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);

respTable(:,'isnap_x') = array2table(uint16(raw_table(:,1)));
respTable(:,'isnap_y') = array2table(uint16(raw_table(:,2)));
respTable(:,'isnap_z') = array2table(uint16(raw_table(:,3)));
respTable(:,'dropout_thresh') = array2table(double(raw_table(:,4)));
respTable(:,'is_true') = array2table(logical(raw_table(:,5)));

end

function dirname = getSetOutputDirName(imgname)
    inparts = split(imgname, '_');
    groupname = inparts{1,1};
    if strcmp(groupname, 'sctc')
        dirname = [groupname filesep inparts{2,1}];
    elseif strcmp(groupname, 'simvarmass')
        if contains(imgname, 'TMRL') | contains(imgname, 'CY5L')
            dirname = 'simytc';
        else
            dirname = groupname;
        end
    else
        dirname = groupname;
    end
end

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end
