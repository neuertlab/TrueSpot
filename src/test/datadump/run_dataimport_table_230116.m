%
%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%BaseDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
addpath('./test');

% ========================== Constants ==========================

START_INDEX = 1;
END_INDEX = 38;

DO_HOMEBREW = true;
DO_BIGFISH = false;
DO_BIGFISHNR = false;
DO_RSFISH = false;
DO_DEEPBLINK = false;

HB_FIXED_TH = 0; %Maybe have a separate script for this?

% ========================== Other paths ==========================

DataFilePath = [BaseDir filesep 'data' filesep 'DataComposite.mat'];

% ========================== Load csv Table ==========================

AllFigDir = [ImgProcDir filesep 'figures' filesep 'curves'];

InputTablePath = [BaseDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Iterate through table entries ==========================
entry_count = size(imgtbl,1);

if ~isfile(DataFilePath)
    image_analyses(entry_count) = ImageResults;
    for r = 1:entry_count
        %Initialize
        image_analyses(r) = ImageResults.initializeNew();
        image_analyses(r).image_name = getTableValue(imgtbl, r, 'IMGNAME');
    end
else
    load(DataFilePath, 'image_analyses');
    %If entry_count is larger, add new entries to end.
    save_count = size(image_analyses,2);
    if entry_count > save_count
        image_analyses(entry_count) = ImageResults; %Expand.
        for r = save_count+1:entry_count
            %Initialize
            image_analyses(r) = ImageResults.initializeNew();
            image_analyses(r).image_name = getTableValue(imgtbl, r, 'IMGNAME');
        end
    end
end

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

for r = START_INDEX:END_INDEX
    myname = getTableValue(imgtbl, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

    %----------------- HB
    if DO_HOMEBREW
        hb_stem_base = getTableValue(imgtbl, r, 'OUTSTEM');
        hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
        [image_analyses(r), bool_okay] = image_analyses(r).importHBResults(hb_stem);
        if ~bool_okay
            fprintf('WARNING: HB Import of %s failed!\n', myname);
        end
    end

    %----------------- BF
    if DO_BIGFISH
        bf_stem_base = replace(getTableValue(imgtbl, r, 'BIGFISH_OUTSTEM'), '/bigfish/', '/bigfish/_rescaled/');
        bf_stem = [BaseDir replace(bf_stem_base, '/', filesep)];
        [image_analyses(r), bool_okay] = image_analyses(r).importBFResults(bf_stem, true);
        if ~bool_okay
            fprintf('WARNING: BF (Rescaled) Import of %s failed!\n', myname);
        end
    end

    %----------------- BFNR
    if DO_BIGFISHNR
        bf_stem_base = getTableValue(imgtbl, r, 'BIGFISH_OUTSTEM');
        bf_stem = [BaseDir replace(bf_stem_base, '/', filesep)];
        [image_analyses(r), bool_okay] = image_analyses(r).importBFResults(bf_stem, false);
        if ~bool_okay
            fprintf('WARNING: BF (Non rescaled) Import of %s failed!\n', myname);
        end
    end

    %----------------- RSFISH
    rsdb_dir_ext = getRSDBGroupOutputDir(myname);
    if DO_RSFISH
        if ~isempty(rsdb_dir_ext)
            rs_stem = [BaseDir filesep 'data' filesep 'rsfish' rsdb_dir_ext myname filesep 'RSFISH_' myname];
            %Try importing if no coord table found
            %Also delete any csvs still lying around
            [rs_dir, ~, ~] = fileparts(rs_stem);
            coord_table_path = [rs_stem '_coordTable.mat'];
            if ~isfile(coord_table_path)
                Main_RSFish2Mat(rs_dir, rs_stem);
            end
            if isfile(coord_table_path)
                deleteCSVs(rs_dir);
            end

            %Then check for spotanno.
            %If no spotanno, check for truthset.
            %Import truthset
            if ~RNA_Threshold_SpotSelector.refsetExists(rs_stem)
                hb_stem_base = getTableValue(imgtbl, r, 'OUTSTEM');
                hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
                if RNA_Threshold_SpotSelector.refsetExists(hb_stem)
                    RefCompare_RSFish(hb_stem, rs_stem, false);
                end
            end

            [image_analyses(r), bool_okay] = image_analyses(r).importRSResults(rs_stem);
            if ~bool_okay
                fprintf('WARNING: RS-FISH Import of %s failed!\n', myname);
            end
        end
    end

    %----------------- DeepBlink
    if DO_DEEPBLINK
        if ~isempty(rsdb_dir_ext)
            db_stem = [BaseDir filesep 'data' filesep 'deepblink' rsdb_dir_ext myname filesep 'DeepBlink_' myname];

            %check for spotanno.
            %If no spotanno, check for truthset.
            %Import truthset
            if ~RNA_Threshold_SpotSelector.refsetExists(db_stem)
                hb_stem_base = getTableValue(imgtbl, r, 'OUTSTEM');
                hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
                if RNA_Threshold_SpotSelector.refsetExists(hb_stem)
                    RefCompare_DeepBlink(hb_stem, db_stem, false);
                end
            end

            [image_analyses(r), bool_okay] = image_analyses(r).importDBResults(db_stem);
            if ~bool_okay
                fprintf('WARNING: DeepBlink Import of %s failed!\n', myname);
            end
        end
    end

end

save(DataFilePath, 'image_analyses');

% ========================== Helper functions ==========================

function outdir = getRSDBGroupOutputDir(imgname)
    outdir = [];
    if isempty(imgname); return; end

    if startsWith(imgname, 'mESC4d_')
        outdir = [filesep 'mESC4d' filesep];
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

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

function deleteCSVs(dir_path)
    dir_contents = dir(dir_path);
    fcount = size(dir_contents,1);
    for i = 1:fcount
        fentry = dir_contents(i,1);
        if fentry.isdir; continue; end
        if endsWith(fentry.name, '.csv')
            filepath = [dir_path filesep fentry.name];
            fprintf('CSV Found: %s - Deleting...\n', filepath);
            delete(filepath);
        end
    end
end