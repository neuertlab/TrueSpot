%
%%  !! UPDATE TO YOUR BASE DIR
BaseDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%BaseDir = 'D:\usr\bghos\labdat\imgproc';

ImgProcDir = 'D:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgProcDir = 'D:\usr\bghos\labdat\imgproc';

addpath('./core');
addpath('./test');

% ========================== Constants ==========================

START_INDEX = 159;
END_INDEX = 190;

DO_HOMEBREW = true;
DO_BIGFISH = true;
DO_BIGFISHNR = true;
DO_RSFISH = true;
DO_DEEPBLINK = true;

OVERWRITE_SPOTANNO_RS = true;
OVERWRITE_SPOTANNO_DB = true;
FORCE_HB_BF_RESNAP = true; %TODO

IMPORT_TS = true;
TS_IDX = 1;

HB_FIXED_TH = 0; %Maybe have a separate script for this?
RS_TH_IVAL = 0.1/250;

% ========================== Other paths ==========================

DataFilePath = [BaseDir filesep 'data' filesep 'DataComposite.mat'];

% ========================== Load csv Table ==========================

AllFigDir = [ImgProcDir filesep 'figures' filesep 'curves'];

InputTablePath = [BaseDir filesep 'test_images.csv'];
imgtbl = testutil_opentable(InputTablePath);

% ========================== Iterate through table entries ==========================
entry_count = size(imgtbl,1);

if ~isfile(DataFilePath)
    image_analyses(entry_count) = struct('imgname', '', 'analysis', []);
    %image_analyses(entry_count) = ImageResults;
    for r = 1:entry_count
        %Initialize
        image_analyses(r).imgname = getTableValue(imgtbl, r, 'IMGNAME');
        image_analyses(r).analysis = ImageResults.initializeNew();
        image_analyses(r).analysis.image_name = image_analyses(r).imgname;
    end
else
    load(DataFilePath, 'image_analyses');
    %If entry_count is larger, add new entries to end.
    save_count = size(image_analyses,2);
    if entry_count > save_count
        %image_analyses(entry_count) = ImageResults; %Expand.
        image_analyses(entry_count) = struct('imgname', '', 'analysis', []);
        for r = save_count+1:entry_count
            %Initialize
            image_analyses(r).imgname = getTableValue(imgtbl, r, 'IMGNAME');
            image_analyses(r).analysis = ImageResults.initializeNew();
            image_analyses(r).analysis.image_name = image_analyses(r).imgname;
        end
    end
end

if START_INDEX < 1; START_INDEX = 1; end
if END_INDEX > entry_count; END_INDEX = entry_count; end

for r = START_INDEX:END_INDEX
    myname = getTableValue(imgtbl, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);
    
    resetSim(imgtbl, r, BaseDir, ImgDir);
    
    %----------------- Truthset
    if IMPORT_TS
        %Check for truthset.
        hb_stem_base = getTableValue(imgtbl, r, 'OUTSTEM');
        hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
        if RNA_Threshold_SpotSelector.refsetExists(hb_stem)
            %Has a substantial refset. Import.
            image_analyses(r).analysis = image_analyses(r).analysis.setTruthset(hb_stem, TS_IDX);
        end
    end

    %----------------- BF
    if DO_BIGFISH
        bf_stem_base = replace(getTableValue(imgtbl, r, 'BIGFISH_OUTSTEM'), '/bigfish/', '/bigfish/_rescaled/');
        bf_stem = [BaseDir replace(bf_stem_base, '/', filesep)];
        if FORCE_HB_BF_RESNAP & RNA_Threshold_SpotSelector.selectorExists(bf_stem)
            RNA_Threshold_SpotSelector.fixPosTable(bf_stem);
            spotanno = RNA_Threshold_SpotSelector.openSelector(bf_stem, true);
            spotanno.toggle_singleSlice = true;
            spotanno.toggle_allz = false;
            
            %See if need to adjust the stop point...
            stopat = 20;
            [bfdir, ~, ~] = fileparts(bf_stem);
            [~, ~, bfthresh] = BigfishCompare.readSummaryTxt([bfdir filesep 'summary.txt']);
            bf_st_path = [bf_stem 'spotTable.mat'];
            if isfile(bf_st_path)
                load(bf_st_path, 'spot_table');
                th_idx = RNAUtils.findThresholdIndex(bfthresh, transpose(spot_table(:,1)));
                clear spot_table;
            else
                %Assume 1 -> 10
                th_idx = bfthresh - 9;
            end
            if th_idx <= 20
                stopat = th_idx - 1;
                if stopat < 1; stopat = 1; end
            end
            
            spotanno = spotanno.refSnapToAutoSpots(stopat);
            spotanno.saveMe();
        end
        [image_analyses(r).analysis, bool_okay] = image_analyses(r).analysis.importBFResults(bf_stem, true, TS_IDX);
        if ~bool_okay
            fprintf('WARNING: BF (Rescaled) Import of %s failed!\n', myname);
        end
    end

    %----------------- BFNR
    if DO_BIGFISHNR
        bf_stem_base = getTableValue(imgtbl, r, 'BIGFISH_OUTSTEM');
        bf_stem = [BaseDir replace(bf_stem_base, '/', filesep)];
        if FORCE_HB_BF_RESNAP & RNA_Threshold_SpotSelector.selectorExists(bf_stem)
            RNA_Threshold_SpotSelector.fixPosTable(bf_stem);
            spotanno = RNA_Threshold_SpotSelector.openSelector(bf_stem, true);
            spotanno.toggle_singleSlice = true;
            spotanno.toggle_allz = false;
            
            stopat = 20;
            [bfdir, ~, ~] = fileparts(bf_stem);
            [~, ~, bfthresh] = BigfishCompare.readSummaryTxt([bfdir filesep 'summary.txt']);
            bf_st_path = [bf_stem 'spotTable.mat'];
            if isfile(bf_st_path)
                load(bf_st_path, 'spot_table');
                th_idx = RNAUtils.findThresholdIndex(bfthresh, transpose(spot_table(:,1)));
                clear spot_table;
            else
                %Assume 1 -> 10
                th_idx = bfthresh - 9;
            end
            if th_idx <= 20
                stopat = th_idx - 1;
                if stopat < 1; stopat = 1; end
            end
            
            spotanno = spotanno.refSnapToAutoSpots(stopat);
            spotanno.saveMe();
        end
        [image_analyses(r).analysis, bool_okay] = image_analyses(r).analysis.importBFResults(bf_stem, false, TS_IDX);
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
            
            if ~isfile(coord_table_path)
                fprintf('WARNING: RS-FISH run for %s not found! Skipping...\n', myname);
                continue;
            end
            
            %Then check for spotanno.
            %If no spotanno, check for truthset.
            %Import truthset
            if OVERWRITE_SPOTANNO_RS | ~RNA_Threshold_SpotSelector.refsetExists(rs_stem)
                hb_stem_base = getTableValue(imgtbl, r, 'OUTSTEM');
                hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
                if RNA_Threshold_SpotSelector.refsetExists(hb_stem)
                    RefCompare_RSFish(hb_stem, rs_stem, false);
                end
            end

            [image_analyses(r).analysis, bool_okay] = image_analyses(r).analysis.importRSResults(rs_stem, RS_TH_IVAL, TS_IDX);
            if ~bool_okay
                fprintf('WARNING: RS-FISH Import of %s failed!\n', myname);
            end
        end
    end

    %----------------- DeepBlink
    if DO_DEEPBLINK
        if ~isempty(rsdb_dir_ext)
            db_stem = [BaseDir filesep 'data' filesep 'deepblink' rsdb_dir_ext myname filesep 'DeepBlink_' myname];

            [db_dir, ~, ~] = fileparts(db_stem);
            coord_table_path = [db_stem '_coordTable.mat'];
            %if ~isfile(coord_table_path)
            try
                Main_DeepBlink2Mat([db_stem '.csv'], db_stem);
            catch
                fprintf('WARNING: DeepBlink run for %s could not be imported! Skipping...\n', myname);
                continue;
            end
            %end
            
            if ~isfile(coord_table_path)
                fprintf('WARNING: DeepBlink run for %s not found! Skipping...\n', myname);
                continue;
            end
            
            %check for spotanno.
            %If no spotanno, check for truthset.
            %Import truthset
            if OVERWRITE_SPOTANNO_DB | ~RNA_Threshold_SpotSelector.refsetExists(db_stem)
                hb_stem_base = getTableValue(imgtbl, r, 'OUTSTEM');
                hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
                if RNA_Threshold_SpotSelector.refsetExists(hb_stem)
                    RefCompare_DeepBlink(hb_stem, db_stem, false);
                end
            end

            [image_analyses(r).analysis, bool_okay] = image_analyses(r).analysis.importDBResults(db_stem, TS_IDX);
            if ~bool_okay
                fprintf('WARNING: DeepBlink Import of %s failed!\n', myname);
            end
        end
    end
    
    %----------------- HB
    if DO_HOMEBREW
        hb_stem_base = getTableValue(imgtbl, r, 'OUTSTEM');
        hb_stem = [BaseDir replace(hb_stem_base, '/', filesep)];
        
        %Resnap truthset if sim.
        if FORCE_HB_BF_RESNAP | startsWith(myname, 'sim_') | startsWith(myname, 'simvar_') | startsWith(myname, 'rsfish_sim_')
            if RNA_Threshold_SpotSelector.refsetExists(hb_stem)
                spotanno = RNA_Threshold_SpotSelector.openSelector(hb_stem, true);
                spotanno.toggle_singleSlice = true;
                spotanno.toggle_allz = false;
                
                %See if need to adjust the stop point...
                stopat = 20;
                spotsrun = RNASpotsRun.loadFrom(hb_stem);
                if ~isempty(spotsrun)
                    if spotsrun.intensity_threshold > 0
                        th_idx = RNAUtils.findThresholdIndex(spotsrun.intensity_threshold, transpose(spotanno.threshold_table(:,1)));
                        if th_idx <= 20
                            stopat = th_idx - 1;
                            if stopat < 1; stopat = 1; end
                        end
                    end
                    clear spotsrun;
                end
                
                spotanno = spotanno.refSnapToAutoSpots(stopat);
                spotanno.saveMe();
            end
        end
        
        [image_analyses(r).analysis, bool_okay] = image_analyses(r).analysis.importHBResults(hb_stem, HB_FIXED_TH, TS_IDX);
        if ~bool_okay
            fprintf('WARNING: HB Import of %s failed!\n', myname);
        end
    end

    %Maybe save after each successful image import?
    save(DataFilePath, 'image_analyses');
end

save(DataFilePath, 'image_analyses');

% ========================== Helper functions ==========================

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

function key_mtx = keyStructs2Mtx(my_key)
    ptcount = size(my_key,2);
    key_mtx = uint16(zeros(ptcount,3));
    key_mtx(:,1) = [my_key.x];
    key_mtx(:,2) = [my_key.y];
    key_mtx(:,3) = [my_key.z];
end

function resetSimSF(image_table, row_index, BaseDir, ImgDir)
    %For simfish sims
    myname = getTableValue(image_table, row_index, 'IMGNAME');
    srcpath_raw = getTableValue(image_table, row_index, 'IMAGEPATH');
    outstem_raw = getTableValue(image_table, row_index, 'OUTSTEM');
    
    srcpath = [ImgDir replace(srcpath_raw, '/', filesep)];
    hb_path = [BaseDir replace(outstem_raw, '/', filesep)];
    
    key = [];
    if endsWith(srcpath, '.mat')
        %Use this one directly.
        load(srcpath, 'key');
    elseif endsWith(srcpath, '.tif')
        %Find mat file and load from that.
        srcpath = replace(srcpath, [filesep 'tif' filesep], filesep);
        srcpath = replace(srcpath, '.tif', '.mat');
        load(srcpath, 'key');
    end
    
    if isempty(key)
        fprintf('ERROR: Could not load sim truthset for %s!\n', myname);
        return; 
    end
    
    %Look for spotsanno in hb dir
    %If not there, create new one.
    if ~RNA_Threshold_SpotSelector.selectorExists(hb_path)
        spotsrun = RNASpotsRun.loadFrom(hb_path);
        selector = RNA_Threshold_SpotSelector;
        th_idx = round((spotsrun.t_max - spotsrun.t_min) ./ 2);
        selector = selector.initializeNew(hb_path, th_idx, []);
        selector.z_min = 1;
        selector.z_max = spotsrun.idims_sample.z;
    else
        selector = RNA_Threshold_SpotSelector.openSelector(hb_path, true);
    end
    
    %Swap out refset and save
    selector.ref_coords = keyStructs2Mtx(key); %Importer already adjusts to 1 based coords.
    selector.ref_last_modified = datetime;
    selector.f_scores_dirty = true;
    selector.saveMe();
end

function resetSimRS(image_table, row_index, BaseDir, ImgDir)
    %For rsfish sims
    myname = getTableValue(image_table, row_index, 'IMGNAME');
    srcpath_raw = getTableValue(image_table, row_index, 'IMAGEPATH');
    outstem_raw = getTableValue(image_table, row_index, 'OUTSTEM');
    
    %I had to convert the locs to csv because MATLAB is a fussbudget
    srcpath = [ImgDir replace(srcpath_raw, '/', filesep)];
    srcpath = replace(srcpath, '.tif', '.csv');
    hb_path = [BaseDir replace(outstem_raw, '/', filesep)];
    
    if ~isfile(srcpath)
        fprintf('ERROR: Could not find sim truthset for %s!\n', myname);
        return;
    end
    
    import_table = table2array(readtable(srcpath,'ReadVariableNames',false));
    import_table = uint16(round(import_table)) + 1; %Make it less spicy
    temp = import_table(:,2);
    import_table(:,2) = import_table(:,1);
    import_table(:,1) = temp;
    
    %Look for spotsanno in hb dir
    %If not there, create new one.
    if ~RNA_Threshold_SpotSelector.selectorExists(hb_path)
        spotsrun = RNASpotsRun.loadFrom(hb_path);
        selector = RNA_Threshold_SpotSelector;
        th_idx = round((spotsrun.t_max - spotsrun.t_min) ./ 2);
        selector = selector.initializeNew(hb_path, th_idx, []);
        selector.z_min = 1;
        selector.z_max = spotsrun.idims_sample.z;
    else
        selector = RNA_Threshold_SpotSelector.openSelector(hb_path, true);
    end
    
    %Swap out refset and save
    selector.ref_coords = import_table(:,1:3);
    selector.ref_last_modified = datetime;
    selector.f_scores_dirty = true;
    selector.saveMe();
end

function resetSim(image_table, row_index, BaseDir, ImgDir)
    %Checks if sim image, and resets truthset if so
    myname = getTableValue(image_table, row_index, 'IMGNAME');
    if startsWith(myname, 'sim_')
        resetSimSF(image_table, row_index, BaseDir, ImgDir);
    elseif startsWith(myname, 'simvar_')
        resetSimSF(image_table, row_index, BaseDir, ImgDir);
    elseif startsWith(myname, 'rsfish_sim')
        resetSimRS(image_table, row_index, BaseDir, ImgDir);
    end
end