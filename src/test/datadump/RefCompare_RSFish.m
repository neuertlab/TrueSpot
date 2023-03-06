%
%%

function RefCompare_RSFish(ref_data_path, rsfish_save_stem, ref_only)
%Generates a spotanno obj for an RSFish run

%Load coordtable
%rsfish_coord_path = [rsfish_save_stem '_coordTable.mat'];
%load(rsfish_coord_path, 'coord_table', 'writerver', 'writer_ver_str');
rsfish_spottbl_path = [rsfish_save_stem '_spotTable.mat'];
if ~isfile(rsfish_spottbl_path)
    fprintf('Run not found for %s. Returning...\n', rsfish_save_stem);
    return;
end

load(rsfish_spottbl_path, 'spot_table'); %For th x values

%Correct coords if RSFish2Mat is an older version
rsfish_coordtbl_path = [rsfish_save_stem '_coordTable.mat'];
mfinfo = who('-file', rsfish_coordtbl_path);
adj_coords = false;
if ~isempty(find(ismember(mfinfo, 'writerver'),1))
    load(rsfish_coordtbl_path, 'writerver');
    if writerver < 23022700
        adj_coords = true;
    end
else
    adj_coords = true;
end

writerver = 23022700;
writer_ver_str = 'v 23.02.27.0';
if adj_coords
    fit_tbl_path = [rsfish_save_stem '_fitTable.mat']; 
    
    load(rsfish_coordtbl_path, 'coord_table');
    tcount = size(coord_table,1);
    for t = 1:tcount
        ct = coord_table{t,1};
        ct = ct+1;
        coord_table{t,1} = ct;
    end
    save(rsfish_coordtbl_path, 'coord_table', 'writer_ver_str', 'writerver');
    clear coord_table;
    
    if isfile(fit_tbl_path)
        load(fit_tbl_path, 'fit_table');
        tcount = size(fit_table,1);
        for t = 1:tcount
            ct = fit_table{t,1};
            ct(:,1) = ct(:,1) + 1;
            ct(:,2) = ct(:,2) + 1;
            ct(:,3) = ct(:,3) + 1;
            fit_table{t,1} = ct;
        end
        save(fit_tbl_path, 'fit_table', 'writer_ver_str', 'writerver');
        clear fit_table;
    end
    
    clear ct;
end

%Load reference spotanno
if ref_only
    %Just a reference .mat file. Load into new spot anno.
    load(ref_data_path, 'ref_coord_tbl');
    myanno = RNA_Threshold_SpotSelector;
    myanno = myanno.initializeNew(rsfish_save_stem, 1, spot_table(:,1));
    myanno.ref_coords = ref_coord_tbl;
    clear ref_coord_tbl;
    clear spot_table;
else
    refanno = RNA_Threshold_SpotSelector.openSelector(ref_data_path,true);
    [~, myanno] = refanno.makeCopy();
    
    %Replace coords with deepblink callset
    load(rsfish_coordtbl_path, 'coord_table');
    myanno = myanno.loadNewSpotset(spot_table(:,1), coord_table);
    myanno.save_stem = rsfish_save_stem;
    
    clear spot_table;
    clear coord_table;
    clear refanno;
end
myanno.toggle_singleSlice = true;
myanno.toggle_allz = false;

%Try snap.
myanno = myanno.refSnapToAutoSpots(1);

myanno.f_scores_dirty = true;
myanno.ref_last_modified = datetime;

myanno = myanno.updateFTable();
myanno.saveMe();

clear myanno;

end