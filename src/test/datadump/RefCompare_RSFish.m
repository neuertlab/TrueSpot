%
%%

function RefCompare_RSFish(ref_data_path, rsfish_save_stem, ref_only)
%Generates a spotanno obj for an RSFish run

%Load coordtable
%rsfish_coord_path = [rsfish_save_stem '_coordTable.mat'];
%load(rsfish_coord_path, 'coord_table', 'writerver', 'writer_ver_str');
rsfish_spottbl_path = [rsfish_save_stem '_spotTable.mat'];
load(rsfish_spottbl_path, 'spot_table'); %For th x values

myanno = RNA_Threshold_SpotSelector;
myanno = myanno.initializeNew(rsfish_save_stem, 1, spot_table(:,1));
clear spot_table;

%Load reference spotanno
if ref_only
    %Just a reference .mat file. Load into new spot anno.
    load(ref_data_path, 'ref_coord_tbl');
    myanno.ref_coords = ref_coord_tbl;
    clear ref_coord_tbl;
else
    refanno = RNA_Threshold_SpotSelector.openSelector(ref_data_path,true);
    
    %Copy mask and ref coords
    myanno.selmcoords = refanno.selmcoords;
    myanno.z_min = refanno.z_min;
    myanno.z_max = refanno.z_max;
    myanno.ref_coords = refanno.ref_coords;
    
    clear refanno;
end

%Try snap.
myanno = myanno.refSnapToAutoSpots(1);

myanno.f_scores_dirty = true;
myanno.ref_last_modified = datetime;

myanno = myanno.updateFTable();
myanno.saveMe();

clear myanno;

end