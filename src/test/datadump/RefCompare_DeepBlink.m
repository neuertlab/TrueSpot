%
%%

function RefCompare_DeepBlink(ref_data_path, deepblink_save_stem, ref_only)
%Generates a spotanno obj for a DeepBlink run

%Load coordtable and correct if needed.
adj_coords = false;
deep_coord_path = [deepblink_save_stem '_coordTable.mat'];
mfinfo = who('-file', deep_coord_path);
if ~isempty(find(ismember(mfinfo, 'writerver'),1))
    load(deep_coord_path, 'writerver');
    if writerver < 23020900
        adj_coords = true;
    end
else
    adj_coords = true;
end
clear mfinfo;

if adj_coords
    writerver = 23020900;
    writer_ver_str = 'v 23.02.09.0';
    load(deep_coord_path, 'coord_table');
    tcount = size(coord_table,1);
    for t = 1:tcount
        ct = coord_table{t,1};
        ct = ct + 1;
        coord_table{t,1} = ct;
    end
    save(deep_coord_path, 'coord_table', 'writerver', 'writer_ver_str');
    clear coord_table;
end


prob_cutoffs = [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];
myanno = RNA_Threshold_SpotSelector;
myanno = myanno.initializeNew(deepblink_save_stem, 1, transpose(prob_cutoffs));
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

%CUSTOM SNAP (DeepBlink only works in 2D, so it looks at each slice!)
addpath('./test');
addpath('./test/datadump');
myanno = SpotAnnoSnap_DeepBlink(myanno);

myanno.f_scores_dirty = true;
myanno.ref_last_modified = datetime;

myanno = myanno.updateFTable();
myanno.saveMe();

clear myanno;

end