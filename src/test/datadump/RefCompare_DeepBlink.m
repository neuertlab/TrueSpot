%
%%

function RefCompare_DeepBlink(ref_data_path, deepblink_save_stem, ref_only)
%Generates a spotanno obj for a DeepBlink run

%Load coordtable and correct if needed.
writerver = 0;
adj_coords = false;
deep_coord_path = [deepblink_save_stem '_coordTable.mat'];
mfinfo = who('-file', deep_coord_path);
xyswap = false;
if ~isempty(find(ismember(mfinfo, 'writerver'),1))
    load(deep_coord_path, 'writerver');
    if writerver == 23022400
        xyswap = true;
    end
    if writerver < 23020900
        adj_coords = true;
    end
else
    adj_coords = true;
end
%clear mfinfo;

writerver = 23022700;
writer_ver_str = 'v 23.02.27.0';
ceilinged_bool = false;
if xyswap
    load(deep_coord_path, 'coord_table');
    
    tcount = size(coord_table,1);
    for t = 1:tcount
        ct = coord_table{t,1};
        if adj_coords
            ct = ct + 1;
        end
        tmp = ct(:,1);
        ct(:,1) = ct(:,2);
        ct(:,2) = tmp;
        coord_table{t,1} = ct;
    end
    
    save(deep_coord_path, 'coord_table', 'writer_ver_str', 'writerver','ceilinged_bool');
    clear coord_table;
    clear ct;
end

if adj_coords
    load(deep_coord_path, 'coord_table');
    tcount = size(coord_table,1);
    for t = 1:tcount
        ct = coord_table{t,1};
        ct = ct+1;
        coord_table{t,1} = ct;
    end
    save(deep_coord_path, 'coord_table', 'writer_ver_str', 'writerver','ceilinged_bool');
    clear coord_table;
    clear ct;
end

if writerver < 23022300
    prob_cutoffs = [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];
else
    spot_table_path = [deepblink_save_stem '_spotTable.mat'];
    load(spot_table_path, 'spot_table');
    prob_cutoffs = transpose(spot_table(:,1));
    clear spot_table;
end

%Load reference spotanno
if ref_only
    %Just a reference .mat file. Load into new spot anno.
    myanno = RNA_Threshold_SpotSelector;
    myanno = myanno.initializeNew(deepblink_save_stem, 1, transpose(prob_cutoffs));
    load(ref_data_path, 'ref_coord_tbl');
    myanno.ref_coords = ref_coord_tbl;
    clear ref_coord_tbl;
else
    refanno = RNA_Threshold_SpotSelector.openSelector(ref_data_path,true);
    [~, myanno] = refanno.makeCopy();
    
    needs_ceiling = true;
    if ~isempty(find(ismember(mfinfo, 'ceilinged_bool'),1))
        load(deep_coord_path, 'ceilinged_bool');
        needs_ceiling = ~ceilinged_bool;
    end
    
    %Replace coords with deepblink callset
    load(deep_coord_path, 'coord_table');
    if needs_ceiling
        ceilinged_bool = true;
        
        %Make sure coords do not exceed image dims
        idims = refanno.getImageDimensions();
        T = size(coord_table,1);
        for t = 1:T
            these_coords = coord_table{t,1};
            these_coords(:,1) = min(these_coords(:,1), idims.x);
            these_coords(:,2) = min(these_coords(:,2), idims.y);
            these_coords(:,3) = min(these_coords(:,3), idims.z);
            coord_table{t,1} = these_coords;
        end
        
        save(deep_coord_path, 'coord_table', 'writer_ver_str', 'writerver','ceilinged_bool');
    end
    
    myanno = myanno.loadNewSpotset(transpose(prob_cutoffs), coord_table);
    myanno.save_stem = deepblink_save_stem;
    
    clear coord_table;
    clear refanno;
end
myanno.toggle_singleSlice = true;
myanno.toggle_allz = false;

%CUSTOM SNAP (DeepBlink only works in 2D, so it looks at each slice!)
addpath('./test');
addpath('./test/datadump');
myanno = SpotAnnoSnap_DeepBlink(myanno);
%myanno = myanno.refSnapToAutoSpots(1);

myanno.f_scores_dirty = true;
myanno.ref_last_modified = datetime;

myanno = myanno.updateFTable();
myanno.saveMe();

clear myanno;

end