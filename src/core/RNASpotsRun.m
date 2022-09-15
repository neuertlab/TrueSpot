%
%%

classdef RNASpotsRun
    
    properties
        img_name;
        tif_path;
        out_dir;
        cellseg_path;
        ctrl_path;
        
        bkg_path;
        bkg_filter_stem;
        out_stem;
        ctrl_stem;
        
        rna_ch;
        light_ch;
        total_ch;
        ctrl_ch;
        ctrl_chcount;
        
        idims_sample;
        idims_ctrl;
        
        type_probe; %ie. 'TMR', 'AF594' 
        type_species; %ie. Mus musculus, Saccharomyces cerevisiae
        type_cell; %ie. 'S. cerevisiae', 'Yeast', 'mESC'
        type_target; %ie. 'Tsix', 'Opy2', 'H3K4me2'
        type_targetmol; %ie. 'lncRNA', 'Histone Mark', 'Protein'
        
        t_min;
        t_max;
        ztrim;
        z_min; %User set boundary
        z_max; %User set boundary
        z_min_apply; %What is actually used
        z_max_apply; %What is actually used
        ztrim_auto;
        ttune_winsize; %DEPRECATED
        %ttune_wscorethresh;
        ttune_madfactor; %DEPRECATED
        intensity_threshold;
        
        threshold_results;
        ttune_madf_min;
        ttune_madf_max;
        ttune_winsz_min;
        ttune_winsz_max;
        ttune_winsz_incr;
        ttune_spline_itr;
        ttune_use_rawcurve;
        ttune_use_diffcurve;
        ttune_fit_to_log;
        ttune_fit_strat; %Enum: 0,1,2,3
        ttune_reweight_fit;
        ttune_thweight_med;
        ttune_thweight_fit;
        ttune_thweight_fisect;
        
        overwrite_output;
    end
    
    methods
        
        function obj = saveMe(obj)
            outpath = [obj.out_stem '_rnaspotsrun.mat'];
            rnaspots_run = obj;
            save(outpath, 'rnaspots_run');
        end
        
        function [obj, thresh_table] = loadThresholdTable(obj)
            thresh_table = [];
            [~, spot_table] = obj.loadSpotsTable();
            if ~isempty(spot_table)
                thresh_table = transpose(spot_table(:,1));
            end
        end
        
        function [obj, spot_table] = loadSpotsTable(obj)
            spot_table = [];
            tbl_path = [obj.out_stem '_spotTable.mat'];
            if isfile(tbl_path)
                load(tbl_path, 'spot_table');
            end
        end
        
        function [obj, spot_table] = loadControlSpotsTable(obj)
            spot_table = [];
            tbl_path = [obj.ctrl_stem '_spotTable.mat'];
            if isfile(tbl_path)
                load(tbl_path, 'spot_table');
            end
        end
        
        function [obj, coord_table] = loadCoordinateTable(obj)
            coord_table = [];
            tbl_path_RNA = [obj.out_stem '_coordTable.mat'];
            if isfile(tbl_path_RNA)
                load(tbl_path_RNA, 'coord_table');
            end
        end
        
        function [obj, coord_table] = loadControlCoordinateTable(obj)
            coord_table = [];
            tbl_path = [obj.ctrl_stem '_coordTable.mat'];
            if isfile(tbl_path)
                load(tbl_path, 'coord_table');
            end
        end
          
        function [obj, background_mask] = loadBackgroundMask(obj)
            background_mask = [];
            if isfile([obj.bkg_path '.mat'])
                load(obj.bkg_path, 'background_mask');
            end
        end
        
        function [obj, my_images] = loadImageViewStructs(obj)
            my_images = [];
            mypath = [obj.out_stem '_imgviewstructs.mat'];
            if isfile(mypath)
                load(mypath, 'my_images');
            end
        end
        
        function [obj, img_filter] = loadFilteredImage(obj)
            img_filter = [];
            mypath = [obj.out_stem '_prefilteredIMG.mat'];
            if isfile(mypath)
                load(mypath, 'img_filter');
            end
        end
        
        function obj = updateBackgroundFilteredCoords(obj)
            [~, background_mask] = obj.loadBackgroundMask();
            [~, coord_table] = obj.loadCoordinateTable();
            [~, th_list] = obj.loadThresholdTable();
            
            if isempty(background_mask)
                return;
            end
            
            [masked_spot_table, masked_coord_table] = RNA_Threshold_Common.mask_spots(background_mask, coord_table, th_list);
            coord_table = masked_coord_table;
            spot_table = masked_spot_table;
            save([obj.bkg_filter_stem '_coordTable'], 'coord_table');
            save([obj.bkg_filter_stem '_spotTable'], 'spot_table');
            
            %Image structs
            [~, my_images] = obj.loadImageViewStructs();
            my_images(1).image = immultiply(my_images(1).image, background_mask);
            my_images(2).image = immultiply(my_images(2).image, background_mask);
            save([obj.bkg_filter_stem '_imgviewstructs'], 'my_images');  
        end
        
        function obj = updateImageDimensions(obj)
            if isempty(obj.idims_sample)
                runinfo_path = [obj.out_stem '_runinfo'];
                if isfile(runinfo_path)
                    %Load the run info.
                    load(runinfo_path, 'idims');
                else
                    %Reload TIF and count manually.
                    [X,Y,Z] = GetTifDims(obj.tif_path, obj.total_ch);
                    idims = struct('x', X, 'y', Y, 'z', Z); 
                end
                idims
                obj.idims_sample = idims;
            end
            
            if isempty(obj.idims_ctrl)
                if ~isempty(obj.ctrl_path)
                    runinfo_path = [obj.ctrl_stem '_runinfo'];
                    if isfile(runinfo_path)
                        %Load the run info.
                        load(runinfo_path, 'idims');
                    else
                        %Reload TIF and count manually.
                        [X,Y,Z] = GetTifDims(obj.ctrl_path, obj.ctrl_chcount);
                        idims = struct('x', X, 'y', Y, 'z', Z);
                    end
                    obj.idims_ctrl = idims;
                else
                    %Assumed to be background
                    obj.idims_sample = idims;
                end
            end
        end
        
        function [obj, tif] = loadSampleTif(obj, verbosity)
            if (nargin < 2); verbosity = 2; end
            [tif, obj.idims_sample] = LoadTif(obj.tif_path, obj.total_ch, [obj.rna_ch, obj.light_ch], verbosity);
        end
        
        function [obj, ch_dat] = loadControlChannel(obj, verbosity)
            ch_dat = [];
            if (nargin < 2); verbosity = 2; end
            if ~isempty(obj.ctrl_path)
                [tif, obj.idims_ctrl] = LoadTif(obj.ctrl_path, obj.ctrl_chcount, [obj.ctrl_ch], verbosity);
                ch_dat = tif{obj.ctrl_ch, 1};
            end
        end
        
        function obj = updateZTrimParams(obj)
            %Trim
            trim_val = obj.ztrim;
            if obj.ztrim_auto > trim_val; trim_val = obj.ztrim_auto; end
            
            %Lower bound. (Index of lowest slice included)
            obj.z_min_apply = 1 + trim_val;
            if obj.z_min > obj.z_min_apply; obj.z_min_apply = obj.z_min; end
            
            %Upper bound. (Index of highest slice included)
            Z = 1;
            if ~isempty(obj.idims_sample)
                Z = obj.idims_sample.z;
            end
            
            obj.z_max_apply = Z - trim_val;
            if obj.z_max < obj.z_max_apply; obj.z_max_apply = obj.z_max; end
            
            %Make sure the min and max don't overlap.
            if obj.z_max_apply < obj.z_min_apply
                obj.z_max_apply = obj.z_min_apply;
                %You've included zero slices for some reason.
            end
        end
        
        function [obj, spots_table, coord_table] = loadZTrimmedTables_Sample(obj)
            [obj, spots_table, coord_table] = obj.loadZTrimmedTables('_spotTablesZTrimmed.mat', false);
        end
        
        function [obj, spots_table, coord_table] = loadZTrimmedTables_Control(obj)
            [obj, spots_table, coord_table] = obj.loadZTrimmedTables('_ctrlTablesZTrimmed.mat', true);
        end
        
        function [obj, spots_table, coord_table] = loadZTrimmedTables(obj, pathsfx, isctrl)
            %Look to see if it's been pre-saved...
            obj = obj.updateZTrimParams();
            tbl_path = [obj.out_stem pathsfx];
            if isfile(tbl_path)
                load(tbl_path, 'zmin', 'zmax');
                if zmin == obj.z_min_apply
                    if zmax == obj.z_max_apply
                        %Same trim. Load and return tables.
                        load(tbl_path, 'spots_table', 'coord_table');
                        return;
                    end
                end
            end
            
            %Need to recalculate.
            if isctrl
                if isempty(obj.ctrl_path)
                    spots_table = [];
                    coord_table = [];
                    return;
                end
                [obj, spots_table] = obj.loadControlSpotsTable();
                [obj, coord_table] = obj.loadControlCoordinateTable();
            else
                [obj, spots_table] = obj.loadSpotsTable();
                [obj, coord_table] = obj.loadCoordinateTable();
            end
            
            T = size(coord_table,1);
            masked_coord_table = cell(T,1);
            masked_spots_table = zeros(T,2);
            masked_spots_table(:,1) = spots_table(:,1);
            for t = 1:T
                tcoords = coord_table{t,1};
                [keeprows, ~] = find(tcoords(:,3) >= obj.z_min_apply);
                tcoords = tcoords(keeprows,:);
                [keeprows, ~] = find(tcoords(:,3) <= obj.z_max_apply);
                tcoords = tcoords(keeprows,:);
                masked_coord_table{t,1} = tcoords;
                masked_spots_table(t,2) = size(tcoords,1);
            end
            
            spots_table = masked_spots_table;
            coord_table = masked_coord_table;
            zmin = obj.z_min_apply;
            zmax = obj.z_max_apply;
            version = 1;
            save(tbl_path, 'version', 'zmin', 'zmax', 'spots_table', 'coord_table');
        end
    end
    
    methods(Static)
        
        function rnaspots_run = initDefault()
            rnaspots_run = RNASpotsRun;
            rnaspots_run.rna_ch = 0;
            rnaspots_run.light_ch = 0;
            rnaspots_run.total_ch = 0;
            rnaspots_run.ctrl_ch = 0;
            rnaspots_run.ctrl_chcount = 0;
            rnaspots_run.t_min = 10;
            rnaspots_run.t_max = 500;
            rnaspots_run.ztrim = 5;
            rnaspots_run.z_min = -1;
            rnaspots_run.z_max = -1;
            rnaspots_run.z_min_apply = -1;
            rnaspots_run.z_max_apply = -1;
            rnaspots_run.ttune_winsize = -1; %DEPR
            rnaspots_run.ttune_madfactor = 0.0; %DEPR
            rnaspots_run.intensity_threshold = 0;
            rnaspots_run.ttune_madf_min = -1.0;
            rnaspots_run.ttune_madf_max = 1.0;
            rnaspots_run.ttune_winsz_min = 3;
            rnaspots_run.ttune_winsz_max = 21;
            rnaspots_run.ttune_winsz_incr = 3;
            rnaspots_run.ttune_spline_itr = 3;
            rnaspots_run.overwrite_output = false;
            rnaspots_run.ttune_use_rawcurve = false;
            rnaspots_run.ttune_use_diffcurve = false;
            rnaspots_run.ttune_fit_to_log = true;
            rnaspots_run.ttune_reweight_fit = false;
            rnaspots_run.ttune_fit_strat = 0;
            rnaspots_run.ttune_thweight_med = 0.0;
            rnaspots_run.ttune_thweight_fit = 0.1;
            rnaspots_run.ttune_thweight_fisect = 0.9;
        end
        
        function rnaspots_run = loadFrom(path)
            if ~endsWith(path, '_rnaspotsrun.mat')
                path = [path '_rnaspotsrun.mat'];
            end
            load(path, 'rnaspots_run');
        end
        
    end
    
end