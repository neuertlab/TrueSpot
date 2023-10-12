%
%%

classdef RNASpotsRun
    
    properties
        img_name;
    
        intensity_threshold;
        threshold_results;

        paths;
        channels;
        th_params;
        dims;
        meta;
        options;
    end
    
    methods
        
        function obj = saveMe(obj)
            outpath = [obj.paths.out_namestem '_rnaspotsrun.mat'];
           
            run_info = struct('RNASpotsRunVersion', 2);

            run_info.img_name = obj.img_name;
            run_info.intensity_threshold = obj.intensity_threshold;
            run_info.threshold_results = obj.threshold_results;
            run_info.paths = obj.paths;
            run_info.channels = obj.channels;
            run_info.th_params = obj.th_params;
            run_info.dims = obj.dims;
            run_info.meta = obj.meta;
            run_info.options = obj.options;

            save(outpath, 'run_info');
        end
        
        function out_stem = getFullOutStem(obj)
            out_stem = [obj.paths.out_dir filesep obj.paths.out_namestem];
        end

        function out_stem = getFullCtrlOutStem(obj)
            out_stem = [obj.paths.out_dir filesep obj.paths.ctrl_out_namestem];
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
            tbl_path = [obj.getFullOutStem() '_spotTable.mat'];
            if isfile(tbl_path)
                load(tbl_path, 'spot_table');
            end
        end
        
        function [obj, spot_table] = loadControlSpotsTable(obj)
            spot_table = [];
            tbl_path = [obj.getFullCtrlOutStem() '_spotTable.mat'];
            if isfile(tbl_path)
                load(tbl_path, 'spot_table');
            end
        end
        
        function [obj, call_table] = loadCallTable(obj)
            tbl_path_RNA = [obj.getFullOutStem() '_callTable.mat'];
            if isfile(tbl_path_RNA)
                load(tbl_path_RNA, 'call_table');
            else
                call_table = [];
            end
        end

        function [obj, call_table] = loadControlCallTable(obj)
            tbl_path_RNA = [obj.getFullCtrlOutStem() '_callTable.mat'];
            if isfile(tbl_path_RNA)
                load(tbl_path_RNA, 'call_table');
            else
                call_table = [];
            end
        end

        function [obj, background_mask] = loadBackgroundMask(obj)
            background_mask = [];
            if isfile([obj.paths.bkg_mask_path '.mat'])
                load(obj.paths.bkg_mask_path, 'background_mask');
            end
        end
        
        function [obj, my_images] = loadImageViewStructs(obj)
            my_images = [];
            mypath = [obj.getFullOutStem() '_imgviewstructs.mat'];
            if isfile(mypath)
                load(mypath, 'my_images');
            end
        end
        
        function [obj, img_filter] = loadFilteredImage(obj)
            img_filter = [];
            mypath = [obj.getFullOutStem() '_prefilteredIMG.mat'];
            if isfile(mypath)
                load(mypath, 'img_filter');
            end
        end
        
        function obj = updateBackgroundFilteredCoords(obj)
            [~, background_mask] = obj.loadBackgroundMask();
            [~, call_table] = obj.loadCallTable();
            [~, th_list] = obj.loadThresholdTable();

            if isempty(background_mask)
                return;
            end

            if ndims(background_mask) < 3
                keep_me = background_mask(call_table{:,'isnap_y'}, call_table{:,'isnap_x'});
            else
                keep_me = background_mask(call_table{:,'isnap_y'}, call_table{:,'isnap_x'}, call_table{:,'isnap_z'});
            end

            T = size(th_list,2);
            spot_table = NaN(T,2);
            spot_table(:,1) = th_list(:);

            if nnz(keep_me) > 0
                keep_idx = find(keep_me);
                call_table = call_table(keep_idx, :);
                for t = 1:T
                    spot_table(t,2) = nnz(call_table{:,'dropout_thresh'} >= spot_table(t,1));
                end
            else
                spot_table(:,2) = 0;
                call_table = table.empty();
            end

            save([obj.paths.bkg_filter_stem '_callTable'], 'call_table');
            save([obj.paths.bkg_filter_stem '_spotTable'], 'spot_table');

            %Image structs
            [~, my_images] = obj.loadImageViewStructs();
            my_images(1).image = immultiply(my_images(1).image, background_mask);
            my_images(2).image = immultiply(my_images(2).image, background_mask);
            save([obj.paths.bkg_filter_stem '_imgviewstructs'], 'my_images');
        end

        function obj = updateImageDimensions(obj)
            out_stem = obj.getFullOutStem();
            ctrl_stem = obj.getFullCtrlOutStem();
            if isempty(obj.dims.idims_sample)
                runinfo_path = [out_stem '_runinfo'];
                if isfile(runinfo_path)
                    %Load the run info.
                    load(runinfo_path, 'idims');
                else
                    %Reload TIF and count manually.
                    [X,Y,Z] = GetTifDims(obj.paths.img_path, obj.channels.total_ch);
                    idims = struct('x', X, 'y', Y, 'z', Z); 
                end
                %idims
                obj.dims.idims_sample = idims;
            end
            
            if isempty(obj.dims.idims_ctrl)
                if ~isempty(obj.paths.ctrl_img_path)
                    runinfo_path = [ctrl_stem '_runinfo'];
                    if isfile(runinfo_path)
                        %Load the run info.
                        load(runinfo_path, 'idims');
                    else
                        %Reload TIF and count manually.
                        [X,Y,Z] = GetTifDims(obj.paths.ctrl_img_path, obj.channels.ctrl_chcount);
                        idims = struct('x', X, 'y', Y, 'z', Z);
                    end
                    obj.dims.idims_ctrl = idims;
                else
                    %Assumed to be background
                    obj.dims.idims_sample = idims;
                end
            end
        end
        
        function [obj, tif] = loadSampleTif(obj, verbosity)
            if (nargin < 2); verbosity = 2; end
            if obj.channels.light_ch > 0
                [tif, obj.dims.idims_sample] = LoadTif(obj.paths.img_path, obj.channels.total_ch, [obj.channels.rna_ch, obj.channels.light_ch], verbosity);
            else
                [tif, obj.dims.idims_sample] = LoadTif(obj.paths.img_path, obj.channels.total_ch, [obj.channels.rna_ch], verbosity);
            end
        end
        
        function [obj, ch_dat] = loadControlChannel(obj, verbosity)
            ch_dat = [];
            if (nargin < 2); verbosity = 2; end
            if ~isempty(obj.paths.ctrl_img_path)
                [tif, obj.dims.idims_ctrl] = LoadTif(obj.paths.ctrl_img_path, obj.channels.ctrl_chcount, [obj.channels.ctrl_ch], verbosity);
                ch_dat = tif{obj.channels.ctrl_ch, 1};
            end
        end
        
        function obj = updateZTrimParams(obj)
            %Upper bound. (Index of highest slice included)
            Z = 1;
            if ~isempty(obj.dims.idims_sample)
                Z = obj.dims.idims_sample.z;
            end
            
            %Update zmin and max if unset.
            if obj.dims.z_max < 1
                obj.dims.z_max = Z;
            end
            if obj.dims.z_max < obj.dims.z_min
                obj.dims.z_max = obj.dims.z_min;
            end
            if obj.dims.z_max > obj.dims.idims_sample.z
                obj.dims.z_max = obj.dims.idims_sample.z;
            end
            if obj.dims.z_min < 1
                obj.dims.z_min = 1;
            end
            
            %Trim
            trim_val = obj.dims.ztrim;
            if obj.dims.ztrim_auto > trim_val; trim_val = obj.dims.ztrim_auto; end
            
            %Lower bound. (Index of lowest slice included)
            obj.dims.z_min_apply = 1 + trim_val;
            if obj.dims.z_min > obj.dims.z_min_apply; obj.dims.z_min_apply = obj.dims.z_min; end
            
            obj.dims.z_max_apply = Z - trim_val;
            if obj.dims.z_max < obj.dims.z_max_apply; obj.dims.z_max_apply = obj.dims.z_max; end
            
            %Make sure the min and max don't overlap.
            if obj.dims.z_max_apply < obj.dims.z_min_apply
                obj.dims.z_max_apply = obj.dims.z_min_apply;
                %You've included zero slices for some reason.
            end
        end

        function obj = deleteSavedZTrimmedTables(obj)
            out_stem = obj.getFullOutStem();
            tbl_path = [out_stem '_spotTablesZTrimmed.mat'];
            if isfile(tbl_path)
                delete(tbl_path);
            end
            tbl_path = [out_stem '_ctrlTablesZTrimmed.mat'];
            if isfile(tbl_path)
                delete(tbl_path);
            end
        end
        
        function [obj, spots_table, call_table] = loadZTrimmedTables_Sample(obj)
            [obj, spots_table, call_table] = obj.loadZTrimmedTables(false);
        end
        
        function [obj, spots_table, call_table] = loadZTrimmedTables_Control(obj)
            [obj, spots_table, call_table] = obj.loadZTrimmedTables(true);
        end
        
        function [obj, spots_table, call_table] = loadZTrimmedTables(obj, isctrl)
            obj = obj.updateZTrimParams();

            if isctrl
                if isempty(obj.paths.ctrl_img_path)
                    spots_table = [];
                    call_table = [];
                    return;
                end
                [obj, spots_table] = obj.loadControlSpotsTable();
                [obj, call_table] = obj.loadControlCallTable();
            else
                [obj, spots_table] = obj.loadSpotsTable();
                [obj, call_table] = obj.loadCallTable();
            end

            zgood = call_table{:, 'isnap_z'} >= obj.dims.z_min_apply;
            zgood = and(zgood, call_table{:, 'isnap_z'} <= obj.dims.z_max_apply);
            if nnz(zgood) > 0
                call_table = call_table{find(zgood), :};
                T = size(spots_table,1);
                for t = 1:T
                    spots_table(t,2) = nnz(call_table{:, 'dropout_thresh'} >= spots_table(t,1));
                end
            else
                spots_table(:,2) = 0;
                call_table = table.empty();
            end

        end

    end
    
    methods(Static)
        
        function rnaspots_run = initDefault()
            rnaspots_run = RNASpotsRun;
            rnaspots_run.img_name = '';
            rnaspots_run.intensity_threshold = 0;
            rnaspots_run.threshold_results = struct.empty();

            rnaspots_run.paths = struct('img_path', '');
            rnaspots_run.paths.out_dir = '';
            rnaspots_run.paths.out_namestem = '';
            rnaspots_run.paths.ctrl_img_path = '';
            rnaspots_run.paths.ctrl_out_namestem = '';
            rnaspots_run.paths.bkg_mask_path = '';
            rnaspots_run.paths.bkg_filter_stem = '';
            rnaspots_run.paths.cellseg_path = '';
            rnaspots_run.paths.csv_out_path = '';
            rnaspots_run.paths.params_out_path = '';

            rnaspots_run.channels = struct('rna_ch', 0);
            rnaspots_run.channels.light_ch = 0;
            rnaspots_run.channels.total_ch = 0;
            rnaspots_run.channels.ctrl_ch = 0;
            rnaspots_run.channels.ctrl_chcount = 0;
            
            rnaspots_run.th_params = RNAThreshold.genEmptyThresholdParamStruct();

            wmin = 3; wmax = 21; wincr = 3;
            rnaspots_run.th_params.window_sizes = [wmin:wincr:wmax];

            rnaspots_run.options = struct('dtune_gaussrad', 7);
            rnaspots_run.options.overwrite_output = false;
            rnaspots_run.options.deadpix_detect = true;
            rnaspots_run.options.csv_zero_based_coords = false;
            rnaspots_run.options.debug_level = 0;
            rnaspots_run.options.use_max_proj = false;
            rnaspots_run.options.threads = 1;
            rnaspots_run.options.winsize_min = wmin;
            rnaspots_run.options.winsize_max = wmax;
            rnaspots_run.options.winsize_incr = wincr;
            rnaspots_run.options.t_min = 10;
            rnaspots_run.options.t_max = 500;

            rnaspots_run.dims = struct('ztrim', 0);
            rnaspots_run.dims.ztrim_auto = 0;
            rnaspots_run.dims.z_min = -1;
            rnaspots_run.dims.z_max = -1;
            rnaspots_run.dims.z_min_apply = -1;
            rnaspots_run.dims.z_max_apply = -1;
            rnaspots_run.dims.idims_sample = struct('x', -1, 'y', -1, 'z', -1);
            rnaspots_run.dims.idims_ctrl = struct('x', -1, 'y', -1, 'z', -1);
      
            rnaspots_run.meta = struct('type_probe', '');
            rnaspots_run.meta.type_species = '';
            rnaspots_run.meta.type_cell = '';
            rnaspots_run.meta.type_target = '';
            rnaspots_run.meta.type_targetmol = '';
            rnaspots_run.meta.idims_voxel = struct('x', -1, 'y', -1, 'z', -1);
            rnaspots_run.meta.idims_expspot = struct('x', -1, 'y', -1, 'z', -1);
        end
        
        function rnaspots_run = loadFrom(path, updateOutDir)
            if nargin < 2
                updateOutDir = false;
            end

            if ~endsWith(path, '_rnaspotsrun.mat')
                path = [path '_rnaspotsrun.mat'];
            end
            if isfile(path)
                %load(path, 'rnaspots_run');
                load(path, 'run_info');

                rnaspots_run = RNASpotsRun;
                rnaspots_run.img_name = run_info.img_name;
                rnaspots_run.intensity_threshold = run_info.intensity_threshold;
                rnaspots_run.threshold_results = run_info.threshold_results;
                rnaspots_run.paths = run_info.paths;
                rnaspots_run.channels = run_info.channels;
                rnaspots_run.th_params = run_info.th_params;
                rnaspots_run.dims = run_info.dims;
                rnaspots_run.meta = run_info.meta;
                rnaspots_run.options = run_info.options;

                %Change output dir
                if updateOutDir
                    [currentdir, ~, ~] = fileparts(path);
                    rnaspots_run.paths.out_dir = currentdir;
                end
            else
                rnaspots_run = [];
            end
        end
        
        function call_table = saveCallTable(call_table, save_stem)
            %Spot coords
            save([save_stem '_callTable.mat'], 'call_table', '-v7.3');
        end

        function [spot_table, call_table] = saveTables(spot_table, call_table, save_stem)
            spot_table = double(spot_table);

            %Spot count table
            save([save_stem '_spotTable.mat'], 'spot_table', '-v7.3');
            %Spot coords
            save([save_stem '_callTable.mat'], 'call_table', '-v7.3');
        end

    end
    
end