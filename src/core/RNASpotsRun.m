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
        
        t_min;
        t_max;
        ztrim;
        ztrim_auto;
        ttune_winsize;
        ttune_wscorethresh;
        intensity_threshold;
        
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
            tbl_path = [obj.ctrl_stem 'spot_table.mat'];
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
        
    end
    
    methods(Static)
        
        function rnaspots_run = loadFrom(path)
            if ~endsWith(path, '_rnaspotsrun.mat')
                path = [path '_rnaspotsrun.mat'];
            end
            load(path, 'rnaspots_run');
        end
        
    end
    
end