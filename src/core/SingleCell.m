%
%%
classdef SingleCell
    
    %%
    properties
        
        cell_number; %Index in cell mask
        cell_loc; %struct describing rectangle encompassing this cell in its original image
        dim_z; %Z dimension of source image
        
        mask_cell; %Boolean filter box.
        mask_nuc; %Same xy dims as cell box.
        mask_cyto; %Same xy dims as cell box
        
        spotcount_nuc;
        spotcount_cyto;
        spotcount_total;
        
        signal_nuc;
        signal_cyto;
        signal_total;
        
        nuc_ellip; %Vector of xy nucR for each z slice
        
        spots; %Obj array of RNASpot
        clouds; %Obj array of RNACloud
        
        %TEMP USAGE
        img_raw;
        img_nobkg;
        coord_list; %Relative to cell
    end
    
    methods
        
        %%
        function obj = preallocSpots(obj, count)
            if count < 1
                obj.spots = RNASpot.empty();
                return;
            end
            
            obj.spots(count) = RNASpot.newRNASpot();
            for i = count-1:-1:1
                obj.spots(i) = RNASpot.newRNASpot();
            end
        end
        
        %%
        function obj = findBoundaries(obj, cell_mask, nuc_mask)
            X = size(cell_mask,2);
            Y = size(cell_mask,1);
            Z = obj.dim_z; %One or both may be either 2D or 3D
            cellmask_3d = ndims(cell_mask) >= 3;
            
            if ~isempty(nuc_mask)
                nucmask_3d = ndims(nuc_mask) >= 3;
            else
                nucmask_3d = false;
            end
            
            cell_filter = logical(cell_mask == obj.cell_number); %Is there any reason this needs to be 16 bits and not boolean?
            if nnz(cell_filter) == 0; return; end
            
            if cellmask_3d
                rp = regionprops3(cell_filter,'BoundingBox','Area');
                cell_box = rp.BoundingBox;
                widx = 4; hidx = 5;
            else
                rp = regionprops(cell_filter,'BoundingBox','Area');
                cell_box = rp.BoundingBox;
                widx = 3; hidx = 4;
            end
                
            X0 = round(cell_box(1)) - 4;
            Y0 = round(cell_box(2)) - 4;
            X1 = round(cell_box(1)+ cell_box(widx)) + 4;
            Y1 = round(cell_box(2)+ cell_box(hidx)) + 4;
                
            if X0 < 1; X0 = 1; end
            if Y0 < 1; Y0 = 1; end
            if X1 > X; X1 = X; end
            if Y1 > Y; Y1 = Y; end
            
            if cellmask_3d
                obj.mask_cell = cell_filter(Y0:Y1,X0:X1,:);
                obj.cell_loc = SingleCell.generateRecPrismStruct(X0,X1,Y0,Y1,1,Z);
                
                if nucmask_3d
                    obj.mask_nuc = logical(nuc_mask(Y0:Y1,X0:X1,1:Z) .* obj.mask_cell);
                    obj.mask_cyto = obj.mask_cell & ~obj.mask_nuc;
                else
                    if ~isempty(nuc_mask)
                        cellmax = max(obj.mask_cell, [], 3);
                        obj.mask_nuc = logical(nuc_mask(Y0:Y1,X0:X1) .* cellmax);
                        obj.mask_cyto = cellmax & ~obj.mask_nuc;
                    end
                end
            else
                obj.mask_cell = cell_filter(Y0:Y1,X0:X1);
                obj.cell_loc = SingleCell.generateRecPrismStruct(X0,X1,Y0,Y1,1,Z);
                
                if nucmask_3d
                    cell_cyl = obj.get3DCellMask();
                    if isinteger(nuc_mask)
                        obj.mask_nuc = logical(nuc_mask(Y0:Y1,X0:X1,1:Z) .* uint16(cell_cyl));
                    else
                        obj.mask_nuc = logical(nuc_mask(Y0:Y1,X0:X1,1:Z) .* cell_cyl);
                    end
                    obj.mask_cyto = cell_cyl & ~obj.mask_nuc;
                else
                    if ~isempty(nuc_mask)
                        obj.mask_nuc = logical(nuc_mask(Y0:Y1,X0:X1) .* obj.mask_cell);
                        obj.mask_cyto = obj.mask_cell & ~obj.mask_nuc;
                    end
                end
            end
        end
        
        %%
        function cell_mask_3d = get3DCellMask(obj)
            if isempty(obj.mask_cell)
                cell_mask_3d = [];
                return;
            end
            if ndims(obj.mask_cell) >= 3
                cell_mask_3d = obj.mask_cell;
                return;
            end
            
            cell_cyl = zeros(obj.cell_loc.height_i, obj.cell_loc.width_i, obj.dim_z);
            for z = 1:obj.dim_z
                cell_cyl(:,:,z) = obj.mask_cell(:,:);
            end
            cell_mask_3d = logical(cell_cyl);
        end
        
        %%
        function obj = calculateNuclearEllipticity(obj)
            obj.nuc_ellip = [];
            if isempty(obj.mask_nuc); return; end
            Z = obj.dim_z;
            obj.nuc_ellip = NaN(1,Z);
            
            for z = 1:Z
                obj.nuc_ellip(z) = RNAQuant.findNucEllipticityAtZ(obj, z);
            end
        end
        
        %%
        function spot_count = getSpotCount(obj)
            if isempty(obj.spots)
                spot_count = 0;
            else
                spot_count = size(obj.spots,2);
            end
        end
        
        %%
        function cell_img = isolateCellBox(obj, input_image)
            cell_img = input_image(obj.cell_loc.top:obj.cell_loc.bottom, ...
               obj.cell_loc.left:obj.cell_loc.right, ...
               obj.cell_loc.z_bottom:obj.cell_loc.z_top);
        end
        
        %%
        function [cell_coord_table, nuc_coord_table] = getCoordsSubset(obj, coord_table_t)
            cell_coord_table = [];
            nuc_coord_table = [];
            if isempty(coord_table_t); return; end
            coord_dims = size(coord_table_t,2);
            
            x_valid = (coord_table_t(:,1) >= obj.cell_loc.left) & (coord_table_t(:,1) <= obj.cell_loc.right);
            y_valid = (coord_table_t(:,2) >= obj.cell_loc.top) & (coord_table_t(:,2) <= obj.cell_loc.bottom);
            incl_mtx = x_valid & y_valid;
            
            if coord_dims >= 3
                z_valid = (coord_table_t(:,3) >= obj.cell_loc.z_bottom) & (coord_table_t(:,3) <= obj.cell_loc.z_top);
                incl_mtx = incl_mtx & z_valid;
            end
            
            [incl_rows, ~] = find(incl_mtx);
            if isempty(incl_rows); return; end
            
            cell_coord_table = coord_table_t(incl_rows,:);
            cell_coord_table(:,1) = cell_coord_table(:,1) - obj.cell_loc.left + 1;
            cell_coord_table(:,2) = cell_coord_table(:,2) - obj.cell_loc.top + 1;
            if coord_dims >= 3
                cell_coord_table(:,3) = cell_coord_table(:,3) - obj.cell_loc.z_bottom + 1;
            end
            
            %Filter through cell mask
            c_count = size(cell_coord_table,1);
            mask_okay = false(c_count,1);
            if ndims(obj.mask_cell) > 2
                for i = 1:c_count
                    mask_okay(i,1) = obj.mask_cell(cell_coord_table(i,2), cell_coord_table(i,1), cell_coord_table(i,3));
                end
            else
                for i = 1:c_count
                    mask_okay(i,1) = obj.mask_cell(cell_coord_table(i,2), cell_coord_table(i,1));
                end
            end
            
            [incl_rows, ~] = find(mask_okay);
            if isempty(incl_rows)
                cell_coord_table = [];
                return;
            end
            cell_coord_table = cell_coord_table(incl_rows,:);

            if ~isempty(cell_coord_table)
                c_count = size(cell_coord_table,1);
                mask_okay = false(c_count,1);
                if ndims(obj.mask_nuc) > 2
                    for i = 1:c_count
                        mask_okay(i,1) = obj.mask_nuc(cell_coord_table(i,2), cell_coord_table(i,1), cell_coord_table(i,3));
                    end
                else
                    for i = 1:c_count
                        mask_okay(i,1) = obj.mask_nuc(cell_coord_table(i,2), cell_coord_table(i,1));
                    end
                end
            end
            [incl_rows, ~] = find(mask_okay);
             if isempty(incl_rows)
                nuc_coord_table = [];
                return;
            end
            nuc_coord_table = cell_coord_table(incl_rows,:);
        end
        
        %%
        function obj = updateSpotAndSignalValues(obj)
            %Reset.
            obj.spotcount_nuc = 0;
            obj.spotcount_cyto = 0;
            obj.spotcount_total = 0;
            obj.signal_nuc = 0.0;
            obj.signal_cyto = 0.0;
            obj.signal_total = 0.0;
            
            %Count spots
            spot_count = obj.getSpotCount();
            if spot_count > 0
                for s = 1:spot_count
                    this_spot = obj.spots(s);
                    if ~isempty(this_spot.gauss_fit)
                        if this_spot.gauss_fit.nucRNA
                            obj.spotcount_nuc = obj.spotcount_nuc + 1;
                        else
                            obj.spotcount_cyto = obj.spotcount_cyto + 1;
                        end
                        
                        %Add to signal.
                        if ~this_spot.in_cloud
                            if this_spot.gauss_fit.nucRNA
                                obj.signal_nuc = obj.signal_nuc + this_spot.gauss_fit.TotFitInt;
                            else
                                obj.signal_cyto = obj.signal_cyto + this_spot.gauss_fit.TotFitInt;
                            end
                        end
                    end
                end
            end
            
            %Clouds
            if ~isempty(obj.clouds)
                cloud_count = size(obj.clouds,2);
                for c = 1:cloud_count
                    this_cloud = obj.clouds(c);
                    if this_cloud.is_nuc
                        obj.signal_nuc = obj.signal_nuc + this_cloud.total_intensity;
                    else
                        obj.signal_cyto = obj.signal_cyto + this_cloud.total_intensity;
                    end
                end
            end
            
            %Update totals
            obj.spotcount_total = obj.spotcount_nuc + obj.spotcount_cyto;
            obj.signal_total = obj.signal_nuc + obj.signal_cyto;
            
        end
        
        %%
        function [obj, nuc_count, cyto_count, nuc_nascent_count] = estimateTargetCounts(obj, inclClouds)
            %TODO
            %Munsky B, Li G, Fox ZR, Shepherd DP, Neuert G. Distribution shapes govern the discovery of predictive models for gene regulation. Proc Natl Acad Sci U S A. 2018;115(29). doi:10.1073/pnas.1804060115

        end

    end
    
    methods(Static)
        
        %%
        function recstruct = generateRectangleStruct(x0, x1, y0, y1)
            recstruct = SingleCell.generateRecPrismStruct(x0, x1, y0, y1, 1, 1);
        end
        
        %%
        function recstruct = generateRecPrismStruct(x0, x1, y0, y1, z0, z1)
            recstruct = struct("left", min(x0,x1));
            recstruct.right = max(x0,x1);
            recstruct.top = min(y0,y1);
            recstruct.bottom = max(y0,y1);
            recstruct.z_bottom = min(z0,z1);
            recstruct.z_top = max(z0,z1);
            recstruct.width = abs(x1 - x0) + 1;
            recstruct.height = abs(y1 - y0) + 1;
            recstruct.depth = abs(z1 - z0) + 1;
            recstruct.width_i = round(recstruct.width);
            recstruct.height_i = round(recstruct.height);
            recstruct.depth_i = round(recstruct.depth);
        end
        
        %%
        function mycell = newRNACell(idx, spot_prealloc)
            if nargin < 1; idx = -1; end
            if nargin < 2; spot_prealloc = 0; end

            mycell = SingleCell;
            mycell.cell_number = idx;
            mycell.cell_loc = SingleCell.generateRectangleStruct(0,0,0,0);
            %mycell.nuc_loc = SingleCell.generateRectangleStruct(0,0,0,0);
            %mycell.nuc_loc_rel = SingleCell.generateRectangleStruct(0,0,0,0);
            
            mycell.spotcount_nuc = 0;
            mycell.spotcount_cyto = 0;
            mycell.spotcount_total = 0;
            
            if spot_prealloc > 0
                mycell = mycell.preallocSpots(spot_prealloc);
            else
                mycell.spots = RNASpot.empty();
            end
            mycell.clouds = RNACloud.empty();
            mycell.img_nobkg = [];
            mycell.coord_list = [];
            mycell.img_raw = [];
        end
        
    end
    
end