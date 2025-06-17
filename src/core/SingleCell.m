%
%%
classdef SingleCell

    %TODO Anything that used the spots field needs to be updated!
    
    %%
    properties
        
        cell_number; %Index in cell mask
        cell_loc; %struct describing rectangle encompassing this cell in its original image
        dim_z; %Z dimension of source image
        cell_stats;
        
        mask_cell; %Boolean filter box.
        mask_nuc; %Same xy dims as cell box.
        mask_cyto; %Same xy dims as cell box
        
        spotcount_nuc;
        spotcount_cyto;
        spotcount_total;
        
        signal_nuc;
        signal_cyto;
        signal_total;

        %Target molecule count estimates
        nucCount;
        nucNascentCount;
        nucCloud;
        nucNascentCloud;
        cytoCount;
        cytoCloud;
        
        nuc_ellip; %Vector of xy nucR for each z slice
        
        spots; %Obj array of RNASpot -- DEPRECATED left here for compatibility
        clouds; %Obj array of RNACloud
        
        spotTable;
        spotZFits; %Cell array stores the fits for each z slice

        %TEMP USAGE
        img_raw;
        img_nobkg;
        coord_list; %Relative to cell
    end
    
    methods
        
        %%
        function obj = convertSpotStorage(obj)
            obj.spotTable = SingleCell.spotObjArr2Table(obj.spots);
            spotList = obj.spots;

            spotCount = size(obj.spotTable, 1);
            obj.spotZFits = cell(1, spotCount);
            for i = 1:spotCount
                obj.spotZFits{i} = spotList(i).gfit_slices;
            end
            obj.spots = [];
        end

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
            if ~isempty(obj.spotTable)
                spot_count = size(obj.spotTable, 1);
            else
                if isempty(obj.spots)
                    spot_count = 0;
                else
                    spot_count = size(obj.spots,2);
                end
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
        function [cyto_calls, nuc_calls] = getCellCalls(obj, call_table)
            if isempty(call_table); return; end
            cyto_calls = [];
            nuc_calls = [];

            %1. convert all coordinates to cell local and toss those
            %   outside box
            xx = int32(call_table{:, 'isnap_x'}) - obj.cell_loc.left + 1;
            yy = int32(call_table{:, 'isnap_y'}) - obj.cell_loc.top + 1;
            zz = int32(call_table{:, 'isnap_z'}) - obj.cell_loc.z_bottom + 1;

            xOkay = and(xx > 0, xx <= obj.cell_loc.width);
            yOkay = and(yy > 0, yy <= obj.cell_loc.height);
            zOkay = and(zz > 0, zz <= obj.cell_loc.depth);
            keep = and(xOkay, and(yOkay, zOkay));
            if nnz(keep) < 1; return; end

            xx = xx(keep,:);
            yy = yy(keep,:);
            zz = zz(keep,:);
            call_table = call_table(keep,:);

            %Put thru masks

            if ndims(obj.mask_nuc) > 2
                inNuc = RNAUtils.isInMask3(obj.mask_nuc, xx, yy, zz);
            else
                inNuc = RNAUtils.isInMask2(obj.mask_nuc, xx, yy);
            end
            if nnz(inNuc > 0)
                nuc_calls = call_table(inNuc, :);
            end

            if ndims(obj.mask_cell) > 2
                inCell = RNAUtils.isInMask3(obj.mask_cell, xx, yy, zz);
            else
                inCell = RNAUtils.isInMask2(obj.mask_cell, xx, yy);
            end
            inCyto = and(inCell, ~inNuc);
            if nnz(inCyto > 0)
                cyto_calls = call_table(inCyto, :);
            end

        end

        %%
        function obj = updateSpotAndSignalValues(obj, globalBrightTh, globalSingleInt, removeLikelyDups, appliedThreshold)
            if nargin < 2; globalBrightTh = 65535.0; end %Used if not enough spots in this cell to work with.
            if nargin < 3; globalSingleInt = 200.0; end
            if nargin < 4; removeLikelyDups = true; end
            if nargin < 5; appliedThreshold = 0; end

            %Reset.
            obj.spotcount_nuc = 0;
            obj.spotcount_cyto = 0;
            obj.spotcount_total = 0;
            obj.signal_nuc = 0.0;
            obj.signal_cyto = 0.0;
            obj.signal_total = 0.0;

            obj.nucCount= 0;
            obj.nucNascentCount = 0;
            obj.nucCloud = 0;
            obj.nucNascentCloud = 0;
            obj.cytoCount = 0;
            obj.cytoCloud = 0;

            %If not stored in table, convert to table
            if isempty(obj.spotTable) & ~isempty(obj.spots)
                obj = obj.convertSpotStorage();
            end

            if isempty(obj.spotTable)
                return;
            end

            %Cluster spots, if applicable
            obj = obj.clusterLikelyDuplicateSpots();
            
            %Count spots
            spot_count = obj.getSpotCount();
            passTh = true(spot_count, 1);
            if appliedThreshold > 0
                passTh = obj.spotTable{:, 'dropout_thresh'} >= appliedThreshold;
            end
            if spot_count > 0
                nucFlag = obj.spotTable{:, 'nucRNA'};
                obj.spotcount_nuc = nnz(nucFlag);
                obj.spotcount_cyto = nnz(~nucFlag);

                cloudFlag = obj.spotTable{:, 'in_cloud'};
                nucNoCloud = and(~cloudFlag, nucFlag);
                nucNoCloud = and(nucNoCloud, passTh);
                cytoNoCloud = and(~cloudFlag, ~nucFlag);
                cytoNoCloud = and(cytoNoCloud, passTh);
                ttfit = obj.spotTable{nucNoCloud, 'TotFitInt'};
                obj.signal_nuc = sum(ttfit, 'all');
                ttfit = obj.spotTable{cytoNoCloud, 'TotFitInt'};
                obj.signal_cyto = sum(ttfit, 'all');

                clear nucFlag cloudFlag nucNoCloud cytoNoCloud ttfit
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
            
            obj = obj.updateCountEstimates(2, globalBrightTh, globalSingleInt, removeLikelyDups, appliedThreshold);
        end
        
        %%
        function obj = updateCountEstimates(obj, brightStDevs, globalBrightTh, globalSingleInt, removeLikelyDups, appliedThreshold)
            if nargin < 2; brightStDevs = 2; end
            if nargin < 3; globalBrightTh = 65535.0; end %Used if not enough spots in this cell to work with.
            if nargin < 4; globalSingleInt = 200.0; end
            if nargin < 5; removeLikelyDups = true; end
            if nargin < 6; appliedThreshold = 0; end

            %Munsky B, Li G, Fox ZR, Shepherd DP, Neuert G. Distribution shapes govern the discovery of predictive models for gene regulation. Proc Natl Acad Sci U S A. 2018;115(29). doi:10.1073/pnas.1804060115
            obj.nucCount = 0;
            obj.nucNascentCount = 0;
            obj.cytoCount = 0;

            %If not stored in table, convert to table
            if isempty(obj.spotTable) & ~isempty(obj.spots)
                obj = obj.convertSpotStorage();
            end

            %Determine which spots to count
            isCounted = true(1, size(obj.spotTable, 1));
            if appliedThreshold > 0
                passTh = obj.spotTable{:, 'dropout_thresh'} >= appliedThreshold;
                isCounted(~passTh) = false;
            end

            %Remove duplicates
            if ~isempty(obj.spotTable)
                if removeLikelyDups
                    dup_flag = obj.spotTable{:, 'likely_dup'} > 0;
                    is_center = (obj.spotTable{:, 'likely_dup'}) == (obj.spotTable{:, 'uid'});
                    is_dup = and(~is_center, dup_flag)';

                    isCounted(is_dup) = false;

                    clear is_center dup_flag
                end
            end
            
            useSpotTable = obj.spotTable(isCounted, :);

            %Get "mature RNA" (single target) intensity
            %Collect spot intensities
            totalSpots = size(useSpotTable, 1);
            if(totalSpots > 0)
                spotints = useSpotTable{:, 'TotFitInt'};
                inNuc = useSpotTable{:, 'nucRNA'};

                if totalSpots > 2
                    iMean = mean(spotints, 'all', 'omitnan');
                    %iMed = median(spotints, 'all', 'omitnan');
                    iStd = std(spotints, 0, 'all', 'omitnan');
                    brightnessThreshold = iMean + (brightStDevs * iStd);
                    %brightnessThreshold = iMed + (brightStDevs * iStd);
                else
                    brightnessThreshold = globalBrightTh;
                end

                isTooBright = (spotints >= brightnessThreshold);
                if totalSpots > 2
                    normSpots = spotints(~isTooBright);
                    singleIntensity = median(normSpots, 'all', 'omitnan');
                    clear normSpots
                else
                    singleIntensity = globalSingleInt;
                end

                %For now, all not-too-bright spots are counted as 1?
                counts = ones(1, totalSpots);
                counts(isTooBright) = round(spotints(isTooBright) ./ singleIntensity);

                obj.nucCount = sum(counts(inNuc));
                obj.nucNascentCount = sum(counts(inNuc & isTooBright));
                obj.cytoCount = sum(counts(~inNuc));

                useSpotTable{:, 'nascent_flag'} = and(isTooBright, inNuc);
            end

            totalClouds = size(obj.clouds, 2);
            if(totalClouds > 0)
                cloudList = obj.clouds;
                cloudInts = [cloudList.total_intensity];

                %Redo spots, omitting spots marked as part of clouds
                if(totalSpots > 0)
                    inCloud = useSpotTable{:, 'in_cloud'};
                    if nnz(inCloud) > 0
                        obj.nucCloud = sum(counts(inNuc & ~inCloud));
                        obj.nucNascentCloud = sum(counts(inNuc & isTooBright & ~inCloud));
                        obj.cytoCloud = sum(counts(~inNuc & ~inCloud));
                    else
                        obj.nucCloud = obj.nucCount;
                        obj.nucNascentCloud = obj.nucNascentCount;
                        obj.cytoCloud = obj.cytoCount;
                    end
                else
                    obj.nucCloud = 0;
                    obj.nucNascentCloud = 0;
                    obj.cytoCloud = 0;

                    %Need a singleIntensity and brightnessThreshold
                    iMean = mean(cloudInts, 'all', 'omitnan');
                    %iMed = median(cloudInts, 'all', 'omitnan');
                    iStd = std(cloudInts, 0, 'all', 'omitnan');
                    brightnessThreshold = iMean + (brightStDevs * iStd);
                    %brightnessThreshold = iMed + (brightStDevs * iStd);

                    %The single will just be the smallest cloud
                    singleIntensity = min(cloudInts, [], 'all', 'omitnan');
                end

                %Try to get target count from clouds.
                cloudNuc = [cloudList.is_nuc];
                cloudTooBright = cloudInts >= brightnessThreshold;
                cloudCounts = round(cloudInts ./ singleIntensity);

                obj.nucCloud = obj.nucCloud + sum(cloudCounts(cloudNuc));
                obj.nucNascentCloud = obj.nucNascentCloud + sum(cloudCounts(cloudNuc & cloudTooBright));
                obj.cytoCloud = obj.cytoCloud + sum(cloudCounts(~cloudNuc));

                for j = 1:totalClouds
                    obj.clouds(j).nascent_flag = cloudTooBright(j);
                end
            else
                obj.nucCloud = obj.nucCount;
                obj.nucNascentCloud = obj.nucNascentCount;
                obj.cytoCloud = obj.cytoCount;
            end

            %Update object table
            obj.spotTable{isCounted, 'nascent_flag'} = useSpotTable{:, 'nascent_flag'};
        end

        %%
        function obj = clusterLikelyDuplicateSpots(obj, mergeRadXY, mergeRadZ)
            if nargin < 2; mergeRadXY = 2; end
            if nargin < 3; mergeRadZ = 1; end

            if isempty(obj.spotTable); return; end

            spotCount = size(obj.spotTable, 1);

            rxysq = double(mergeRadXY) ^ 2;
            rzsq = double(mergeRadZ) ^ 2;
            mergeRad3 = sqrt(rxysq + rxysq + rzsq);

            %Update spot table if needed
            colnames = obj.spotTable.Properties.VariableNames;
            if ~ismember('likely_dup', colnames)
                obj.spotTable{:, 'uid'} = [1:1:spotCount]';
            end
            obj.spotTable{:, 'likely_dup'} = uint32(zeros(spotCount, 1)); %Reset all to 0

            %Matrix
            callmtx = NaN(spotCount, 3);
            callmtx(:,1) = obj.spotTable{:, 'xinit'};
            callmtx(:,2) = obj.spotTable{:, 'yinit'};
            callmtx(:,3) = obj.spotTable{:, 'zinit'};

            cand_table = RNACoords.findMatchCandidates(callmtx, callmtx, mergeRad3, mergeRadZ);
            
            %Remove any self-matches
            %dist3 row col dist2 distz
            cand_table = cand_table((cand_table(:,1) > 0), :);
            if isempty(cand_table); return; end
        
            %Find clusters
            gg = digraph(cand_table(:,2), cand_table(:,3));
            gbins = conncomp(gg);
            binnedSpots = size(gbins, 2);
            clusterCount = max(gbins, [], 'all', 'omitnan');
            for c = 1:clusterCount
                in_cluster = (gbins == c);
                nodeCount = nnz(in_cluster);
                if nodeCount > 1
                    %find the biggest member to be "center" (using max
                    %total intensity)
                    totExpInt = immultiply(obj.spotTable{1:binnedSpots, 'TotExpInt'}', in_cluster);
                    [~, idx] = max(totExpInt, [], 'all', 'omitnan');
                    nodeId = obj.spotTable{idx, 'uid'};
                    obj.spotTable{in_cluster, 'likely_dup'} = nodeId;
                end
            end
        end

        %%
        function cloudList = getAllClouds(obj, nucFlag, nascentFlag)
            if nargin < 2; nucFlag = -1; end
            if nargin < 3; nascentFlag = -1; end

            if isempty(obj.clouds)
                cloudList = struct.empty();
                return;
            end

            myclouds = obj.clouds;
            ccount = size(myclouds, 2);

            nn = [myclouds.nascent_flag];
            if nascentFlag == 0
                nascentOkay = ~nn;
            elseif nascentFlag == 1
                nascentOkay = nn;
            else
                nascentOkay = true(1, ccount);
            end
            clear nn;

            nn = [myclouds.is_nuc];
            if nucFlag == 0
                nucOkay = ~nn;
            elseif nucFlag == 1
                nucOkay = nn;
            else
                nucOkay = true(1, ccount);
            end
            clear nn;

            keepSet = and(nascentOkay, nucOkay);
            cloudList = myclouds(keepSet);

        end

        %%
        function obj = getBasicStats(obj, cell_img)
            obj.cell_stats = struct();
            cell_img = double(cell_img);

            %Nucleus
            if ndims(obj.mask_nuc) > 2
                % 3D
                cnuc3 = immultiply(cell_img, obj.mask_nuc);
                cnuc3(cnuc3 == 0) = NaN;
                obj.cell_stats.nuc_vol = nnz(obj.mask_nuc);
                obj.cell_stats.nuc_intensity_total = sum(cnuc3, 'all', 'omitnan');
                obj.cell_stats.nuc_intensity_mean = mean(cnuc3, 'all', 'omitnan');
                obj.cell_stats.nuc_intensity_median = median(cnuc3, 'all', 'omitnan');
                obj.cell_stats.nuc_intensity_stdev = std(cnuc3, 0, 'all', 'omitnan');

                % Max proj
                cmax = max(cell_img, [], 3, 'omitnan');
                nuc2 = max(obj.mask_nuc, [], 3, 'omitnan');
                cnucmax = immultiply(cmax, nuc2);
                cnucmax(cnucmax == 0) = NaN;
                obj.cell_stats.nuc_max_area = nnz(nuc2);
                obj.cell_stats.nuc_max_intensity_total = sum(cnucmax, 'all', 'omitnan');
                obj.cell_stats.nuc_max_intensity_mean = mean(cnucmax, 'all', 'omitnan');
                obj.cell_stats.nuc_max_intensity_median = median(cnucmax, 'all', 'omitnan');
                obj.cell_stats.nuc_max_intensity_stdev = std(cnucmax, 0, 'all', 'omitnan');

                % Per slice
                slice_count = size(obj.mask_nuc, 3);
                obj.cell_stats.nuc_slices_stats = cell(1, slice_count);
                for z = 1:slice_count
                    stat_struct = struct();
                    cnucz = immultiply(cell_img(:,:,z), obj.mask_nuc(:,:,z));
                    cnucz(cnucz == 0) = NaN;

                    stat_struct.nuc_slice_area = nnz(cnucz);
                    stat_struct.nuc_slice_intensity_total = sum(cnucz, 'all', 'omitnan');
                    stat_struct.nuc_slice_intensity_mean = mean(cnucz, 'all', 'omitnan');
                    stat_struct.nuc_slice_intensity_median = median(cnucz, 'all', 'omitnan');
                    stat_struct.nuc_slice_intensity_stdev = std(cnucz, 0, 'all', 'omitnan');

                    obj.cell_stats.nuc_slices_stats{z} = stat_struct;
                end

            else
                %Just the 2D
                cmax = max(cell_img, [], 3, 'omitnan');
                cnuc2 = immultiply(cmax, obj.mask_nuc);
                cnuc2(cnuc2 == 0) = NaN;
                obj.cell_stats.nuc_area = nnz(obj.mask_nuc);
                obj.cell_stats.nuc_intensity_total = sum(cnuc2, 'all', 'omitnan');
                obj.cell_stats.nuc_intensity_mean = mean(cnuc2, 'all', 'omitnan');
                obj.cell_stats.nuc_intensity_median = median(cnuc2, 'all', 'omitnan');
                obj.cell_stats.nuc_intensity_stdev = std(cnuc2, 0, 'all', 'omitnan');
            end

            %Cyto/Full cell
            if ndims(obj.mask_cyto) > 2
                if ndims(obj.mask_cell) > 2
                    obj.cell_stats.cell_vol = nnz(obj.mask_cell);
                    cell_mask2 = max(obj.mask_cell, [], 3, 'omitnan');
                else
                    obj.cell_stats.cell_area = nnz(obj.mask_cell);
                    cell_mask2 = obj.mask_cell;
                end
                obj.cell_stats.cyto_vol = nnz(obj.mask_cyto);

                ccyto3 = immultiply(cell_img, obj.mask_cyto);
                ccyto3(ccyto3 == 0) = NaN;
                obj.cell_stats.cyto_intensity_total = sum(ccyto3, 'all', 'omitnan');
                obj.cell_stats.cyto_intensity_mean = mean(ccyto3, 'all', 'omitnan');
                obj.cell_stats.cyto_intensity_median = median(ccyto3, 'all', 'omitnan');
                obj.cell_stats.cyto_intensity_stdev = std(ccyto3, 0, 'all', 'omitnan');

                % Max proj (use inverse of nuc mask max proj)
                if ndims(obj.mask_nuc) > 2
                    nuc_mask2 = max(obj.mask_nuc, [], 3, 'omitnan');
                else
                    nuc_mask2 = obj.mask_nuc;
                end
                cmax = max(cell_img, [], 3, 'omitnan');
                cytomask2 = and(cell_mask2, ~nuc_mask2);
                ccmax = immultiply(cmax, cytomask2);
                ccmax(ccmax == 0) = NaN;
                obj.cell_stats.cyto_max_area = nnz(ccmax);
                obj.cell_stats.cyto_max_intensity_total = sum(ccmax, 'all', 'omitnan');
                obj.cell_stats.cyto_max_intensity_mean = mean(ccmax, 'all', 'omitnan');
                obj.cell_stats.cyto_max_intensity_median = median(ccmax, 'all', 'omitnan');
                obj.cell_stats.cyto_max_intensity_stdev = std(ccmax, 0, 'all', 'omitnan');

                % Per slice
                slice_count = size(obj.mask_cyto, 3);
                obj.cell_stats.cyto_slices_stats = cell(1, slice_count);
                for z = 1:slice_count
                    stat_struct = struct();
                    ccz = immultiply(cell_img(:,:,z), obj.mask_cyto(:,:,z));
                    ccz(ccz == 0) = NaN;

                    stat_struct.cyto_slice_area = nnz(ccz);
                    stat_struct.cyto_slice_intensity_total = sum(ccz, 'all', 'omitnan');
                    stat_struct.cyto_slice_intensity_mean = mean(ccz, 'all', 'omitnan');
                    stat_struct.cyto_slice_intensity_median = median(ccz, 'all', 'omitnan');
                    stat_struct.cyto_slice_intensity_stdev = std(ccz, 0, 'all', 'omitnan');

                    obj.cell_stats.cyto_slices_stats{z} = stat_struct;
                end
            else
                obj.cell_stats.cell_area = nnz(obj.mask_cell);
                obj.cell_stats.cyto_area = nnz(obj.mask_cyto);
                
                cmax = max(cell_img, [], 3, 'omitnan');
                ccyto2 = immultiply(cmax, obj.mask_cyto);
                ccyto2(ccyto2 == 0) = NaN;
                obj.cell_stats.cyto_intensity_total = sum(ccyto2, 'all', 'omitnan');
                obj.cell_stats.cyto_intensity_mean = mean(ccyto2, 'all', 'omitnan');
                obj.cell_stats.cyto_intensity_median = median(ccyto2, 'all', 'omitnan');
                obj.cell_stats.cyto_intensity_stdev = std(ccyto2, 0, 'all', 'omitnan');
            end

        end

        %%
        function pkg = packageForSave(obj)
            pkg = struct();
            pkg.version = 4;
            pkg.cell_number = obj.cell_number;
            pkg.cell_loc = obj.cell_loc;
            pkg.dim_z = obj.dim_z;
            pkg.mask_cell = obj.mask_cell;
            pkg.mask_nuc = obj.mask_nuc;
            pkg.spotcount_nuc = obj.spotcount_nuc;
            pkg.spotcount_total = obj.spotcount_total;
            pkg.signal_nuc = obj.signal_nuc;
            pkg.signal_cyto = obj.signal_cyto;
            pkg.signal_total = obj.signal_total;

            pkg.nucCount = obj.nucCount;
            pkg.nucNascentCount = obj.nucNascentCount;
            pkg.nucCloud = obj.nucCloud;
            pkg.nucNascentCloud = obj.nucNascentCloud;
            pkg.cytoCount = obj.cytoCount;
            pkg.cytoCloud = obj.cytoCloud;
            pkg.cell_stats = obj.cell_stats;

            pkg.nuc_ellip = obj.nuc_ellip;
            
            if ~isempty(obj.clouds)
                pkg.clouds = arrayfun(@(cloud) cloud.packageForSave(), obj.clouds);
            else
                pkg.clouds = [];
            end
            
            %Spots
            if isempty(obj.spotTable) & ~isempty(obj.spots)
                obj.spotTable = SingleCell.spotObjArr2Table(obj.spots);
                spotList = obj.spots;

                spotCount = size(obj.spotTable, 1);
                obj.spotZFits = cell(1, spotCount);
                for i = 1:spotCount
                    obj.spotZFits{i} = spotList(i).gfit_slices;
                end
            end

            pkg.spotTable = obj.spotTable;
            pkg.spotZFits = obj.spotZFits;
        end

    end
    
    methods(Static)
        
        %%
        function [fieldNames, fieldTypes] = getSpotTableFields()
            fieldNames = {'uid' 'xinit' 'yinit' 'zinit' 'dropout_thresh' ...
                'likely_dup' 'in_cloud' 'fit_volume' 'nascent_flag' ...
                'xfit' 'yfit' 'xgw' 'ygw' 'xFWHM' 'yFWHM' ...
                'expMInt' 'fitMInt' 'TotExpInt' 'TotFitInt' ...
                'r' 'rFit' 'back' 'xsem' 'ysem' 'zxfit' 'zyfit' ...
                'zint' 'zrel' 'zabs' 'zstd' 'zqFit' ...
                'nucRNA' 'cytoRNA' 'distRNA' 'szNUC' ...
                'xabsloc' 'yabsloc' 'normdist' 'nucR'};

            fieldTypes = {'uint32' 'uint16' 'uint16' 'uint16' 'uint16' ...
                'uint32' 'logical' 'single' 'logical' ...
                'single' 'single' 'single' 'single' 'single' 'single'...
                'uint32' 'single' 'uint32' 'single' ...
                'single' 'single' 'single' 'single' 'single' 'single' 'single' ...
                'single' 'single' 'single' 'single' 'single'...
                'logical' 'logical' 'single' 'single' ...
                'single' 'single' 'single' 'single'};
        end

        %%
        function mycell = readFromSavePackage(pkg)

            mycell = SingleCell;

            mycell.cell_number = pkg.cell_number;
            mycell.cell_loc = pkg.cell_loc;
            mycell.dim_z = pkg.dim_z;
            mycell.mask_cell = pkg.mask_cell;
            mycell.mask_nuc = pkg.mask_nuc;
            mycell.spotcount_nuc = pkg.spotcount_nuc;
            mycell.spotcount_total = pkg.spotcount_total;
            mycell.signal_nuc = pkg.signal_nuc;
            mycell.signal_cyto = pkg.signal_cyto;
            mycell.signal_total = pkg.signal_total;

            mycell.nucCount = pkg.nucCount;
            mycell.nucNascentCount = pkg.nucNascentCount;
            mycell.nucCloud = pkg.nucCloud;
            mycell.nucNascentCloud = pkg.nucNascentCloud;
            mycell.cytoCount = pkg.cytoCount;
            mycell.cytoCloud = pkg.cytoCloud;

            mycell.nuc_ellip = pkg.nuc_ellip;
            mycell.spotTable = pkg.spotTable;
            mycell.spotZFits = pkg.spotZFits;

            if pkg.version < 3
                %Gen UIDs and add dup/uid table fields
                spotcount = size(mycell.spotTable, 1);
                if spotcount > 0
                    mycell.spotTable{:, 'uid'} = uint32([1:1:spotcount]');
                    mycell.spotTable{:, 'likely_dup'} = uint32(zeros(spotcount, 1));
                end

%                 for i = 1:spotCount
%                     %TODO
%                     %Eventually, tag the zfits with UID as well. RN would
%                     %have to update anything accessing zfits though so do
%                     %later.
%                 end
            end
            if pkg.version >= 4
                mycell.cell_stats = pkg.cell_stats;
            end

            if ~isempty(pkg.clouds)
                mycell.clouds = arrayfun(@(cloudStruct) RNACloud.readFromSavePackage(cloudStruct), pkg.clouds);
            else
                mycell.clouds = [];
            end

        end

        %%
        function spotTable = spotObjArr2Table(spotArr)
            if isempty(spotArr)
                spotTable = table.empty();
                return;
            end

            [varNames, varTypes] = SingleCell.getSpotTableFields();

            alloc = size(spotArr, 2);
            table_size = [alloc size(varNames,2)];
            spotTable = table('Size', table_size, 'VariableTypes',varTypes, 'VariableNames',varNames);

            spotTable{:, 'in_cloud'} = [spotArr.in_cloud]';
            spotTable{:, 'fit_volume'} = single([spotArr.fit_volume]');
            spotTable{:, 'nascent_flag'} = [spotArr.nascent_flag]';

            fits = [spotArr.gauss_fit];
            spotTable{:, 'xinit'} = uint16([fits.xinit]');
            spotTable{:, 'yinit'} = uint16([fits.yinit]');
            spotTable{:, 'zinit'} = uint16([fits.zinit]');
            spotTable{:, 'dropout_thresh'} = uint16([spotArr.dropout_thresh]');

            spotTable{:, 'uid'} = [1:1:alloc]';
            spotTable{:, 'likely_dup'} = uint32(zeros(alloc, 1));

            spotTable{:, 'xfit'} = single([fits.xfit]');
            spotTable{:, 'yfit'} = single([fits.yfit]');
            spotTable{:, 'xgw'} = single([fits.xgw]');
            spotTable{:, 'ygw'} = single([fits.ygw]');
            spotTable{:, 'xFWHM'} = single([fits.xFWHM]');
            spotTable{:, 'yFWHM'} = single([fits.yFWHM]');
            spotTable{:, 'expMInt'} = uint32([fits.expMInt]');
            spotTable{:, 'fitMInt'} = single([fits.fitMInt]');
            spotTable{:, 'TotExpInt'} = uint32([fits.TotExpInt]');
            spotTable{:, 'TotFitInt'} = single([fits.TotFitInt]');
            spotTable{:, 'r'} = single([fits.r]');
            spotTable{:, 'rFit'} = single([fits.rFit]');
            spotTable{:, 'back'} = single([fits.back]');
            spotTable{:, 'xsem'} = single([fits.xsem]');
            spotTable{:, 'ysem'} = single([fits.ysem]');
            spotTable{:, 'zxfit'} = single([fits.zxfit]');
            spotTable{:, 'zyfit'} = single([fits.zyfit]');
            spotTable{:, 'zint'} = single([fits.zint]');
            spotTable{:, 'zrel'} = single([fits.zrel]');
            spotTable{:, 'zabs'} = single([fits.zabs]');
            spotTable{:, 'zstd'} = single([fits.zstd]');
            spotTable{:, 'zqFit'} = single([fits.zqFit]');
            spotTable{:, 'nucRNA'} = [fits.nucRNA]';
            spotTable{:, 'cytoRNA'} = [fits.cytoRNA]';
            spotTable{:, 'distRNA'} = single([fits.distRNA]');
            spotTable{:, 'szNUC'} = single([fits.szNUC]');
            spotTable{:, 'xabsloc'} = single([fits.xabsloc]');
            spotTable{:, 'yabsloc'} = single([fits.yabsloc]');
            spotTable{:, 'normdist'} = single([fits.normdist]');
            spotTable{:, 'nucR'} = single([fits.nucR]');
        end

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