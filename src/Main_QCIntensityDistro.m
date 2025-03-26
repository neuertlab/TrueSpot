%
%%
function Main_QCIntensityDistro(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2025.03.25.00';
VERSION_STRING = 'v1.1.2';

% ========================== Process args ==========================
arg_debug = true; %CONSTANT used for debugging arg parser.

opStruct = genOptionsStruct();

lastkey = [];
for i = 1:nargin
    argval = varargin{i};
    if ischar(argval) & startsWith(argval, "-")
        %Key
        if size(argval,2) >= 2
            lastkey = argval(2:end);
        else
            lastkey = [];
        end

        if strcmp(lastkey, "wspotmask")
            opStruct.outputSpotMasks = true;
            if arg_debug; fprintf("Spot Mask Output: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "nodpc")
            opStruct.nodpc = true;
            if arg_debug; fprintf("Dead Pixel Removal: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "noprobe")
            opStruct.forceNoProbe = true;
            if arg_debug; fprintf("Force No Probe: On\n"); end
            lastkey = [];
        end
        
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "input")
            opStruct.input_dir = argval;
            if arg_debug; fprintf("Input Directory Set: %s\n", opStruct.input_dir); end
        elseif strcmp(lastkey, "output")
            opStruct.output_dir = argval;
            if arg_debug; fprintf("Output Directory Set: %s\n", opStruct.output_dir); end
        elseif strcmp(lastkey, "zmin")
            opStruct.zmin = Force2Num(argval);
            if arg_debug; fprintf("Trim Z Min Set: %s\n", num2str(opStruct.zmin)); end
        elseif strcmp(lastkey, "zmax")
            opStruct.zmax = Force2Num(argval);
            if arg_debug; fprintf("Trim Z Max Set: %s\n", num2str(opStruct.zmax)); end
        elseif strcmp(lastkey, "localxy")
            opStruct.blockXY = Force2Num(argval);
            if arg_debug; fprintf("Local Block XY Size Set: %s\n", num2str(opStruct.blockXY)); end
        elseif strcmp(lastkey, "localz")
            opStruct.blockZ = Force2Num(argval);
            if arg_debug; fprintf("Local Block Z Size Set: %s\n", num2str(opStruct.blockZ)); end
        elseif strcmp(lastkey, "namecontains")
            opStruct.nameincl = argval;
            if arg_debug; fprintf("Run Names Must Include String: %s\n", opStruct.nameincl); end
        elseif strcmp(lastkey, "correctionmtx")
            opStruct.correctionmtx = argval;
            if arg_debug; fprintf("Correction matrix MAT path: %s\n", opStruct.correctionmtx); end
        elseif strcmp(lastkey, "workers")
            opStruct.workerCount = Force2Num(argval);
            if arg_debug; fprintf("Requested parallel workers: %d\n", opStruct.workerCount); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

%--- Check args (Fill in defaults based on inputs)

if isempty(opStruct.input_dir)
    fprintf('Please provide an input directory!\n');
    return;
end

if isempty(opStruct.output_dir)
    opStruct.output_dir = [opStruct.input_dir filesep 'IntensityDistro'];
end

fprintf('Main_QCIntensityDistro\n');
fprintf('Script Version: %s\n', BUILD_STRING);
fprintf('TrueSpot Version: %s\n', VERSION_STRING);
fprintf('Run Start: %s\n', datetime);
fprintf('Input: %s\n', opStruct.input_dir);
fprintf('Output: %s\n', opStruct.output_dir);

if ~isfolder(opStruct.output_dir); mkdir(opStruct.output_dir); end

% ========================== Recursive Scan ==========================

try
    runList = scanDirRec(opStruct.input_dir, opStruct);
    if ~isempty(runList)
        if opStruct.workerCount > 1
            processRuns_P(runList, opStruct);
        else
            processRuns_S(runList, opStruct);
        end
    else
        fprintf('No valid runs found! Exiting...\n');
    end

    
catch MEx
    fprintf('There was an error processing %s. Skipped remainder. See below for details.\n', opStruct.input_dir);
    %disp(MEx.message);
    disp(MEx.getReport('extended'));
end

end

% ========================== Additional Functions ==========================

function opStruct = genOptionsStruct()
    opStruct = struct();
    opStruct.input_dir = [];
    opStruct.output_dir = []; %Defaults to input dir/IntensityDistro. Place to put intensity stat profiles.
    opStruct.zmin = 0; %Min z to use for trimmed stats
    opStruct.zmax = 0; %Max z to use for trimmed stats
    
    opStruct.blockZ = 8;
    opStruct.blockXY = 256;
    opStruct.outputSpotMasks = false;
    opStruct.nodpc = false;
    opStruct.nameincl = [];
    opStruct.forceNoProbe = false;
    opStruct.correctionmtx = [];

    opStruct.workerCount = 1;
end

function runStruct = genRunCandidateStruct()
    runStruct = struct();
    runStruct.srPath = [];
    runStruct.qdPath = [];
end

function [fieldNames, fieldTypes] = getSpotProfileTableFields()
    fieldNames = {'x' 'y' 'z' 'cell' 'dropout_thresh' 'nascent_flag' 'nuc_flag' ...
        'xgw' 'ygw' 'xFWHM' 'yFWHM' ...
        'expMInt' 'fitMInt' 'TotExpInt' 'TotFitInt' 'mFiltInt' ... 
        'bkgLevelq' 'localBkgStats' 'noBkgSubStats'};
    
    fieldTypes = {'single' 'single' 'single' 'uint16' 'uint16' 'logical' 'logical' ...
        'single' 'single' 'single' 'single' ...
        'uint16' 'single' 'uint32' 'single' 'uint16' ...
        'single' 'cell' 'cell'};
end

function [fieldNames, fieldTypes] = getLocalRegionStatsFields()
    fieldNames = {'x' 'y' 'z' 'width' 'height' 'depth' ...
        'min' 'max' 'median' 'mad' 'mean' 'stdev' ...
        'histo_x' 'histo_y'};
    
    fieldTypes = {'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' ...
        'uint16' 'uint16' 'single' 'single' 'single' 'single'...
        'cell' 'cell'};
end

function [fieldNames, fieldTypes] = getCellTableFields()
    fieldNames = {'cellNo' 'centroid_x' 'centroid_y' 'box' 'spot_count' ...
        'local_img_bkg' 'cell_bkg' 'nuc_bkg' 'cyto_bkg'};
    
    fieldTypes = {'uint16' 'single' 'single' 'cell' 'int32' ...
        'cell' 'cell' 'cell' 'cell'};
end

function statsStruct = takeStats(data, mask, statsStruct, localXY, localZ, auto2D)
    if nargin < 4; localXY = 0; end
    if nargin < 5; localZ = 0; end
    if nargin < 6; auto2D = false; end

    if auto2D & (ndims(data) > 2)
        %Make life easy and just recurse...
        Z = size(data, 3);
        data_maxproj = max(data, [], 3, 'omitnan');
        mask_maxproj = max(mask, [], 3, 'omitnan');
        statsStruct.maxproj = takeStats(data_maxproj, mask_maxproj, struct(), localXY, localZ, false);
        statsStruct.slices = cell(1,Z);
        for z = 1:Z
            statsStruct.slices{z} = takeStats(data(:,:,z), mask(:,:,z), struct(), localXY, localZ, false);
        end
    end

%Apply mask if nonempty
    if ~isempty(mask) & (nnz(~mask) > 0)
        data = double(data);
        %fprintf('DEBUG -- Expected NaNs: %d\n', nnz(~mask));
        data(~mask) = NaN;
        %fprintf('DEBUG -- Added NaNs: %d\n', nnz(isnan(data)));
    end

    if nnz(isfinite(data)) < 1
        statsStruct.histo_y = [];
        statsStruct.histo_x = [];
        statsStruct.min = 0;
        statsStruct.max = 0;
        statsStruct.median = NaN;
        statsStruct.mad = NaN;
        statsStruct.mean = NaN;
        statsStruct.stdev = NaN;
        return;
    end

%Histo
    data_all = uint16(data(isfinite(data)));
    statsStruct.min = min(data_all, [], 'all', 'omitnan');
    statsStruct.max = max(data_all, [], 'all', 'omitnan');

    edges = [0:1:statsStruct.max];
    [bins, edges] = histcounts(data_all, edges);
    statsStruct.histo_y = uint32(bins);
    statsStruct.histo_x = uint16(edges(1:size(bins,2)));
    clear edges bins

    data_all = double(data_all);

    statsStruct.median = median(data_all, 'all', 'omitnan');
    statsStruct.mad = mad(data_all, 1, 'all');
    statsStruct.mean = mean(data_all, 'all', 'omitnan');
    statsStruct.stdev = std(data_all, 0, 'all', 'omitnan');
    clear data_all

%Local blocks, if applicable
    if (localXY > 0) & (localZ > 0)
        X = size(data, 2);
        Y = size(data, 1);
        Z = size(data, 3);
        zBlocks = ceil(Z ./ localZ);
        xBlocks = ceil(X ./ localXY);
        yBlocks = ceil(Y ./ localXY);
        totalBlocks = zBlocks * xBlocks * yBlocks;

        [fieldNames, fieldTypes] = getLocalRegionStatsFields();
        table_size = [totalBlocks size(fieldNames,2)];
        statsStruct.localBlocks = table('Size', table_size, 'VariableTypes',fieldTypes, 'VariableNames',fieldNames);

        i = 1;
        z0 = 1;
        for zb = 1:zBlocks
            z1 = min((z0 + localZ) - 1, Z);
            x0 = 1;
            for xb = 1:xBlocks
                x1 = min((x0 + localXY) - 1, X);
                y0 = 1;
                for yb = 1:yBlocks
                    y1 = min((y0 + localXY) - 1, Y);

                    %Do actual block.
                    statsStruct.localBlocks{i, 'x'} = x0;
                    statsStruct.localBlocks{i, 'y'} = y0;
                    statsStruct.localBlocks{i, 'z'} = z0;
                    statsStruct.localBlocks{i, 'width'} = x1 - x0 + 1;
                    statsStruct.localBlocks{i, 'height'} = y1 - y0 + 1;
                    statsStruct.localBlocks{i, 'depth'} = z1 - z0 + 1;

                    %Don't need mask again because it's already been
                    %applied.
                    localSStruct = ...
                        takeStats(data(y0:y1,x0:x1,z0:z1), [], struct(), 0, 0);

                    statsStruct.localBlocks{i, 'min'} = localSStruct.min;
                    statsStruct.localBlocks{i, 'max'} = localSStruct.max;
                    statsStruct.localBlocks{i, 'median'} = single(localSStruct.median);
                    statsStruct.localBlocks{i, 'mad'} = single(localSStruct.mad);
                    statsStruct.localBlocks{i, 'mean'} = single(localSStruct.mean);
                    statsStruct.localBlocks{i, 'stdev'} = single(localSStruct.stdev);

                    ccx = cell(1,1);
                    ccx{1} = localSStruct.histo_x;
                    ccy = cell(1,1);
                    ccy{1} = localSStruct.histo_y;
                    statsStruct.localBlocks{i, 'histo_x'} = ccx;
                    statsStruct.localBlocks{i, 'histo_y'} = ccy;

                    y0 = y1 + 1;
                    i = i + 1;
                end
                x0 = x1 + 1;
            end
            z0 = z1+1;
        end
    end

end

function [spotMask, indivSpotMasks] = genSpotMaskQuant(quant_results, dims, gaussThresh)
    if(nargin < 3); gaussThresh = 0.125; end

    spotMask = false(dims);
    X = size(spotMask, 2);
    Y = size(spotMask, 1);
    Z = size(spotMask, 3);
    myCells = quant_results.cell_rna_data;
    cellCount = size(myCells, 2);
    indivSpotMasks = cell(1, cellCount);

    for c = 1:cellCount
        myCell = myCells(c);
        if ~isempty(myCell.spotTable)
            cellSpots = size(myCell.spotTable, 1);
            cellSpotMasks = cell(1, cellSpots);

            cell_x = myCell.cell_loc.left;
            cell_y = myCell.cell_loc.top;
            cell_z = myCell.cell_loc.z_bottom;
            for s = 1:cellSpots
                %First, determine xyz widths of box
                %Just use double FWHM
                xFWHM = myCell.spotTable{s, 'xFWHM'};
                yFWHM = myCell.spotTable{s, 'yFWHM'};
                xyrad = max(xFWHM, yFWHM);
                xyrad = ceil(xyrad * 2);

                %Expected center of output image
                x_snap = myCell.spotTable{s, 'xinit'} + cell_x - 1;
                y_snap = myCell.spotTable{s, 'yinit'} + cell_y - 1;
                z_snap = myCell.spotTable{s, 'zinit'} + cell_z - 1;
                sim_spot = RNASpot.generateSimSpotFromFit_Table(myCell.spotTable, s, myCell.spotZFits{s}, xyrad);
                mOvrl = max(sim_spot, [], 'all', 'omitnan');
                sim_mask = (sim_spot >= (gaussThresh * mOvrl));

                [x_min, x_max, xtrim_lo, xtrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(x_snap, X, xyrad);
                [y_min, y_max, ytrim_lo, ytrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(y_snap, Y, xyrad);
                [z_min, z_max, ztrim_lo, ztrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(z_snap, Z, 2);

                y0 = 1 + ytrim_lo;
                y1 = size(sim_mask, 1) - ytrim_hi;
                x0 = 1 + xtrim_lo;
                x1 = size(sim_mask, 2) - xtrim_hi;
                z0 = 1 + ztrim_lo;
                z1 = size(sim_mask, 3) - ztrim_hi;

                spotMask(y_min:y_max, x_min:x_max, z_min:z_max) =...
                    or(sim_mask(y0:y1, x0:x1, z0:z1), ...
                    spotMask(y_min:y_max, x_min:x_max, z_min:z_max));

                %Also save mask and top right corner coords
                mySpotMask = struct();
                mySpotMask.corner_x = x_min;
                mySpotMask.corner_y = y_min;
                mySpotMask.corner_z = z_min;
                mySpotMask.mask = sim_mask(y0:y1, x0:x1, z0:z1);

                cellSpotMasks{s} = mySpotMask;
            end
            indivSpotMasks{c} = cellSpotMasks;
        end
    end

end

function [spotMask, call_table] = genSpotMaskCall(spotsrun, gaussrad, gaussThresh)
    if(nargin < 3); gaussThresh = 0.125; end

    [~, call_table] = spotsrun.loadCallTable();
    th = spotsrun.intensity_threshold;
    call_table = call_table((call_table{:, 'dropout_thresh'} >= th), :);

    Z = spotsrun.dims.idims_sample.z;
    Y = spotsrun.dims.idims_sample.y;
    X = spotsrun.dims.idims_sample.x;

    spotMask = false(Y,X,Z);

    gaussrad_z = 2; %Fixed;
    xydim = (gaussrad * 2) + 1;
    zdim = (gaussrad_z * 2) + 1;
    xy_mid = round(xydim/2);
    z_mid = round(zdim/2);
    fwhm_xy = gaussrad/2;
    fwhm_z = gaussrad_z/2;
    xy_w = fwhm_xy ./ (2 .* sqrt(2.*log(2)));
    z_w = fwhm_z ./ (2 .* sqrt(2.*log(2)));

    gauss_spot = RNAUtils.generateGaussian3D(xydim, xydim, zdim, xy_mid-1, xy_mid-1, z_mid-1, xy_w, xy_w, z_w, 1);
    gauss_spot_mask = (gauss_spot >= gaussThresh);

    spotCount = size(call_table,1);
    for s = 1:spotCount
        [x_min, x_max, xtrim_lo, xtrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(call_table{s,'isnap_x'}, X, gaussrad);
        [y_min, y_max, ytrim_lo, ytrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(call_table{s,'isnap_y'}, Y, gaussrad);
        [z_min, z_max, ztrim_lo, ztrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(call_table{s,'isnap_z'}, Z, gaussrad_z);

%         [x_min, x_max, xtrim_lo, xtrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(1, X, gaussrad);
%         [y_min, y_max, ytrim_lo, ytrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(1, Y, gaussrad);
%         [z_min, z_max, ztrim_lo, ztrim_hi, ~] = RNAUtils.getDimSpotIsolationParams(1, Z, gaussrad_z);

        spotMask(y_min:y_max, x_min:x_max, z_min:z_max) = ...
            or(spotMask(y_min:y_max, x_min:x_max, z_min:z_max), ...
            gauss_spot_mask((1+ytrim_lo):(xydim-ytrim_hi), ...
                            (1+xtrim_lo):(xydim-xtrim_hi), ...
                            (1+ztrim_lo):(zdim-ztrim_hi)));
    end

end

function spotProfile = getSpotProfileQuant(spotsrun, quant_results, imageData, spotMask, indivSpotMasks)
    myCells = quant_results.cell_rna_data;
    cellCount = size(myCells, 2);

    spotCount = 0;
    for c = 1:cellCount
        myCell = myCells(c);
        if ~isempty(myCell.spotTable)
            spotCount = spotCount + size(myCell.spotTable, 1);
        end
    end

    X = size(imageData, 2);
    Y = size(imageData, 1);
    Z = size(imageData, 3);

    [fieldNames, fieldTypes] = getSpotProfileTableFields();
    table_size = [spotCount size(fieldNames,2)];
    spTable = table('Size', table_size, 'VariableTypes',fieldTypes, 'VariableNames',fieldNames);
    %Merge together all spot tables
    spotPos = 1;
    for c = 1:cellCount
        myCell = myCells(c);
        cellSpotMasks = indivSpotMasks{c};
        if ~isempty(myCell.spotTable)
            cspotCount = size(myCell.spotTable, 1);
            lastSpot = spotPos + cspotCount - 1;
            spTable{spotPos:lastSpot, 'cell'} = c;
            spTable{spotPos:lastSpot, 'x'} = myCell.spotTable{:, 'xabsloc'};
            spTable{spotPos:lastSpot, 'y'} = myCell.spotTable{:, 'yabsloc'};
            spTable{spotPos:lastSpot, 'z'} = myCell.spotTable{:, 'zabs'};
            spTable{spotPos:lastSpot, 'dropout_thresh'} = myCell.spotTable{:, 'dropout_thresh'};
            spTable{spotPos:lastSpot, 'nascent_flag'} = myCell.spotTable{:, 'nascent_flag'};
            spTable{spotPos:lastSpot, 'nuc_flag'} = myCell.spotTable{:, 'nucRNA'};
            spTable{spotPos:lastSpot, 'xgw'} = myCell.spotTable{:, 'xgw'};
            spTable{spotPos:lastSpot, 'ygw'} = myCell.spotTable{:, 'ygw'};
            spTable{spotPos:lastSpot, 'xFWHM'} = myCell.spotTable{:, 'xFWHM'};
            spTable{spotPos:lastSpot, 'yFWHM'} = myCell.spotTable{:, 'yFWHM'};
            spTable{spotPos:lastSpot, 'expMInt'} = myCell.spotTable{:, 'expMInt'};
            spTable{spotPos:lastSpot, 'fitMInt'} = myCell.spotTable{:, 'fitMInt'};
            spTable{spotPos:lastSpot, 'TotExpInt'} = myCell.spotTable{:, 'TotExpInt'};
            spTable{spotPos:lastSpot, 'TotFitInt'} = myCell.spotTable{:, 'TotFitInt'};
            spTable{spotPos:lastSpot, 'bkgLevelq'} = myCell.spotTable{:, 'back'};

            %Clean up any NaN or negative (bugged) z positions...
            zbad = isnan(spTable{spotPos:lastSpot, 'z'});
            zbad = or(zbad, spTable{spotPos:lastSpot, 'z'} < 1);
            if nnz(zbad) > 0
                zint = myCell.spotTable{:, 'zinit'} + myCell.cell_loc.z_bottom - 1;
                rows = find(zbad) + spotPos - 1;
                spTable{rows, 'z'} = zint(zbad);
            end

            x0 = int32(max(floor(spTable{spotPos:lastSpot, 'x'} - spTable{spotPos:lastSpot, 'xFWHM'} - 5), 1));
            x1 = int32(min(ceil(spTable{spotPos:lastSpot, 'x'} + spTable{spotPos:lastSpot, 'xFWHM'} + 5), X));
            y0 = int32(max(floor(spTable{spotPos:lastSpot, 'y'} - spTable{spotPos:lastSpot, 'yFWHM'} - 5), 1));
            y1 = int32(min(ceil(spTable{spotPos:lastSpot, 'y'} + spTable{spotPos:lastSpot, 'yFWHM'} + 5), Y));
            z0 = int32(max(floor(spTable{spotPos:lastSpot, 'z'} - 3), 1));
            z1 = int32(min(ceil(spTable{spotPos:lastSpot, 'z'} + 3), Z));

            %Local background for each spot.
            for s = 1:cspotCount
                idat = imageData(y0(s):y1(s), x0(s):x1(s), z0(s):z1(s));
                lmask = spotMask(y0(s):y1(s), x0(s):x1(s), z0(s):z1(s));
                cll = cell(1,1);
                cll{1} = takeStats(idat, ~lmask, struct(), 0, 0);
                spTable{spotPos + s - 1, 'localBkgStats'} = cll;

                %Also get exp maximums without bkg subtraction
                mySpotMask = cellSpotMasks{s};
                xX = size(mySpotMask.mask, 2);
                yY = size(mySpotMask.mask, 1);
                zZ = size(mySpotMask.mask, 3);

                xx0 = mySpotMask.corner_x;
                yy0 = mySpotMask.corner_y;
                zz0 = mySpotMask.corner_z;
                xx1 = xx0 + xX - 1;
                yy1 = yy0 + yY - 1;
                zz1 = zz0 + zZ - 1;

                idat = imageData(yy0:yy1, xx0:xx1, zz0:zz1);

                cll = cell(1,1);
                cll{1} = takeStats(idat, mySpotMask.mask, struct(), 0, 0);
                spTable{spotPos + s - 1, 'noBkgSubStats'} = cll;
            end

            spotPos = lastSpot + 1;
        end
    end

    %Get filtered values
    fImg = RNA_Threshold_SpotDetector.run_spot_detection_pre(imageData, [spotsrun.paths.out_dir filesep], true, spotsrun.options.dtune_gaussrad, false);
    fImg = uint16(fImg);
    xr = int32(min(max(round(spTable{:, 'x'}),1), X));
    yr = int32(min(max(round(spTable{:, 'y'}),1), Y));
    zr = int32(min(max(round(spTable{:, 'z'}),1), Z));
    coords1 = sub2ind([Y X Z], yr, xr, zr);
    spTable{:, 'mFiltInt'} = fImg(coords1);
    clear fImg

    %Gen output struct and calculate overall stats of interest
    spotProfile = struct();
    spotProfile.statsExpMax = takeStats(spTable{:, 'expMInt'}, [], struct(), 0, 0);
    spotProfile.statsFitMax = takeStats(spTable{:, 'fitMInt'}, [], struct(), 0, 0);
    spotProfile.statsExpTotal = takeStats(spTable{:, 'TotExpInt'}, [], struct(), 0, 0);
    spotProfile.statsFitTotal = takeStats(spTable{:, 'TotFitInt'}, [], struct(), 0, 0);

    spotProfile.spotData = spTable;
end

function spotProfile = getSpotProfileCall(spotsrun, imageData, call_table, gaussrad)
    spotProfile = struct();

    if isempty(call_table)
        [~, call_table] = spotsrun.loadCallTable();
        if isempty(call_table); return; end
    end
   
    ithresh = spotsrun.intensity_threshold;
    iokay = call_table{:, 'dropout_thresh'} >= ithresh;
    call_table = call_table(iokay, :);

    spotCount = size(call_table, 1);
    [fieldNames, fieldTypes] = getSpotProfileTableFields();
    table_size = [spotCount size(fieldNames,2)];
    spTable = table('Size', table_size, 'VariableTypes',fieldTypes, 'VariableNames',fieldNames);

    spTable{:, 'cell'} = 0;
    spTable{:, 'x'} = single(call_table{:, 'isnap_x'});
    spTable{:, 'y'} = single(call_table{:, 'isnap_y'});
    spTable{:, 'z'} = single(call_table{:, 'isnap_z'});
    spTable{:, 'nascent_flag'} = false;
    spTable{:, 'nuc_flag'} = false;

    spTable{:, 'expMInt'} = single(call_table{:, 'intensity'});
    spTable{:, 'fitMInt'} = spTable{:, 'expMInt'};

    %Use gaussrad for width...
    spTable{:, 'xFWHM'} = single(gaussrad - 2);
    spTable{:, 'yFWHM'} = single(gaussrad - 2);

    spTable{:, 'xgw'} =  spTable{:, 'xFWHM'}./(2*sqrt(2*log(2)));
    spTable{:, 'ygw'} =  spTable{:, 'yFWHM'}./(2*sqrt(2*log(2)));

    %Not immediately obvious to me how to get around looping...
    xyrad = 2 .* (gaussrad + 2);
    zrad = 4;
    zw = 1/sqrt(2*log(2));
    for s = 1:spotCount
        x = spTable{s, 'x'};
        y = spTable{s, 'y'};
        z = spTable{s, 'z'};
        spot_data = RNAUtils.isolateSpotData(imageData, x, y, z, xyrad, zrad);
        [spot_mask, gauss_sim] = RNASpot.genSpotMask(size(spot_data), spTable{s, 'xgw'}, spTable{s, 'ygw'}, zw);
        spTable{s, 'TotExpInt'} = sum(spot_data(spot_mask), 'all');
        spTable{s, 'TotFitInt'} = sum(gauss_sim(spot_mask), 'all');
    end

    %Finish up struct...
    spotProfile.statsExpMax = takeStats(spTable{:, 'expMInt'}, [], struct(), 0, 0);
    spotProfile.statsFitMax = takeStats(spTable{:, 'fitMInt'}, [], struct(), 0, 0);
    spotProfile.statsExpTotal = takeStats(spTable{:, 'TotExpInt'}, [], struct(), 0, 0);
    spotProfile.statsFitTotal = takeStats(spTable{:, 'TotFitInt'}, [], struct(), 0, 0);

    spotProfile.spotData = spTable;
end

function runList = scanDirRec(dirPath, opStruct)
    fprintf('Scanning %s ...\n', dirPath);
    dirContents = dir(dirPath);
    childCount = size(dirContents, 1);

    srPath = [];
    qdPath = [];

    runsFound = 0;
    rawList = cell(1, childCount);
    runList = [];
    for i = 1:childCount
        fname = dirContents(i,1).name;
        isdir = dirContents(i,1).isdir;
        if isdir
            if ~strcmp(fname, '.') & ~strcmp(fname, '..')
                subdirPath = [dirPath filesep fname];
                try
                    subRuns = scanDirRec(subdirPath, opStruct);
                    if ~isempty(subRuns)
                        rawList{i} = subRuns;
                        runsFound = runsFound + size(subRuns, 2);
                    end
                catch MEx
                    fprintf('There was an error processing %s. Skipped remainder. See below for details.\n', subdirPath);
                    %disp(MEx.message);
                    disp(MEx.getReport('extended'));
                end
            end
        else
            if ~isempty(opStruct.nameincl)
                if contains(fname, opStruct.nameincl)
                    if endsWith(fname, '_rnaspotsrun.mat')
                        srPath = [dirPath filesep fname];
                    elseif endsWith(fname, '_quantData.mat')
                        qdPath = [dirPath filesep fname];
                    end
                end
            else
                if endsWith(fname, '_rnaspotsrun.mat')
                    srPath = [dirPath filesep fname];
                elseif endsWith(fname, '_quantData.mat')
                    qdPath = [dirPath filesep fname];
                end
            end
        end
    end

    if ~isempty(srPath)
        myRun = genRunCandidateStruct();
        myRun.srPath = srPath;
        myRun.qdPath = qdPath;
        runsFound = runsFound + 1;
    else
        myRun = [];
    end

    %Condense raw list into list of valid run structs
    if runsFound > 0
        runList = cell(1,runsFound);
        index = 1;
        if ~isempty(myRun) 
            runList{1} = myRun;
            index = 2;
        end

        for i = 1:childCount
            rawRes = rawList{i};
            if ~isempty(rawRes)
                subList = rawRes;
                subCount = size(subList, 2);
                for j = 1:subCount
                    runList{index} = subList{j};
                    index = index + 1;
                end
            end
        end
    end

end

function processRuns_S(runList, opStruct)
    if isempty(runList); return; end
    runCount = size(runList, 2);
    for r = 1:runCount
        runStruct = runList{r};
        processRun(runStruct, opStruct);
    end
end

function initParallel(threads, workers, workdir)
    mkdir(workdir);
    local_cluster = parcluster();
    local_cluster.NumThreads = threads;
    local_cluster.NumWorkers = workers;
    local_cluster.JobStorageLocation = workdir;
    parpool(local_cluster);
end

function processRuns_P(runList, opStruct)
    if isempty(runList); return; end

    opStruct_clean.input_dir = opStruct.input_dir;
    opStruct_clean.output_dir = opStruct.output_dir;
    opStruct_clean.zmin = opStruct.zmin;
    opStruct_clean.zmax = opStruct.zmax;
    opStruct_clean.blockZ = opStruct.blockZ;
    opStruct_clean.blockXY = opStruct.blockXY;
    opStruct_clean.outputSpotMasks = opStruct.outputSpotMasks;
    opStruct_clean.nodpc = opStruct.nodpc;
    opStruct_clean.nameincl = opStruct.nameincl;
    opStruct_clean.forceNoProbe = opStruct.forceNoProbe;
    opStruct_clean.correctionmtx = opStruct.correctionmtx;

    rlistDup = runList;
    
    fprintf("Now testing thresholds using %d workers...\n", opStruct.workerCount);
    tic;
    
    %parpool(common_ctx.threads);
    pardir_sfx = [];
    if ~isempty(opStruct.nameincl)
        pardir_sfx = ['_' replace(opStruct.nameincl, ' ', '_')];
    end
    pardir = [opStruct.output_dir filesep 'parallel' pardir_sfx];
    initParallel(1, opStruct.workerCount, pardir);

    runCount = size(runList, 2);
    parfor r = 1:runCount
        runStruct = rlistDup{r};
        processRun(runStruct, opStruct_clean);
    end
    delete(gcp('nocreate'));
    rmres = rmdir(pardir, 's');
    toc;

end

function processRun(runStruct, opStruct)
    try
        intensityStats = struct();
        spotsrun = RNASpotsRun.loadFrom(runStruct.srPath, false);
        fprintf('\tRun found: %s\n', spotsrun.img_name);
    
        %Update threshold stats
        tres = [];
        if ~isempty(spotsrun.threshold_results)
            fprintf('\t[%s] Updating threshold pool stats...\n', intensityStats.img_name);
            spotsrun.threshold_results = RNAThreshold.scoreThresholdSuggestions(spotsrun.threshold_results);
            tres = spotsrun.threshold_results;
            spotsrun.saveMeTo(runStruct.srPath);
        end
    
        intensityStats.img_name = spotsrun.img_name;
        intensityStats.channel = spotsrun.channels.rna_ch;
        intensityStats.target = spotsrun.meta.type_target;
        intensityStats.probe = spotsrun.meta.type_probe;
        intensityStats.t_min = spotsrun.options.t_min;
        intensityStats.t_max = spotsrun.options.t_max;
        intensityStats.threshold = spotsrun.intensity_threshold;
        intensityStats.gaussrad = spotsrun.options.dtune_gaussrad;
        intensityStats.noprobe = false;
        intensityStats.trim_z_min = opStruct.zmin;
        intensityStats.trim_z_max = opStruct.zmax;
        intensityStats.threshold_results = tres;
        intensityStats.nodpc = opStruct.nodpc;
    
        if opStruct.forceNoProbe
            intensityStats.noprobe = true;
        else
            if isfield(spotsrun.meta, 'noProbe_flag')
                intensityStats.noprobe = spotsrun.meta.noProbe_flag;
            end
            if isfield(spotsrun.options, 'noProbe_flag')
                %From a bug in the Spots main. Overrides.
                intensityStats.noprobe = spotsrun.options.noProbe_flag;
            end
        end
    
        %Load source image
        tifPath = spotsrun.paths.img_path;
        if isempty(tifPath) | ~isfile(tifPath)
            fprintf('\t[%s] ERROR: Spotsrun was found in directory, but TIF file could not be linked! Cannot calculate image stats...\n', intensityStats.img_name);
            return;
        end
    
        [channels, ~] = LoadTif(tifPath, spotsrun.channels.total_ch, [spotsrun.channels.rna_ch], 1);
        sampleCh = channels{spotsrun.channels.rna_ch, 1};
        sampleCh = uint16(sampleCh); %Reduce memory size
        clear channels;
    
        if ~isempty(opStruct.correctionmtx)
            fprintf('\t[%s] Loading intensity correction matrix from: %s\n', intensityStats.img_name, opStruct.correctionmtx);
            load(opStruct.correctionmtx, 'BaseLineControl_ImageCorrectionFactor');
    
            sampleCh = double(sampleCh);
            sampleCh = sampleCh ./ BaseLineControl_ImageCorrectionFactor;
            sampleCh = uint16(round(sampleCh));
    
            clear BaseLineControl_ImageCorrectionFactor
        end
    
        if ~opStruct.nodpc
            %Clean out dead pixels
            dead_pix_path = [spotsrun.getFullOutStem() '_deadpix_idistro.mat'];
            fprintf('\t[%s] Cleaning dead pixels. Dead pixel annotation saving to: %s\n', intensityStats.img_name, dead_pix_path);
            RNA_Threshold_Common.saveDeadPixels(sampleCh, dead_pix_path, true);
            sampleCh = RNA_Threshold_Common.cleanDeadPixels(sampleCh, dead_pix_path, true);
        end
    
        zmin = opStruct.zmin;
        zmax = opStruct.zmax;
    
        Z = size(sampleCh, 3);
        Y = size(sampleCh, 1);
        X = size(sampleCh, 2);
        if zmin < 1; zmin = 1; end
        if zmax < 1; zmax = Z; end
        doTrimmed = false;
        if zmin > 1; doTrimmed = true; end
        if zmax < Z; doTrimmed = true; end
    
        intensityStats.statsRawFull = takeStats(sampleCh, [], struct(), 0, 0);
    
        if(doTrimmed)
            intensityStats.statsRawZTrim = takeStats(sampleCh(:,:,zmin:zmax), [], struct(), 0, 0);
        end
    
        %Load cellseg and quant
        cspath = spotsrun.paths.cellseg_path;
        cellMask = [];
        nucMask = [];
        if ~isempty(cspath)
            if isfile(cspath)
                cellMask = CellSeg.openCellMask(cspath);
                nucMask = CellSeg.openNucMask(cspath);
    
                %Convert cellMask to 3D...
                if ndims(cellMask) < 3
                    cellMask = repmat(cellMask, [1, 1, Z]);
                end
                cellMaskBool = (cellMask ~= 0);
                cellCount = max(cellMask, [], 'all', 'omitnan');
            end
        end
    
        quant_results = [];
        if ~intensityStats.noprobe
            if ~isempty(runStruct.qdPath)
                load(runStruct.qdPath, 'quant_results');
                if ~isfield(quant_results, 'cell_rna_data')
                    quant_results = RNAQuant.readResultsSavePackage(quant_results);
                end
            end
    
            if ~isempty(quant_results)
                [spotMask, indivSpotMasks] = genSpotMaskQuant(quant_results, size(sampleCh));
            else
                [spotMask, call_table] = genSpotMaskCall(spotsrun, intensityStats.gaussrad);
                call_table = RNACoords.applyCellSegMask(call_table, cellMask);
            end
        else
            spotMask = false(Y,X,Z);
        end
    
        fprintf('\t[%s] > Working on image background...\n', intensityStats.img_name);
        currentMask = and(~cellMaskBool, ~spotMask);
        intensityStats.statsRawFull.ibkg = takeStats(sampleCh, currentMask, struct(), opStruct.blockXY, opStruct.blockZ);
        if(doTrimmed)
            intensityStats.statsRawZTrim.ibkg = takeStats(sampleCh(:,:,zmin:zmax), currentMask(:,:,zmin:zmax), struct(), 0, 0);
        end
    
        fprintf('\t[%s] > Working on cell background...\n', intensityStats.img_name);
        currentMask = and(cellMaskBool, ~spotMask);
        intensityStats.statsRawFull.cbkg = takeStats(sampleCh, currentMask, struct(), opStruct.blockXY, opStruct.blockZ);
        if(doTrimmed)
            intensityStats.statsRawZTrim.cbkg = takeStats(sampleCh(:,:,zmin:zmax), currentMask(:,:,zmin:zmax), struct(), 0, 0);
        end
    
        if ~intensityStats.noprobe
            fprintf('\t[%s] > Working on overall spots...\n', intensityStats.img_name);
            intensityStats.statsRawFull.spotsTotal = takeStats(sampleCh, spotMask, struct(), opStruct.blockXY, opStruct.blockZ);
            if(doTrimmed)
                intensityStats.statsRawZTrim.spotsTotal = takeStats(sampleCh(:,:,zmin:zmax), spotMask(:,:,zmin:zmax), struct(), 0, 0);
            end
    
            fprintf('\t[%s] > Working on spot profiles...\n', intensityStats.img_name);
            if ~isempty(quant_results)
                intensityStats.spotProfile = getSpotProfileQuant(spotsrun, quant_results, sampleCh, spotMask, indivSpotMasks);
            else
                intensityStats.spotProfile = getSpotProfileCall(spotsrun, sampleCh, call_table, intensityStats.gaussrad);
            end
        end
    
        fprintf('\t[%s] > Working on individual cells...\n', intensityStats.img_name);
        [fieldNames, fieldTypes] = getCellTableFields();
        table_size = [cellCount size(fieldNames,2)];
        cellTable = table('Size', table_size, 'VariableTypes',fieldTypes, 'VariableNames',fieldNames);
    
        for c = 1:cellCount
            cellTable{c, 'cellNo'} = uint16(c);
            oneCellMask = (cellMask == c);
            rp3 = regionprops3(oneCellMask, 'BoundingBox', 'Centroid');
            cellTable{c, 'centroid_x'} = single(rp3.Centroid(1));
            cellTable{c, 'centroid_y'} = single(rp3.Centroid(2));
    
            if ~isempty(quant_results)
                cobj = quant_results.cell_rna_data(c);
                cll = cell(1,1);
                cll{1} = cobj.cell_loc;
                cellTable{c, 'box'} = cll;
                cellTable{c, 'spot_count'} = size(cobj.spotTable, 1);
    
                x0 = max(cobj.cell_loc.left - 5,1);
                x1 = min(cobj.cell_loc.right + 5,X);
                y0 = max(cobj.cell_loc.top - 5, 1);
                y1 = min(cobj.cell_loc.bottom + 5,Y);
                z0 = max(cobj.cell_loc.z_bottom - 5,1);
                z1 = min(cobj.cell_loc.z_top + 5,Z);
            else
                %Cell bounds from regionprops
                x0 = max(round(rp3.BoundingBox(1)),1);
                x1 = min(round(x0 + rp3.BoundingBox(4)),X);
                y0 = max(round(rp3.BoundingBox(2)),1);
                y1 = min(round(y0 + rp3.BoundingBox(5)),Y);
                z0 = max(round(rp3.BoundingBox(3)),1);
                z1 = min(round(z0 + rp3.BoundingBox(6)),Z);
    
                cll = cell(1,1);
                cll{1} = SingleCell.generateRecPrismStruct(x0, x1, y0, y1, z0, z1);
                cellTable{c, 'box'} = cll;
    
                x0 = max(x0 - 5,1);
                x1 = min(x1 + 5,X);
                y0 = max(y0 - 5, 1);
                y1 = min(y1 + 5,Y);
                z0 = max(z0 - 5,1);
                z1 = min(z1 + 5,Z);
    
                %spot count from call table
                if ~intensityStats.noprobe
                    cellTable{c, 'spot_count'} = nnz(call_table{:,'cell'} == c);
                end
            end
        
            cellData = sampleCh(y0:y1,x0:x1,z0:z1);
            oneCellMask = oneCellMask(y0:y1,x0:x1,z0:z1);
            cellLocMask = ~oneCellMask;

            cll = cell(1,1);
            cll{1} = takeStats(cellData, cellLocMask, struct(), 0, 0, true);
            cellTable{c, 'local_img_bkg'} = cll;
            clear cellLocMask
    
            cellSpotMask = spotMask(y0:y1,x0:x1,z0:z1);
            cbkgMask = and(oneCellMask, ~cellSpotMask);
            cll{1} = takeStats(cellData, cbkgMask, struct(), 0, 0, true);
            cellTable{c, 'cell_bkg'} = cll;
    
            cellNucMask = nucMask(y0:y1,x0:x1,z0:z1);
            nucBkgMask = and(cbkgMask, cellNucMask);
            cll{1} = takeStats(cellData, nucBkgMask, struct(), 0, 0, true);
            cellTable{c, 'nuc_bkg'} = cll;

            cytoBkgMask = and(cbkgMask, ~cellNucMask);
            cll{1} = takeStats(cellData, cytoBkgMask, struct(), 0, 0, true);
            cellTable{c, 'cyto_bkg'} = cll;
    
            clear cellData oneCellMask cellNucMask
            clear cbkgMask cytoBkgMask nucBkgMask
        end
    
        intensityStats.cellProfile = cellTable;
    
        %Save
        savepath = [opStruct.output_dir, filesep, spotsrun.img_name '_istats.mat'];
        save(savepath, 'intensityStats');
    
        if opStruct.outputSpotMasks & ~intensityStats.noprobe
            savepath = [opStruct.output_dir, filesep, spotsrun.img_name '_spotmask.tif'];
            tifop = struct('overwrite', true);
            spotMask = uint8(spotMask);
            saveastiff(spotMask, savepath, tifop);
        end
    
    catch MEx
        fprintf('There was an error processing %s. Skipped remainder. See below for details.\n', runStruct.srPath);
        %disp(MEx.message);
        disp(MEx.getReport('extended'));
    end
end

