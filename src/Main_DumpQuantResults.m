function Main_DumpQuantResults(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2024.07.05.00';
VERSION_STRING = 'v1.1.0';

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.

input_dir = []; %Give a dir for a batch and it will look recursively through it for TS results.
output_path = [];

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
        
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "input")
            input_dir = argval;
            if arg_debug; fprintf("Input Directory Set: %s\n", input_dir); end
        elseif strcmp(lastkey, "output")
            output_path = argval;
            if arg_debug; fprintf("Output Path Set: %s\n", output_path); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

%--- Check args (Fill in defaults based on inputs)

if isempty(input_dir)
    fprintf('Please provide an input directory!\n');
    return;
end

if isempty(output_path)
    output_path = [input_dir filesep 'cellCounts.tsv'];
end

fprintf('Main_DumpQuantResults\n');
fprintf('Script Version: %s\n', BUILD_STRING);
fprintf('TrueSpot Version: %s\n', VERSION_STRING);
fprintf('Run Time: %s\n', datetime);
fprintf('Input: %s\n', input_dir);
fprintf('Output: %s\n', output_path);

tableHandle = openOutTable(output_path);
doDir(input_dir, tableHandle);

fclose(tableHandle);
end

% ========================== Additional Functions ==========================

function tableHandle = openOutTable(tablePath)
    tableHandle = fopen(tablePath, 'w');

    outfields = {'#SRCIMGNAME' 'TARGET' 'PROBE' 'VOXDIMS'...
        'CELLNO' 'CELLAREA(PIX)' 'NUCVOL(VOX)'...
        'SPOTS_NUC' 'SPOTS_CYTO' 'SIGNAL_NUC' 'SIGNAL_CYTO'...
        'EST_COUNT_NUC' 'EST_NASCENT_COUNT_NUC' 'EST_COUNT_CYTO'...
        'EST_COUNT_NUC_CLOUD' 'EST_NASCENT_COUNT_NUC_CLOUD' 'EST_COUNT_CYTO_CLOUD'};
    field_count = size(outfields, 2);
    for ii = 1:field_count
        if ii > 1; fprintf(tableHandle, '\t'); end
        fprintf(tableHandle, outfields{ii});
    end
    fprintf(tableHandle, '\n');
end

function doResultsSet(quantFilePath, tableHandle)
    %Look for RNASpotsRun in same directory...
    [mydir, fname, ~] = fileparts(quantFilePath);
    fstem = replace(fname, '_quantData', '');

    dirContents = dir(mydir);
    spotsRunPath = [];
    childCount = size(dirContents, 1);
    for c = 1:childCount
        child = dirContents(c,1);
        if ~child.isdir
            if endsWith(child.name, '_rnaspotsrun.mat')
                if startsWith(child.name, fstem)
                    spotsRunPath = [mydir filesep child.name];
                else
                    if isempty(spotsRunPath)
                        spotsRunPath = [mydir filesep child.name];
                    end
                end
            end
        end
    end

    if ~isempty(spotsRunPath)
        fprintf('\t\t-> Run found!\n');
        fprintf('\t\tQuant Data: %s\n', quantFilePath);
        fprintf('\t\tRun Summary: %s\n', spotsRunPath);
        spotsRun = RNASpotsRun.loadFrom(spotsRunPath, true);

        if ~isempty(spotsRun.img_name)
            iname = spotsRun.img_name;
        else
            iname = fstem;
        end

        if ~isempty(spotsRun.meta.type_target)
            targetName = spotsRun.meta.type_target;
        else
            targetName = '[UNKNOWN]';
        end

        if ~isempty(spotsRun.meta.type_probe)
            probeName = spotsRun.meta.type_probe;
        else
            probeName = '[UNKNOWN]';
        end

        if ~isempty(spotsRun.meta.idims_voxel)
            voxDims = spotsRun.meta.idims_voxel;
        else
            voxDims = [];
        end

        clear spotsRun

        load(quantFilePath, 'quant_results');
        cellDat = quant_results.cell_rna_data;
        cellCount = size(cellDat, 2);

        for c = 1:cellCount
            myCell = cellDat(c);
            fprintf(tableHandle, '%s\t%s\t%s', iname, targetName, probeName);
            if ~isempty(voxDims)
                fprintf(tableHandle, '\t(%d,%d,%d)', voxDims.x, voxDims.y, voxDims.z);
            else
                fprintf(tableHandle, '\t(0,0,0)');
            end
            fprintf(tableHandle, '\t%d\t%d', c, nnz(myCell.mask_cell));
            fprintf(tableHandle, '\t%d', nnz(myCell.mask_nuc));
            fprintf(tableHandle, '\t%d\t%d', myCell.spotcount_nuc, myCell.spotcount_cyto);
            fprintf(tableHandle, '\t%.3f\t%.3f', myCell.signal_nuc, myCell.signal_cyto);
            [nucCount, nucNascentCount, nucCloud, nucNascentCloud, cytoCount, cytoCloud] = estimateTargetCounts(myCell);
            fprintf(tableHandle, '\t%d\t%d', nucCount, nucNascentCount);
            fprintf(tableHandle, '\t%d', cytoCount);
            fprintf(tableHandle, '\t%d\t%d', nucCloud, nucNascentCloud);
            fprintf(tableHandle, '\t%d', cytoCloud);

            fprintf(tableHandle, '\n');
            clear myCell
        end

        clear quant_results cellDat
    else
        fprintf('\t\tQuant data found for %s, but run summary could not be found!\n', fname);
    end

end

function [nucCount, nucNascentCount, nucCloud, nucNascentCloud, cytoCount, cytoCloud] = estimateTargetCounts(myCell)
    %TODO Maybe just move this function to RNAQuant and add fields in
    %SingleCell for these values?
    nucCount = 0;
    nucNascentCount = 0;
    cytoCount = 0;

    BRIGHT_STDEVS = 2; %Number of standard deviations above which a spot is considered too bright to be a single target. 

    %Get "mature RNA" (single target) intensity
    %Collect spot intensities
    totalSpots = size(myCell.spots, 2);
    if(totalSpots > 0)
        spots = myCell.spots;
        fits = [spots.gauss_fit];
        %spotints = [fits.fitMInt];
        spotints = [fits.TotFitInt];
        inNuc = [fits.nucRNA];
        clear spots fits

%         figure(1);
%         histogram(spotints);

        iMean = mean(spotints, 'all', 'omitnan');
        iStd = std(spotints, 0, 'all', 'omitnan');
        brightnessThreshold = iMean + (BRIGHT_STDEVS * iStd);
        isTooBright = spotints >= brightnessThreshold;
        normSpots = spotints(~isTooBright);
        singleIntensity = median(normSpots, 'all', 'omitnan');
        clear normSpots

        %For now, all not-too-bright spots are counted as 1?
        counts = ones(1, totalSpots);
        counts(isTooBright) = round(spotints(isTooBright) ./ singleIntensity);

        nucCount = sum(counts(inNuc));
        nucNascentCount = sum(counts(inNuc & isTooBright));
        cytoCount = sum(counts(~inNuc));
    end

    totalClouds = size(myCell.clouds, 2);
    if(totalClouds > 0)
        clouds = myCell.clouds;
        cloudInts = [clouds.total_intensity];

        %Redo spots, omitting spots marked as part of clouds
        if(totalSpots > 0)
            spots = myCell.spots;
            inCloud = [spots.in_cloud];
            if nnz(inCloud) > 0
                nucCloud = sum(counts(inNuc & ~inCloud));
                nucNascentCloud = sum(counts(inNuc & isTooBright & ~inCloud));
                cytoCloud = sum(counts(~inNuc & ~inCloud));
            else
                nucCloud = nucCount;
                nucNascentCloud = nucNascentCount;
                cytoCloud = cytoCount;
            end
        else
            nucCloud = 0;
            nucNascentCloud = 0;
            cytoCloud = 0;

            %Need a singleIntensity and brightnessThreshold
            iMean = mean(cloudInts, 'all', 'omitnan');
            iStd = std(cloudInts, 0, 'all', 'omitnan');
            brightnessThreshold = iMean + (BRIGHT_STDEVS * iStd);
            
            %The single will just be the smallest cloud
            singleIntensity = min(cloudInts, [], 'all', 'omitnan');
        end

        %Try to get target count from clouds.
        cloudNuc = [clouds.is_nuc];
        cloudTooBright = cloudInts >= brightnessThreshold;
        cloudCounts = round(cloudInts ./ singleIntensity);

        nucCloud = nucCloud + sum(cloudCounts(cloudNuc));
        nucNascentCloud = nucNascentCloud + sum(cloudCounts(cloudNuc & cloudTooBright));
        cytoCloud = cytoCloud + sum(cloudCounts(~cloudNuc));
    else
        nucCloud = nucCount;
        nucNascentCloud = nucNascentCount;
        cytoCloud = cytoCount;
    end

end

function doDir(dirPath, tableHandle)
    %Search for files ending with _quantData.mat
    %Then look for one in the same dir ending in _rnaspotsrun.mat
    fprintf('\tProcessing %s...\n', dirPath);
    dirContents = dir(dirPath);
    childCount = size(dirContents, 1);
    for c = 1:childCount
        child = dirContents(c,1);
        if child.isdir
            if ~strcmp(child.name, '.') & ~strcmp(child.name, '..')
                doDir([dirPath filesep child.name], tableHandle);
            end
        else
            if endsWith(child.name, '_quantData.mat')
                doResultsSet([dirPath filesep child.name], tableHandle);
            end
        end
    end
end
