function Main_DumpQuantResults(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2025.03.10.02';
VERSION_STRING = 'v1.1.2';

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.

input_dir = []; %Give a dir for a batch and it will look recursively through it for TS results.
output_path = [];

skipFlaggedDups = false;

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

        if strcmp(lastkey, "skipdups")
            skipFlaggedDups = true;
            if arg_debug; fprintf("Skip Flagged Duplicates: On\n"); end
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
    if skipFlaggedDups
        output_path = [input_dir filesep 'cellCounts_skipdups.tsv'];
    else
        output_path = [input_dir filesep 'cellCounts.tsv'];
    end
end

fprintf('Main_DumpQuantResults\n');
fprintf('Script Version: %s\n', BUILD_STRING);
fprintf('TrueSpot Version: %s\n', VERSION_STRING);
fprintf('Run Time: %s\n', datetime);
fprintf('Input: %s\n', input_dir);
fprintf('Output: %s\n', output_path);

tableHandle = openOutTable(output_path);
doDir(input_dir, tableHandle, skipFlaggedDups);

fclose(tableHandle);
end

% ========================== Additional Functions ==========================

function tableHandle = openOutTable(tablePath)
    tableHandle = fopen(tablePath, 'w');

    outfields = {'#SRCIMGNAME' 'TARGET' 'PROBE' 'VOXDIMS'...
        'CELLNO' 'CELLAREA_PIX' 'NUCVOL_VOX'...
        'SPOTS_NUC' 'SPOTS_CYTO' 'SIGNAL_NUC' 'SIGNAL_CYTO'...
        'EST_COUNT_NUC' 'EST_NASCENT_COUNT_NUC' 'EST_COUNT_CYTO'...
        'EST_COUNT_NUC_CLOUD' 'EST_NASCENT_COUNT_NUC_CLOUD' 'EST_COUNT_CYTO_CLOUD' ...
        'NUC_TOT_INT' 'NUC_MED_INT' 'NUC_INT_STDEV' ...
        'NUC_MAXPROJ_AREA' 'NUC_MAXPROJ_TOT_INT' 'NUC_MAXPROJ_MED_INT' 'NUC_MAXPROJ_INT_STDEV' ...
        'CYTOVOL_VOX' 'CYTO_TOT_INT' 'CYTO_MED_INT' 'CYTO_INT_STDEV' ...
        'CYTO_MAXPROJ_AREA' 'CYTO_MAXPROJ_TOT_INT' 'CYTO_MAXPROJ_MED_INT' 'CYTO_MAXPROJ_INT_STDEV'};
    field_count = size(outfields, 2);
    for ii = 1:field_count
        if ii > 1; fprintf(tableHandle, '\t'); end
        fprintf(tableHandle, outfields{ii});
    end
    fprintf(tableHandle, '\n');
end

function doResultsSet(quantFilePath, tableHandle, skipFlaggedDups)
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
        if isfield(quant_results, 'cell_rna_data')
            cellDat = quant_results.cell_rna_data;
        elseif isfield(quant_results, 'cellData')
            quant_results = RNAQuant.readResultsSavePackage(quant_results);
            cellDat = quant_results.cell_rna_data;
        end
        cellCount = size(cellDat, 2);

        globalBrightTh = quant_results.globalBrightTh;
        globalSingleInt = quant_results.globalSingleInt;

        CSTATS_COL_COUNT = 15;

        for c = 1:cellCount
            myCell = cellDat(c);
            fprintf(tableHandle, '%s\t%s\t%s', iname, targetName, probeName);
            if ~isempty(voxDims)
                fprintf(tableHandle, '\t(%d,%d,%d)', voxDims.x, voxDims.y, voxDims.z);
            else
                fprintf(tableHandle, '\t(0,0,0)');
            end
            myCell = myCell.updateSpotAndSignalValues(globalBrightTh, globalSingleInt, skipFlaggedDups);

            fprintf(tableHandle, '\t%d\t%d', c, nnz(myCell.mask_cell));
            fprintf(tableHandle, '\t%d', nnz(myCell.mask_nuc));
            fprintf(tableHandle, '\t%d\t%d', myCell.spotcount_nuc, myCell.spotcount_cyto);
            fprintf(tableHandle, '\t%.3f\t%.3f', myCell.signal_nuc, myCell.signal_cyto);
            %[nucCount, nucNascentCount, nucCloud, nucNascentCloud, cytoCount, cytoCloud] = estimateTargetCounts(myCell);
            %myCell = myCell.updateCountEstimates();
            fprintf(tableHandle, '\t%d\t%d', myCell.nucCount, myCell.nucNascentCount);
            fprintf(tableHandle, '\t%d', myCell.cytoCount);
            fprintf(tableHandle, '\t%d\t%d', myCell.nucCloud, myCell.nucNascentCloud);
            fprintf(tableHandle, '\t%d', myCell.cytoCloud);

            if ~isempty(myCell.cell_stats)

                suffixes = {'intensity_total' 'intensity_median' 'intensity_stdev'};
                fmtwrite = {'\t%d' '\t%d' '\t%.3f'};
                ifieldcount = size(suffixes, 2);

                %Nuc3
                for i = 1:ifieldcount
                    myFieldName = ['nuc_' suffixes{i}];
                    if isfield(myCell.cell_stats, myFieldName)
                        fprintf(tableHandle, fmtwrite{i}, myCell.cell_stats.(myFieldName));
                    else
                        fprintf(tableHandle, '\tNaN');
                    end
                end

                %Nuc Max
                if isfield(myCell.cell_stats, 'nuc_max_area')
                    fprintf(tableHandle, '\t%d', myCell.cell_stats.nuc_max_area);
                else
                    fprintf(tableHandle, '\tNaN');
                end
                for i = 1:ifieldcount
                    myFieldName = ['nuc_max_' suffixes{i}];
                    if isfield(myCell.cell_stats, myFieldName)
                        fprintf(tableHandle, fmtwrite{i}, myCell.cell_stats.(myFieldName));
                    else
                        fprintf(tableHandle, '\tNaN');
                    end
                end

                %Cyto3
                if isfield(myCell.cell_stats, 'cyto_vol')
                    fprintf(tableHandle, '\t%d', myCell.cell_stats.cyto_vol);
                else
                    fprintf(tableHandle, '\tNaN');
                end
                for i = 1:ifieldcount
                    myFieldName = ['cyto_' suffixes{i}];
                    if isfield(myCell.cell_stats, myFieldName)
                        fprintf(tableHandle, fmtwrite{i}, myCell.cell_stats.(myFieldName));
                    else
                        fprintf(tableHandle, '\tNaN');
                    end
                end

                %Cyto Max
                if isfield(myCell.cell_stats, 'cyto_max_area')
                    fprintf(tableHandle, '\t%d', myCell.cell_stats.cyto_max_area);
                else
                    fprintf(tableHandle, '\tNaN');
                end
                for i = 1:ifieldcount
                    myFieldName = ['cyto_max_' suffixes{i}];
                    if isfield(myCell.cell_stats, myFieldName)
                        fprintf(tableHandle, fmtwrite{i}, myCell.cell_stats.(myFieldName));
                    else
                        fprintf(tableHandle, '\tNaN');
                    end
                end
            else
                for i = 1:CSTATS_COL_COUNT
                    fprintf(tableHandle, '\tNaN');
                end
            end

            fprintf(tableHandle, '\n');
            clear myCell
        end

        clear quant_results cellDat
    else
        fprintf('\t\tQuant data found for %s, but run summary could not be found!\n', fname);
    end

end

function doDir(dirPath, tableHandle, skipFlaggedDups)
    %Search for files ending with _quantData.mat
    %Then look for one in the same dir ending in _rnaspotsrun.mat
    fprintf('\tProcessing %s...\n', dirPath);
    dirContents = dir(dirPath);
    childCount = size(dirContents, 1);
    for c = 1:childCount
        child = dirContents(c,1);
        if child.isdir
            if ~strcmp(child.name, '.') & ~strcmp(child.name, '..')
                doDir([dirPath filesep child.name], tableHandle, skipFlaggedDups);
            end
        else
            if endsWith(child.name, '_quantData.mat')
                doResultsSet([dirPath filesep child.name], tableHandle, skipFlaggedDups);
            end
        end
    end
end
