function Main_DumpQuantResults(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2025.04.10.00';
VERSION_STRING = 'v1.2.0';

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.

input_dir = []; %Give a dir for a batch and it will look recursively through it for TS results.
output_path = [];
xml_path = [];

opsStruct = genOpsStruct();
opsStruct.buildString = BUILD_STRING;

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
            opsStruct.skipFlaggedDups = true;
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
        elseif strcmp(lastkey, "chdef")
            xml_path = argval;
            if arg_debug; fprintf("Channel/Sample Definition Doc Set: %s\n", xml_path); end
        elseif strcmp(lastkey, "mthresh")
            opsStruct.mThresh = parseNumberList(argval);
            if arg_debug
                fprintf("Manual threshold(s) set: ");
                dbgPrintNumberList(opsStruct.mThresh);
                fprintf("\n");
            end
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
    if opsStruct.skipFlaggedDups
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
fprintf('XML Def: %s\n', xml_path);

%Determine threshold values
thGroupNames = [];
if ~isempty(opsStruct.mThresh)
    thCount = size(opsStruct.mThresh, 2);
    thGroupNames = cell(1, thCount);
    for i = 1:thCount
        thGroupNames{i} = num2str(opsStruct.mThresh(i));
    end
else
    if ~isempty(xml_path)
        opsStruct.channelDefs = MultichParams.readChannelDefXML(xml_path);
    end
    thGroupNames = {'MID' 'LO' 'HI'};
end

tableHandle = openOutTable(output_path, thGroupNames);
doDir(input_dir, tableHandle, opsStruct);

fclose(tableHandle);
end

% ========================== Additional Functions ==========================

function opsStruct = genOpsStruct()
    opsStruct = struct();
    opsStruct.channelDefs = [];
    opsStruct.skipFlaggedDups = false;
    opsStruct.mThresh = [];
    opsStruct.buildString = [];
end

function tableHandle = openOutTable(tablePath, thNames)
    tableHandle = fopen(tablePath, 'w');

    outfields_cellcmn = {'SRCIMGNAME' 'TARGET' 'PROBE' 'VOXDIMS'...
        'CELLNO' 'CELLAREA_PIX' 'NUCVOL_VOX' ...
        'NUC_TOT_INT' 'NUC_MED_INT' 'NUC_INT_STDEV' ...
        'NUC_MAXPROJ_AREA' 'NUC_MAXPROJ_TOT_INT' 'NUC_MAXPROJ_MED_INT' 'NUC_MAXPROJ_INT_STDEV' ...
        'CYTOVOL_VOX' 'CYTO_TOT_INT' 'CYTO_MED_INT' 'CYTO_INT_STDEV' ...
        'CYTO_MAXPROJ_AREA' 'CYTO_MAXPROJ_TOT_INT' 'CYTO_MAXPROJ_MED_INT' 'CYTO_MAXPROJ_INT_STDEV'};

    thfields = {'THVAL' 'EST_COUNT_NUC' 'EST_NASCENT_COUNT_NUC' 'EST_COUNT_CYTO' ...
        'EST_COUNT_NUC_CLOUD' 'EST_NASCENT_COUNT_NUC_CLOUD' 'EST_COUNT_CYTO_CLOUD'};

%     outfields = {'SRCIMGNAME' 'TARGET' 'PROBE' 'VOXDIMS'...
%         'CELLNO' 'CELLAREA_PIX' 'NUCVOL_VOX'...
%         'SPOTS_NUC' 'SPOTS_CYTO' 'SIGNAL_NUC' 'SIGNAL_CYTO'...
%         'EST_COUNT_NUC' 'EST_NASCENT_COUNT_NUC' 'EST_COUNT_CYTO'...
%         'EST_COUNT_NUC_CLOUD' 'EST_NASCENT_COUNT_NUC_CLOUD' 'EST_COUNT_CYTO_CLOUD' ...
%         'NUC_TOT_INT' 'NUC_MED_INT' 'NUC_INT_STDEV' ...
%         'NUC_MAXPROJ_AREA' 'NUC_MAXPROJ_TOT_INT' 'NUC_MAXPROJ_MED_INT' 'NUC_MAXPROJ_INT_STDEV' ...
%         'CYTOVOL_VOX' 'CYTO_TOT_INT' 'CYTO_MED_INT' 'CYTO_INT_STDEV' ...
%         'CYTO_MAXPROJ_AREA' 'CYTO_MAXPROJ_TOT_INT' 'CYTO_MAXPROJ_MED_INT' 'CYTO_MAXPROJ_INT_STDEV'};
    field_count = size(outfields_cellcmn, 2);
    for ii = 1:field_count
        if ii > 1; fprintf(tableHandle, '\t'); end
        fprintf(tableHandle, outfields_cellcmn{ii});
    end

    thfieldCount = size(thfields, 2);
    if ~isempty(thNames)
        thgroupCount = size(thNames, 2);
        for gg = 1:thgroupCount
            gName = thNames{gg};
            for ii = 1:thfieldCount
                if ii > 1; fprintf(tableHandle, '\t'); end
                fprintf(tableHandle, [thfields{ii} '_' gName]);
            end
        end
    else
        for ii = 1:thfieldCount
            if ii > 1; fprintf(tableHandle, '\t'); end
            fprintf(tableHandle, thfields{ii});
        end
    end

    fprintf(tableHandle, '\n');
end

function doResultsSet(quantFilePath, tableHandle, opsStruct)
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

        load(quantFilePath, 'quant_results', 'runMeta');
        if isfield(quant_results, 'cell_rna_data')
            cellDat = quant_results.cell_rna_data;
        elseif isfield(quant_results, 'cellData')
            quant_results = RNAQuant.readResultsSavePackage(quant_results);
            cellDat = quant_results.cell_rna_data;
        end
        cellCount = size(cellDat, 2);

        %globalBrightTh = quant_results.globalBrightTh;
        %globalSingleInt = quant_results.globalSingleInt;

        CSTATS_COL_COUNT = 15;

        thList = [];
        if ~isempty(opsStruct.mThresh)
            thList = opsStruct.mThresh;
        else
            if ~isempty(opsStruct.channelDefs)
                %Figure out which channel it belongs to.
                sampleInfo = MultichParams.matchRunToChannel(spotsRun, opsStruct.channelDefs);
                if ~isempty(sampleInfo)
                    thList = [sampleInfo.thMid sampleInfo.thLow sampleInfo.thHigh];
                end
            end

            %If no channel defs or couldn't match...
            if isempty(thList)
                %Pull from run
                thMid = spotsRun.intensity_threshold;
                if ~isempty(spotsRun.th_alt) & isfield(spotsRun.th_alt, 'thPresetSugg')
                    spool = spotsRun.th_alt.thPresetSugg(:, 1);
                    thLow = round(prctile(spool, 25, 'all'));
                    thHigh = round(prctile(spool, 75, 'all'));
                    clear spool
                else
                    thLow = max(thMid - 20, 10);
                    thHigh = thMid + 50;
                end
                thList = [thMid thLow thHigh];
                clear thMid thLow thHigh
            end
        end
        T = size(thList, 2);

        globalBrightTh = NaN(1, T);
        globalSingleInt = NaN(1, T);
        for t = 1:T
            quant_results = RNAQuant.updateQuantCounts(quant_results, opsStruct.skipFlaggedDups, thList(t));
            globalBrightTh(t) = quant_results.globalBrightTh;
            globalSingleInt(t) = quant_results.globalSingleInt;
        end
        quant_results = RNAQuant.results2SavePackage(quant_results);
        runMeta.modifiedDate = datetime;
        runMeta.tsDumpCountsBuild = opsStruct.buildString;
        save(quantFilePath, 'quant_results', 'runMeta');
        clear quant_results

        for c = 1:cellCount
            myCell = cellDat(c);

            %Common to cell
            fprintf(tableHandle, '%s\t%s\t%s', iname, targetName, probeName);
            if ~isempty(voxDims)
                fprintf(tableHandle, '\t(%d,%d,%d)', voxDims.x, voxDims.y, voxDims.z);
            else
                fprintf(tableHandle, '\t(0,0,0)');
            end
            fprintf(tableHandle, '\t%d\t%d', c, nnz(myCell.mask_cell));
            fprintf(tableHandle, '\t%d', nnz(myCell.mask_nuc));

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

            %Counts per tested threshold
            for t = 1:T
                thVal = thList(t);
                myCell = myCell.updateSpotAndSignalValues(globalBrightTh(t), globalSingleInt(t), opsStruct.skipFlaggedDups, thVal);

                %fprintf(tableHandle, '\t%d\t%d\t%d', thVal, myCell.spotcount_nuc, myCell.spotcount_cyto);
                %fprintf(tableHandle, '\t%.3f\t%.3f', myCell.signal_nuc, myCell.signal_cyto);
                %[nucCount, nucNascentCount, nucCloud, nucNascentCloud, cytoCount, cytoCloud] = estimateTargetCounts(myCell);
                %myCell = myCell.updateCountEstimates();
                fprintf(tableHandle, '\t%d\t%d', myCell.nucCount, myCell.nucNascentCount);
                fprintf(tableHandle, '\t%d', myCell.cytoCount);
                fprintf(tableHandle, '\t%d\t%d', myCell.nucCloud, myCell.nucNascentCloud);
                fprintf(tableHandle, '\t%d', myCell.cytoCloud);
            end

            fprintf(tableHandle, '\n');
            clear myCell
        end

        clear cellDat
    else
        fprintf('\t\tQuant data found for %s, but run summary could not be found!\n', fname);
    end

end

function doDir(dirPath, tableHandle, opsStruct)
    %Search for files ending with _quantData.mat
    %Then look for one in the same dir ending in _rnaspotsrun.mat
    fprintf('\tProcessing %s...\n', dirPath);
    dirContents = dir(dirPath);
    childCount = size(dirContents, 1);
    for c = 1:childCount
        child = dirContents(c,1);
        if child.isdir
            if ~strcmp(child.name, '.') & ~strcmp(child.name, '..')
                doDir([dirPath filesep child.name], tableHandle, opsStruct);
            end
        else
            if endsWith(child.name, '_quantData.mat')
                doResultsSet([dirPath filesep child.name], tableHandle, opsStruct);
            end
        end
    end
end

function numberList = parseNumberList(inputString)
    %Expects comma delimited decimal numbers - may or may not be surrounded
    %by parentheses
    inputString = replace(inputString, '(', '');
    inputString = replace(inputString, ')', '');
    rawList = split(inputString, ',');
    count = size(rawList, 1);
    numberList = zeros(1, count);
    for i = 1:count
        numberList(i) = str2double(rawList{i, 1});
    end
    numberList = int32(numberList);
end

function dbgPrintNumberList(numberList)
    count = size(numberList, 2);
    if count > 1
        fprintf('(');
    end
    for i = 1:count
        if i > 1; fprintf(','); end
        fprintf('%d', numberList(i));
    end
    if count > 1
        fprintf(')');
    end
end
