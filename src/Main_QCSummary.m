function Main_QCSummary(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2024.07.01.00';
VERSION_STRING = 'v1.1.0';

% ========================== Process args ==========================
arg_debug = true; %CONSTANT used for debugging arg parser.

input_dir = []; %Give a dir for a batch and it will look recursively through it for TS results.
output_dir = []; %Defaults to input dir. Place to put batch summaries.

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
            output_dir = argval;
            if arg_debug; fprintf("Output Directory Set: %s\n", output_dir); end
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

if isempty(output_dir)
    output_dir = input_dir;
end

fprintf('Main_QCSummary\n');
fprintf('Script Version: %s\n', BUILD_STRING);
fprintf('TrueSpot Version: %s\n', VERSION_STRING);
fprintf('Run Time: %s\n', datetime);
fprintf('Input: %s\n', input_dir);
fprintf('Output: %s\n', output_dir);

% ========================== Recursive Scan ==========================

thTablePath = [output_dir filesep 'batchThresholdInfo.tsv'];
tableHandle = openThTable(thTablePath);
doDir(tableHandle, input_dir);
fclose(tableHandle);

end

% ========================== Additional Functions ==========================

function tableHandle = openThTable(tablePath)
    tableHandle = fopen(tablePath, 'w');

    outfields = {'IMGNAME' 'CHANNEL_IDX' 'TARGET_NAME' 'PROBE_NAME' ...
        'TH_SCAN_MIN' 'TH_SCAN_MAX' 'TH_SELECTED' 'TH_SPOT_COUNT' ...
        'CELL_COUNT'};
    field_count = size(outfields, 2);
    for i = 1:field_count
        if i > 1; fprintf(tableHandle, '\t'); end
        fprintf(tableHandle, outfields{i});
    end
    fprintf(tableHandle, '\n');
end

function doDir(thTableHandle, dirPath)
    %Look for cellseg, rnaspotsrun/callTable, or quant results
    %For cellseg, render cell and nuc mask to png.
    %For spots, grab th info and render spot count curve to png
    %For quant, grab cell count. Not sure what else yet.
    fprintf('Scanning %s ...\n', dirPath);
    dirContents = dir(dirPath);
    childCount = size(dirContents, 1);

    csPath = [];
    srPath = [];
    qdPath = [];

    for i = 1:childCount
        fname = dirContents(i,1).name;
        isdir = dirContents(i,1).isdir;
        if isdir
            if ~strcmp(fname, '.') & ~strcmp(fname, '..')
                doDir(thTableHandle, [dirPath filesep fname]);
            end
        else
            if startsWith(fname, 'CellSeg_')
                csPath = [dirPath filesep fname];
            elseif endsWith(fname, '_rnaspotsrun.mat')
                srPath = [dirPath filesep fname];
            elseif endsWith(fname, '_quantData.mat')
                qdPath = [dirPath filesep fname];
            end
        end
    end

    if ~isempty(csPath)
        load(csPath, 'cellSeg', 'nucleiSeg');
        
        figHandle = figure(1);
        clf;
        imshow(cellSeg.cell_mask, []);
        saveas(figHandle, [dirPath filesep 'cellMask.png']);
        close(figHandle);

        figHandle = figure(1);
        clf;
        imshow(nucleiSeg.results.nuc_label, []);
        saveas(figHandle, [dirPath filesep 'nucLabel.png']);
        close(figHandle);

        clear cellSeg nucleiSeg
    end

    if ~isempty(srPath)
        spotsrun = RNASpotsRun.loadFrom(srPath, true);
        fprintf(thTableHandle, '%s', spotsrun.img_name);
        fprintf(thTableHandle, '\t%d', spotsrun.channels.rna_ch);
        fprintf(thTableHandle, '\t%s', spotsrun.meta.type_target);
        fprintf(thTableHandle, '\t%s', spotsrun.meta.type_probe);

        fprintf(thTableHandle, '\t%d', spotsrun.options.t_min);
        fprintf(thTableHandle, '\t%d', spotsrun.options.t_max);
        fprintf(thTableHandle, '\t%d', spotsrun.intensity_threshold);

        %Now need to load call table...
        [~, callTable] = spotsrun.loadCallTable();
        figHandle = renderSpotPlot(spotsrun, callTable);
        saveas(figHandle, [dirPath filesep 'scLog.png']);
        close(figHandle);

        if spotsrun.intensity_threshold > 0
            thSpotCount = nnz(callTable{:, 'dropout_thresh'} >= spotsrun.intensity_threshold);
            fprintf(thTableHandle, '\t%d', thSpotCount);
            clear thSpotCount
        else
            fprintf(thTableHandle, '\tNaN');
        end
        clear callTable spotsrun

        if ~isempty(qdPath)
            load(qdPath, 'quant_results');
            cellCount = size(quant_results.cell_rna_data, 2);
            fprintf(thTableHandle, '\t%d', cellCount);
            clear quant_results
        else
            fprintf(thTableHandle, '\tNOQUANT');
        end

        fprintf(thTableHandle, '\n');
    end

end

function figHandle = renderSpotPlot(spotsrun, callTable)

    color = [0.667 0.220 0.220];

    [x, y] = RNAUtils.spotCountFromCallTable(callTable, false, spotsrun.options.t_min, spotsrun.options.t_max);

    figHandle = figure(1);
    clf;

    y = double(y);
    y = log10(y);
    ymax = max(y, [], 'all', 'omitnan');
    thval = 0;

    hold on;
    if ~isempty(spotsrun.threshold_results)
        score_list = RNAThreshold.getAllThresholdSuggestions(spotsrun.threshold_results);
        score_mean = mean(score_list, 'all', 'omitnan');
        score_std = std(score_list, 0, 'all', 'omitnan');
        xmin = score_mean - score_std;
        xmax = score_mean + score_std;
        if xmin < 0; xmin = 0; end

        cdiff = [1 1 1] - color;
        cdiff = cdiff ./ 2;
        boxcolor = color + cdiff;

        rectangle('Position', [xmin 0 (xmax - xmin) ymax],...
            'FaceColor', boxcolor, 'LineStyle', 'none');

        thval = spotsrun.threshold_results.threshold;
    end

    plot(x, y, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', color);

    if thval > 0
        xline(thval, 'LineStyle', '--', 'LineWidth', 1.5);
    end

    if ~isempty(spotsrun.threshold_results)
        score_min = min(score_list, [], 'all', 'omitnan');
        score_max = max(score_list, [], 'all', 'omitnan');
        if(score_min < 0); score_min = 0; end
        xline(score_min, 'LineStyle', ':', 'LineWidth', 1);
        xline(score_max, 'LineStyle', ':', 'LineWidth', 1);
    end

end

