%
%%
function Main_AnalyzeBatchThresholds(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2025.04.04.00';
VERSION_STRING = 'v1.2.0';

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.

input_dir = []; %Give a dir for a batch and it will look recursively through it for TS results.
table_out_path = [];
stats_out_path = [];

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
            input_dir = char(argval);
            if arg_debug; fprintf("Input Directory Set: %s\n", input_dir); end
        elseif strcmp(lastkey, "tableout")
            table_out_path = argval;
            if arg_debug; fprintf("Table Output Path Set: %s\n", table_out_path); end
        elseif strcmp(lastkey, "statsout")
            stats_out_path = argval;
            if arg_debug; fprintf("Stats Output Path Set: %s\n", stats_out_path); end
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

if isempty(table_out_path)
    table_out_path = [input_dir filesep 'thTable.tsv'];
end
if isempty(stats_out_path)
    stats_out_path = [input_dir filesep 'thStats.txt'];
end

fprintf('Main_AnalyzeBatchThresholds\n');
fprintf('Script Version: %s\n', BUILD_STRING);
fprintf('TrueSpot Version: %s\n', VERSION_STRING);
fprintf('Run Time: %s\n', datetime);
fprintf('Input: %s\n', input_dir);
fprintf('Table Output: %s\n', table_out_path);
fprintf('Stats Output: %s\n', stats_out_path);

% ========================== Scan for and sort runs ==========================
GLALLOC = 8;
runList = scanDirRec(input_dir);
groupList = cell(1,GLALLOC);
groupCount = 0;

runCount = size(runList, 2);
for r = 1:runCount
    runPath = runList{r};
    spotsRun = RNASpotsRun.loadFrom(runPath);
    gno = 0;
    if groupCount > 0
        for g = 1:groupCount
            groupStruct = groupList{g};
            if runIsInGroup(spotsRun, groupStruct)
                gno = g;
                groupList{g} = addToGroup(spotsRun, groupStruct);
                break;
            end
        end
    end

    if gno == 0
        %New group
        groupStruct = newGroup(spotsRun);
        capacity = size(groupList, 2);
        if groupCount >= capacity
            oldList = groupList;
            groupList = cell(1, capacity + GLALLOC);
            groupList(1:groupCount) = oldList(:);
        end
        groupCount = groupCount + 1;
        groupList{groupCount} = groupStruct;
    end
end

% ========================== Output ==========================
tableHandle = fopen(table_out_path, 'w');
statsTxtHandle = fopen(stats_out_path, 'w');

%Table header
tblHeaderFields = {'IMGNAME' 'TARGET' 'PROBE' 'CHANNEL' ...
    'AUTO_TH' 'PRESET_TEST_THS' 'THPOOL_MEAN' 'THPOOL_MEDIAN' ...
    'THPOOL_STD' 'THPOOL_MAD'};
field_count = size(tblHeaderFields, 2);
for ii = 1:field_count
    if ii > 1; fprintf(tableHandle, '\t'); end
    fprintf(tableHandle, tblHeaderFields{ii});
end
fprintf(tableHandle, '\n');
clear field_count ii tblHeaderFields

%Go through groups
for g = 1:groupCount
    groupStruct = groupList{g};
    groupThs = NaN(1, groupStruct.memberCount);
    for m = 1:groupStruct.memberCount
        spotsRun = groupStruct.members{m};
        fprintf(tableHandle, '%s', spotsRun.img_name);
        if ~isempty(spotsRun.meta.type_target)
            fprintf(tableHandle, '\t%s', spotsRun.meta.type_target);
        else
            fprintf(tableHandle, '\t<None>');
        end
        if ~isempty(spotsRun.meta.type_probe)
            fprintf(tableHandle, '\t%s', spotsRun.meta.type_probe);
        else
            fprintf(tableHandle, '\t<None>');
        end
        fprintf(tableHandle, '\t%d', spotsRun.channels.rna_ch);

        fprintf(tableHandle, '\t%d', spotsRun.intensity_threshold);
        groupThs(m) = spotsRun.intensity_threshold;

        if ~isempty(spotsRun.th_alt)
            if isfield(spotsRun.th_alt, 'thPresetSugg')
                presetCount = size(spotsRun.th_alt.thPresetSugg, 1); 
                for p = 1:presetCount
                    if p > 1
                        fprintf(tableHandle, ';');
                    else
                        fprintf(tableHandle, '\t');
                    end
                    fprintf(tableHandle, '%d', spotsRun.th_alt.thPresetSugg(p,1));
                end

                fprintf(tableHandle, '\t%.3f', mean(spotsRun.th_alt.thPresetSugg(:,1), 'all', 'omitnan'));
                fprintf(tableHandle, '\t%.1f', median(spotsRun.th_alt.thPresetSugg(:,1), 'all', 'omitnan'));
                fprintf(tableHandle, '\t%.3f', std(spotsRun.th_alt.thPresetSugg(:,1), 0, 'all', 'omitnan'));
                fprintf(tableHandle, '\t%.3f', mad(spotsRun.th_alt.thPresetSugg(:,1), 1, 'all'));
            else
                fprintf(tableHandle, '\t<Not evaluated>\tNaN\tNaN\tNaN\tNaN');
            end
        else
            fprintf(tableHandle, '\t<Not evaluated>\tNaN\tNaN\tNaN\tNaN');
        end
        fprintf(tableHandle, '\n');
    end

    %Group summary
    groupStruct.stats_mean = mean(groupThs, 'all', 'omitnan');
    groupStruct.stats_median = median(groupThs, 'all', 'omitnan');
    groupStruct.stats_stdev = std(groupThs, 0, 'all', 'omitnan');
    groupStruct.stats_mad = mad(groupThs, 1, 'all');
    groupStruct.stats_min = min(groupThs, [], 'all', 'omitnan');
    groupStruct.stats_max = max(groupThs, [], 'all', 'omitnan');
    fprintf(statsTxtHandle, '================================================================\n');
    fprintf(statsTxtHandle, 'Target: %s\n', groupStruct.targetName);
    fprintf(statsTxtHandle, 'Probe: %s\n', groupStruct.probeName);
    fprintf(statsTxtHandle, 'Channel: %d\n', groupStruct.channelNumber);
    fprintf(statsTxtHandle, 'Image Count: %d\n', groupStruct.memberCount);
    fprintf(statsTxtHandle, '\tThreshold Mean: %.3f\n', groupStruct.stats_mean);
    fprintf(statsTxtHandle, '\tThreshold Median: %.1f\n', groupStruct.stats_median);
    fprintf(statsTxtHandle, '\tThreshold StDev: %.3f\n', groupStruct.stats_stdev);
    fprintf(statsTxtHandle, '\tThreshold MAD: %.3f\n', groupStruct.stats_mad);
    fprintf(statsTxtHandle, '\tThreshold Min: %d\n', groupStruct.stats_min);
    fprintf(statsTxtHandle, '\tThreshold Max: %d\n', groupStruct.stats_max);

    fprintf(statsTxtHandle, '\n');
end

fclose(tableHandle);
fclose(statsTxtHandle);
end

% ========================== Additional Functions ==========================

function groupStruct = newGroup(spotsRun)
    groupStruct = genRunGroupStruct();
    groupStruct.probeName = spotsRun.meta.type_probe;
    groupStruct.targetName = spotsRun.meta.type_target;
    groupStruct.channelNumber = spotsRun.channels.rna_ch;
    groupStruct = addToGroup(spotsRun, groupStruct);
end

function groupStruct = addToGroup(spotsRun, groupStruct)
    if isempty(groupStruct.members)
        groupStruct.members = cell(1,64);
        groupStruct.memberCount = 1;
        groupStruct.members{1} = spotsRun;
        return;
    end

    capacity = size(groupStruct.members, 2);
    if groupStruct.memberCount >= capacity
        %Reallocate, yay!
        newList = cell(1, capacity*2);
        newList(1:capacity) = groupStruct.members(:);
        groupStruct.members = newList;
        clear newList
    end

    groupStruct.memberCount = groupStruct.memberCount + 1;
    groupStruct.members{groupStruct.memberCount} = spotsRun;
end

function boolres = runIsInGroup(spotsRun, groupStruct)
    boolres = false;
    if isempty(groupStruct); return; end
    if isempty(spotsRun); return; end

    %If metadata are there, check that.
    %If neither probe nor target name are there, just match channel
    if ~isempty(groupStruct.targetName) & ~isempty(spotsRun.meta.type_target)
        if strcmp(groupStruct.targetName, spotsRun.meta.type_target)
            if ~isempty(groupStruct.probeName) & ~isempty(spotsRun.meta.type_probe)
                if strcmp(groupStruct.probeName, spotsRun.meta.type_probe)
                    %Channel numbers do not need to match. Though they probably
                    %SHOULD.
                    boolres = true;
                    return;
                end
            end
        end
    end

    if ~isempty(groupStruct.probeName) & ~isempty(spotsRun.meta.type_probe)
        %If probe but no target name, then channel numbers should match
        if strcmp(groupStruct.probeName, spotsRun.meta.type_probe)
            if groupStruct.channelNumber == spotsRun.channels.rna_ch
                boolres = true;
                return;
            end
        end
    end
    
    if groupStruct.channelNumber == spotsRun.channels.rna_ch
        boolres = true;
        return;
    end
end

function rgrpStruct = genRunGroupStruct()
    rgrpStruct = struct();
    rgrpStruct.probeName = [];
    rgrpStruct.targetName = [];
    rgrpStruct.channelNumber = 0;
    rgrpStruct.members = [];
    rgrpStruct.memberCount = 0;

    rgrpStruct.stats_mean = NaN;
    rgrpStruct.stats_median = NaN;
    rgrpStruct.stats_stdev = NaN;
    rgrpStruct.stats_mad = NaN;
    rgrpStruct.stats_min = NaN;
    rgrpStruct.stats_max = NaN;
end

function runList = scanDirRec(dirPath)
    fprintf('Scanning %s ...\n', dirPath);
    dirContents = dir(dirPath);
    childCount = size(dirContents, 1);

    srPath = [];

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
                    subRuns = scanDirRec(subdirPath);
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
            if endsWith(fname, '_rnaspotsrun.mat')
                srPath = [dirPath filesep fname];
            end
        end
    end

    if ~isempty(srPath)
        runsFound = runsFound + 1;
    end

    %Condense raw list into list of valid run structs
    if runsFound > 0
        runList = cell(1,runsFound);
        index = 1;
        if ~isempty(srPath) 
            runList{1} = srPath;
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
