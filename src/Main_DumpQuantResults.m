function Main_DumpQuantResults(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2024.05.20.00';
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

end

% ========================== Additional Functions ==========================

function tableHandle = openOutTable(tablePath)
    tableHandle = fopen(tablePath, 'w');

    outfields = {'#SRCIMGNAME' 'TARGET' 'PROBE' ...
        'CELLNO' 'CELLSIZE' 'NUCSIZE'...
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