%
%%

function Main_SpotVisGUI(varargin)

addpath('./core');
addpath('./thirdparty');

% ========================== Process args ==========================
arg_debug = true; %CONSTANT used for debugging arg parser.

inputPath = [];
tifPath = []; %In case tif not linked properly in runinfo
cellSegPath = []; %Optional

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
            inputPath = argval;
            if arg_debug; fprintf("Input Path Set: %s\n", inputPath); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

%--- Check args (Fill in defaults based on inputs)

if isempty(inputPath)
    fprintf('Please provide an input runinfo file!\n');
    return;
end

end