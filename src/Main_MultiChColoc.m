%
%%

function Main_MultiChColoc(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2024.10.17.00';
VERSION_STRING = 'v1.1.1';

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.

param_struct.chARunPath = [];
param_struct.chBRunPath = [];

param_struct.no3d = false;

%Manual radius
param_struct.rad3 = 0;
param_struct.radZ = 0;

%Auto radius
param_struct.scanrad3 = 0; %Maximum radius to scan up to
param_struct.scanradZ = 0;
param_struct.randWeightZ = true;

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
        
        if strcmp(lastkey, "flat")
            param_struct.no3d = true;
            if arg_debug; fprintf("Collapse 3D Stacks: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "randzuni")
            param_struct.randWeightZ = false;
            if arg_debug; fprintf("Randomizer Weight Z: Off\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "randzw")
            param_struct.randWeightZ = true;
            if arg_debug; fprintf("Randomizer Weight Z: On\n"); end
            lastkey = [];
        end
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "chrunA")
            param_struct.chARunPath = argval;
            if arg_debug; fprintf("Channel A RunInfo Path Set: %s\n", param_struct.chARunPath); end
        elseif strcmp(lastkey, "chrunB")
            param_struct.chBRunPath = argval;
            if arg_debug; fprintf("Channel B RunInfo Path Set: %s\n", param_struct.chBRunPath); end
        end
    end
end %End of argin for loop


end