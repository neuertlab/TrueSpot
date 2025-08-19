%
%%

function Main_ImportCellSeg(varargin)
addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2025.07.31.00';
VERSION_STRING = 'v1.3.1';

% ========================== Process args ==========================

arg_debug = true; %CONSTANT used for debugging arg parser.
cellMaskPath = [];
nucMaskPath = [];
outputPath = [];
nucZMin = 1;

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
        if strcmp(lastkey, "cellmask")
            cellMaskPath = argval;
            if arg_debug; fprintf("Cell Mask Path Set: %s\n", cellMaskPath); end
        elseif strcmp(lastkey, "nucmask")
            nucMaskPath = argval;
            if arg_debug; fprintf("Nuclear Mask Path Set: %s\n", nucMaskPath); end
        elseif strcmp(lastkey, "out")
            outputPath = argval;
            if arg_debug; fprintf("Output Path Set: %s\n", outputPath); end
        elseif strcmp(lastkey, "nzmin")
            nucZMin = Force2Num(argval);
            if arg_debug; fprintf("Nuclear Mask Min Z Set: %d\n", nucZMin); end
        end
    end
end %End of argin for loop

% ========================== Run ==========================

fprintf("-------------- TRUESPOT Import Cell Segmentation Data --------------\n");
fprintf("TrueSpot Version: %s\n", VERSION_STRING);
fprintf("ImportCellSeg Version: %s\n", BUILD_STRING);
fprintf("Timestamp: %s\n", datetime);

if isempty(cellMaskPath) & isempty(nucMaskPath)
    fprintf("At least one input file is required! Exiting... \n");
    return;
end

if isempty(outputPath)
    if ~isempty(cellMaskPath)
        [outdir, cmfile, ~] = fileparts(cellMaskPath);
        outputPath = [outdir filesep cmfile '.mat'];
    elseif ~isempty(nucMaskPath)
        [outdir, nmfile, ~] = fileparts(nucMaskPath);
        outputPath = [outdir filesep nmfile '.mat'];
    else
        fprintf("At least one input file is required! Exiting... \n");
        return;
    end
    fprintf("WARNING: Output path was not provided. Set to %s \n", outputPath);
end

CellSeg.importCellSegData(cellMaskPath, nucMaskPath, outputPath, nucZMin);

fprintf("Done: %s\n", datetime);

end