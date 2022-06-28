%
%%

function Main_RNASpotsMetaAdd(dirpath, varargin)

    %Look for a valid RNASpotsRun file in the given directory, and load.
% ========================== Scan directory ==========================

DirContents = dir(dirpath);
fcount = size(DirContents,1);

rundatpath = [];
for i = 1:fcount
    fname = DirContents(i,1).name;
    if ~DirContents(i,1).isdir && endsWith(fname, "_rnaspotsrun.mat")
        rundatpath = fname;
        break;
    end
end

if isempty(rundatpath)
    fprintf("No run data found in %s. Exiting...\n", dirpath);
    return;
end

fprintf("Run data found: %s\n", rundatpath);
spotsrun = RNASpotsRun.loadFrom([dirpath filesep fname]);

% ========================== Handle remaining args ==========================

lastkey = [];
for i = 1:nargin-1
    argval = varargin{i};
    if startsWith(argval, "-")
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
        if strcmp(lastkey, "probetype")
            spotsrun.type_probe = argval;
            fprintf("Probe Set: %s\n", spotsrun.type_probe);
        elseif strcmp(lastkey, "target")
            spotsrun.type_target = argval;
            fprintf("Target Set: %s\n", spotsrun.type_target);
        elseif strcmp(lastkey, "targettype")
            spotsrun.type_targetmol = argval;
            fprintf("Target Type Set: %s\n", spotsrun.type_targetmol);
        elseif strcmp(lastkey, "species")
            spotsrun.type_species = argval;
            fprintf("Species Set: %s\n", spotsrun.type_species);
        elseif strcmp(lastkey, "celltype")
            spotsrun.type_cell = argval;
            fprintf("Cell Type Set: %s\n", spotsrun.type_cell);
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

spotsrun.saveMe();

end