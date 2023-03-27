%
%%

function Test_CleanupBigfish(dirpath)
    %This cleans up any bigfish runs that finished, but didn't import.

    addpath('..');
    versionString = "2023.03.27.02";
    fprintf('Test_CleanupBigfish | Version %s\n', versionString);
    fprintf('Input directory: %s\n', dirpath);
    
    doDirectory(dirpath);
end

function doDirectory(dirpath)

    dirContents = dir(dirpath);
    contentCount = size(dirContents,1);

    fprintf('> Processing directory: %s\n', dirpath);
    asBFDir = false;
    slurmOut = [];
    %Look for summary.txt
    for i = 1:contentCount
        dirEntry = dirContents(i);
        
        cPath = [dirpath filesep dirEntry.name];
        if dirEntry.isdir
            if strcmp(dirEntry.name, '.'); continue; end
            if strcmp(dirEntry.name, '..'); continue; end
            doDirectory(cPath);
        else
            if strcmp(dirEntry.name, 'summary.txt')
                asBFDir = true;
            end
            
            if endsWith(dirEntry.name, '_slurm.out')
                slurmOut = dirEntry.name;
            end
        end
    end
    
    if ~asBFDir; return; end
    fprintf('\tsummary.txt found!\n');
    
    %Look for stem...
    if ~isempty(slurmOut)
        %Use this file name
        bfstem = replace(slurmOut, '_slurm.out', '');
    else
        %Use BIGFISH_dirname
        [~, dname, ~] = fileparts(dirpath);
        bfstem = ['BIGFISH_' dname];
    end
    
    %Run import
    bfstem_full = [dirpath filesep bfstem];
    ctpath = [bfstem_full '_coordTable.mat'];
    %Check for existing run
    fprintf('\tChecking for %s ...\n', ctpath);
    if isfile(ctpath)
        finfo = who('-file', ctpath);
        if ~isempty(find(ismember(finfo, 'coord_table'),1))
            fprintf('\tRun already found for %s! No import needed.\n', bfstem_full);
            return;
        end
    end
    
    Main_Bigfish2Mat(dirpath, bfstem_full);
    
    %Delete csv files if successful
    if isfile(ctpath)
        finfo = who('-file', ctpath);
        if ~isempty(find(ismember(finfo, 'coord_table'),1))
            load(ctpath, 'coord_table');
            if ~isempty(coord_table)
                if size(coord_table,1) > 1
                    fprintf('\tImport successful! Removing raw csvs...\n');
                    for i = 1:contentCount
                        dirEntry = dirContents(i);
                        
                        if dirEntry.isdir; continue; end
                        cPath = [dirpath filesep dirEntry.name];
                        if endsWith(dirEntry.name, '.csv')
                            if startsWith(dirEntry.name, 'spots_')
                                delete(cPath);
                            end
                        end
                    end
                end
            end
            clear coord_table;
        end
    end

end