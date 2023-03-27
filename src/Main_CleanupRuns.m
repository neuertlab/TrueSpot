%
%%

function Main_CleanupRuns(directoryPath, deleteDebug)
    versionStr = '2023.03.27.03';
    
    addpath('./core');
    deleteDebug = Force2Bool(deleteDebug);

    fprintf('Main_CleanupRuns | Version %s\n', versionStr);
    fprintf('------------------------------------------------\n');

    fprintf('Input directory: %s\n', directoryPath);
    if ~isfolder(directoryPath)
        fprintf('[ERROR] %s does not exist! Exiting...\n', directoryPath);
        return;
    end

    cleanupDir(directoryPath, deleteDebug);

end

function cleanupDir(dirPath, delDebug)
    dirContents = dir(dirPath);
    contentCount = size(dirContents,1);

    fprintf('> Processing directory: %s\n', dirPath);
    for i = 1:contentCount
        dirEntry = dirContents(i);
        
        cPath = [dirPath filesep dirEntry.name];
        if dirEntry.isdir
            if strcmp(dirEntry.name, '.'); continue; end
            if strcmp(dirEntry.name, '..'); continue; end
            cleanupDir(cPath, delDebug);
        else
            if endsWith(dirEntry.name, 'coordTable.mat')
                %See if it's intact
                finfo = who('-file', cPath);
                if ~isempty(find(ismember(finfo, 'coord_table'),1))
                    load(cPath, 'coord_table');
                    if isempty(coord_table)
                        delete(cPath);
                        fprintf('\tMalformed coord table removed: %s\n', cPath);
                    end
                    clear coord_table;
                else
                    %Bad file. Delete.
                    delete(cPath);
                    fprintf('\tMalformed coord table removed: %s\n', cPath);
                end
            end
            if endsWith(dirEntry.name, 'ZTrimmed.mat')
                %Can be regenned when needed. Clear for now.
                delete(cPath);
                fprintf('\tRegeneratable ZTrimmed table removed: %s\n', cPath);
            end
            if delDebug & endsWith(dirEntry.name, '.png')
                %Can be regenned when needed. Clear for now.
                delete(cPath);
                fprintf('\tDebug PNG output removed: %s\n', cPath);
            end
            if delDebug & endsWith(dirEntry.name, '.eps')
                %Can be regenned when needed. Clear for now.
                delete(cPath);
                fprintf('\tDebug EPS output removed: %s\n', cPath);
            end
        end
    end
    
    %If dir is empty at the end, delete the dir
    dirContents = dir(dirPath);
    contentCountNew = size(dirContents,1);
    
    if contentCountNew <= 2
        rmdir(dirPath);
        fprintf('\tEmpty directory removed: %s\n', dirPath);
    end
end
