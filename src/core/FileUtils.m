%
%%
classdef FileUtils

    methods(Static)

        %---------------------------- Internal ---------------------------

        %%
        function tempList = internal_initTempList(initCapacity)
            tempList = struct();
            tempList.nextIndex = 1;
            tempList.capacity = initCapacity;
            tempList.list = cell(1, initCapacity);
        end

        %%
        function tempList = internal_expandTempList(tempList, additionalCapacity)
            newCap = tempList.capacity + additionalCapacity;
            oldList = tempList.list;
            lsize = tempList.nextIndex - 1;
            tempList.list = cell(1, newCap);
            if lsize > 0
                tempList.list(1:lsize) = oldList(1:lsize);
            end
            tempList.capacity = newCap;
        end

        %%
        function tempList = internal_appendToTempList(tempList, item, additionalCapacity)
            if nargin < 3; additionalCapacity = 64; end
            if (tempList.nextIndex > tempList.capacity)
                tempList = FileUtils.internal_expandTempList(tempList, additionalCapacity);
            end
            tempList.list{tempList.nextIndex} = item;
            tempList.nextIndex = tempList.nextIndex + 1;
        end

        %%
        function tempList = internal_scanForFilesEndingWith(dirPath, ptnString, boolRecursive, tempList)
            dirList = dir(dirPath);
            childCount = size(dirList, 1);
            for i = 1:childCount
                childItem = dirList(i);
                childName = childItem.name;
                childPath = [dirPath filesep childName];
                if childItem.isdir
                    if boolRecursive
                        if ~strcmp(childName, '.') & ~strcmp(childName, '..')
                            tempList = FileUtils.internal_scanForFilesEndingWith(childPath, ptnString, boolRecursive, tempList);
                        end
                    end
                else
                    if endsWith(childName, ptnString)
                        tempList = FileUtils.internal_appendToTempList(tempList, childPath);
                    end
                end
            end
        end

        %%
        function tempList = internal_scanForFilesWithName(dirPath, fileName, boolRecursive, tempList)
            dirList = dir(dirPath);
            childCount = size(dirList, 1);
            for i = 1:childCount
                childItem = dirList(i);
                childName = childItem.name;
                childPath = [dirPath filesep childName];
                if childItem.isdir
                    if boolRecursive
                        if ~strcmp(childName, '.') & ~strcmp(childName, '..')
                            tempList = FileUtils.internal_scanForFilesWithName(childPath, fileName, boolRecursive, tempList);
                        end
                    end
                else
                    if strcmp(childName, fileName)
                        tempList = FileUtils.internal_appendToTempList(tempList, childPath);
                    end
                end
            end
        end

        %---------------------------- Interface ---------------------------

        %%
        function pathList = scanForFilesEndingWith(dirPath, ptnString, boolRecursive)
            pathList = [];
            tempList = FileUtils.internal_initTempList(64);
            tempList = FileUtils.internal_scanForFilesEndingWith(dirPath, ptnString, boolRecursive, tempList);
            lsize = tempList.nextIndex - 1;

            if lsize < 1; return; end

            pathList = tempList.list(1:lsize);
        end

        %%
        function pathList = scanForFilesWithName(dirPath, fileName, boolRecursive)
            pathList = [];
            tempList = FileUtils.internal_initTempList(64);
            tempList = FileUtils.internal_scanForFilesWithName(dirPath, fileName, boolRecursive, tempList);
            lsize = tempList.nextIndex - 1;

            if lsize < 1; return; end

            pathList = tempList.list(1:lsize);
        end

    end
end