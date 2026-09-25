%%
function Main_CZI2TIFBatch(varargin)
addpath('./core');
addpath('./thirdparty');

BF_DIR = './thirdparty/bfmatlab';
if isfolder(BF_DIR)
    addpath(BF_DIR);
else
    fprintf(['Bio-Formats was not found! Please download MATLAB version from ' ...
        'https://www.openmicroscopy.org/bio-formats/downloads and copy to ' BF_DIR]);
end

BUILD_STRING = '2026.09.25.02';
VERSION_STRING = 'v1.3.3';

% ========================== Process args ======================

fprintf('Running Main_CZI2TIFBatch build %s (TrueSpot %s)!\n', BUILD_STRING, VERSION_STRING);

arg_debug = true; %CONSTANT used for debugging arg parser.
opStruct = genOptionsStruct();

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
        
        %Account for boolean keys...
        if strcmp(lastkey, "skip2d")
            opStruct.skip2D = true;
            if arg_debug; fprintf("Skip 2D Inputs: On\n"); end
            lastkey = [];
        elseif strcmp(lastkey, "pngrender")
            opStruct.pngRender = true;
            if arg_debug; fprintf("Output PNG MIPs: On\n"); end
            lastkey = [];
        end
  
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "indir")
            opStruct.inDirPath = argval;
            if arg_debug; fprintf("Input Path Set: %s\n", opStruct.inDirPath); end
        elseif strcmp(lastkey, "outdir")
            opStruct.outDirPath = argval;
            if arg_debug; fprintf("Output Directory Set: %s\n", opStruct.outDirPath); end
        elseif strcmp(lastkey, "xmlout")
            opStruct.outXmlPath = argval;
            if arg_debug; fprintf("Output XML Path Set: %s\n", opStruct.outXmlPath); end
        elseif strcmp(lastkey, "incl")
            opStruct.inPatternMatch = argval;
            if arg_debug; fprintf("Input File Name Requirement Set: %s\n", opStruct.inPatternMatch); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

%--- Check args (Fill in defaults based on inputs)
if isempty(opStruct.inDirPath)
    fprintf("ERROR: Input directory required! Exiting...\n");
    return;
end

if ~isfolder(opStruct.inDirPath)
    fprintf("ERROR: Please provide a directory as input. Exiting...\n");
    return;
end

if isempty(opStruct.outDirPath)
    opStruct.outDirPath = opStruct.inDirPath;
    fprintf("Output directory was not explicitly provided. Set to %s...\n", opStruct.outDirPath);
end

if ~isfolder(opStruct.outDirPath)
    mkdir(opStruct.outDirPath);
end

if isempty(opStruct.outXmlPath)
    opStruct.outXmlPath = [opStruct.outDirPath filesep 'batchSettings.xml'];
    fprintf("Output XML path was not explicitly provided. Set to %s...\n", opStruct.outXmlPath);
end

%--- Gather a list of image files that meet input specifications
%   1. They must be files
%   2. Their name must end in .czi (case insensitive)
%   3. Their name must contain whatever string the user input, if
%   applicable

dlist = dir(opStruct.inDirPath);
allocAmt = size(dlist, 1);
fnames = cell(1, allocAmt);
ct = 0;
for i = 1:allocAmt
    fitem = dlist(i);
    if ~fitem.isdir
        if endsWith(lower(fitem.name), '.czi')
            if ~isempty(opStruct.inPatternMatch)
                if contains(fitem.name, opStruct.inPatternMatch, 'IgnoreCase', true)
                    ct = ct + 1;
                    fnames{ct} = fitem.name;
                end
            else
                ct = ct + 1;
                fnames{ct} = fitem.name;
            end
        end
    end
end

if (ct < 1)
    fprintf("No valid czi files were found in the input directory! Exiting...\n");
    return;
end

fnames = fnames(1:ct);
clear ct allocAmt fitem dlist i

iCount = size(fnames, 2);
fprintf("%d czi images found!\n", iCount);

%--- Read images, convert, and sort into batches by metadata

batchCount = 0;
batchList = cell(1,iCount);
cdList = cell(1,iCount);
outTable = cell(iCount*8, 2); %Output image file name, batch index
tifCount = 0;

for i = 1:iCount
    fprintf("Working on %s...\n", fnames{i});
    fstem = replace(fnames{i}, '.czi', '');
    fstem = replace(fstem, ' ', '_');
    fstem = replace(fstem, '.', '_');
    fstem = replace(fstem, '(', '');
    fstem = replace(fstem, ')', '');
    fstem = replace(fstem, '+', '_');
    inpath = [opStruct.inDirPath filesep fnames{i}];

    fprintf("\tReading %s...\n", inpath);
    try
        idatRaw = bfopen(inpath);
    catch ME
        fprintf("\tWARNING: Failed to read %s! File will be skipped!\n", inpath);
        continue;
    end

    seriesCount = size(idatRaw, 1);
    for s = 1:seriesCount
        stackDat = idatRaw{s, 1};
        omeMeta = idatRaw{s, 4};
        subImgCount = omeMeta.getImageCount();
        stackPos = 1;
        for j = 1:subImgCount
            if (subImgCount == 1) & (seriesCount == 1)
                fostem = fstem;
            else
                fostem = [fstem '_s' num2str(s) 'i' num2str(j)];
            end
            foname = [fostem '.tif'];
            outpath = [opStruct.outDirPath filesep foname];
            idims = omeGetStackDims(omeMeta, j-1);
            if opStruct.skip2D & (idims.z < 2)
                continue;
            end

            chCount = omeMeta.getChannelCount(j-1);

            myStack = zeros(idims.y, idims.x, idims.z, chCount, class(stackDat{stackPos,1}));
            % for z = 1:idims.z
            %     for c = 1:chCount
            %         myPlane = stackDat{stackPos, 1};
            %         myStack(:,:,z,c) = myPlane(:,:);
            %         stackPos = stackPos + 1;
            %     end
            % end

            for c = 1:chCount
                for z = 1:idims.z
                    myPlane = stackDat{stackPos, 1};
                    myStack(:,:,z,c) = myPlane(:,:);
                    stackPos = stackPos + 1;
                end
            end

            fprintf("\tGenerating %s...\n", outpath);
            if isfile(outpath)
                delete(outpath);
            end
            bfsave(myStack, outpath);

            %PNG renders, if applicable
            if opStruct.pngRender
                stackPos = 1;
                for z = 1:idims.z
                    for c = 1:chCount
                        myPlane = stackDat{stackPos, 1};
                        myStack(:,:,z,c) = myPlane(:,:);
                        stackPos = stackPos + 1;
                    end
                end

                for c = 1:chCount
                    pngOut = [fostem '_c' num2str(c) '.png'];
                    pngPath = [opStruct.outDirPath filesep pngOut];
                    chDat = myStack(:,:,:,c);
                    cmip = max(chDat, [], 3, 'omitnan');

                    if isfile(pngPath)
                        delete(pngPath);
                    end

                    fh = figure(10);
                    clf;
                    imshow(cmip, []);
                    saveas(fh, pngPath);
                    close(fh);
                end
                clear pngOut c cmip fh chDat pngPath
            end

            %Clean up
            clear myStack myPlane c z

            %Assign to batch
            assignToBatch = -1;
            if batchCount > 0
                for b = 1:batchCount
                    if imageBelongsInBatch(omeMeta, j-1, batchList{b}, cdList{b})
                        assignToBatch = b;
                        break;
                    end
                end
            end

            if assignToBatch < 1
                %New batch
                batchCount = batchCount + 1;
                assignToBatch = batchCount;
                myBatch = TrueSpotProjSettings.newBatchSettings();
                myChannels(chCount) = TrueSpotChannelSettings.newChannelSettings();
                for c = 1:(chCount-1)
                    myChannels(c) = TrueSpotChannelSettings.newChannelSettings();
                end

                myBatch.metadata.voxelSize = omeGetVoxelSizeNano(omeMeta, j-1);
                myBatch.metadata.batchName = ['batch' num2str(assignToBatch)];
                myBatch.channelInfo.channelCount = chCount;

                for c = 1:chCount
                    myChannels(c).channelIndex = c;
                    myChannels(c).metadata.targetName = omeMeta.getChannelName(j-1, c-1);
                    myChannels(c).metadata.probeName = omeMeta.getChannelFluor(j-1, c-1);
                end

                batchList{assignToBatch} = myBatch;
                cdList{assignToBatch} = myChannels;

                clear myBatch myChannels
            end

            tifCount = tifCount + 1;
            outTable{tifCount, 1} = foname;
            outTable{tifCount, 2} = assignToBatch;

            clear assignToBatch outpath idims chCount foname
        end
    end
end

batchList = batchList(1:batchCount);
cdList = cdList(1:batchCount);
outTable = outTable(1:tifCount, :);

if tifCount < 1
    fprintf("No valid stacks found! Exiting...\n");
    return;
end

fprintf("%d valid stacks successfully converted!\n", tifCount);

%--- Move output TIFs to correct batch subdirectories (if applicable)
%Also update the batch settings input/output info

outBatches = zeros(1, tifCount);
for i = 1:tifCount
    outBatches(i) = outTable{i, 2};
end
clear i

if batchCount > 1
    for b = 1:batchCount
        myBatch = batchList{b};
        batchDir = [opStruct.outDirPath filesep myBatch.metadata.batchName];
        myBatch.paths.inputPath = batchDir;
        myBatch.paths.outputPath = batchDir;
        batchList{b} = myBatch;

        if ~isfolder(batchDir)
            mkdir(batchDir);
        end

        ibool = (outBatches == b);
        ilist = outTable(ibool, 1);
        icount = size(ilist, 1);
        for i = 1:icount
            srcPath = [opStruct.outDirPath filesep ilist{i}];
            dstPath = [batchDir filesep ilist{i}];
            if isfile(dstPath)
                delete(dstPath);
            end
            movefile(srcPath, dstPath);
        end

        clear ibool ilist icount myBatch batchDir srcPath dstPath
    end
else
    myBatch = batchList{1};
    myBatch.paths.inputPath = opStruct.outDirPath;
    myBatch.paths.outputPath = opStruct.outDirPath;
    batchList{1} = myBatch;

    clear myBatch
end


%--- Output XML file from metadata gathered by batches
fprintf("\tOutputting XML...\n");
if isfile(opStruct.outXmlPath)
    delete(opStruct.outXmlPath);
end
TrueSpotXML.writeSettingsXML(opStruct.outXmlPath, batchList, cdList);

end

% ========================== Helper functions ======================

function opStruct = genOptionsStruct()
    opStruct = struct();

    opStruct.inDirPath = [];
    opStruct.outDirPath = [];
    opStruct.outXmlPath = [];

    opStruct.inPatternMatch = []; %File name must contain this pattern to be included
    opStruct.skip2D = false;
    opStruct.pngRender = false;
end

function resbool = imageBelongsInBatch(omeMeta, imageIndex, batchStruct, batchChannels)
    % Must match:
    %   - Channel count
    %   - Pixel/Voxel dims
    %   - All channel names
    %   - All channel fluor

    resbool = false;

    cCount = omeMeta.getChannelCount(imageIndex);
    if cCount ~= batchStruct.channelInfo.channelCount
        return;
    end

    vdims = omeGetVoxelSizeNano(omeMeta, imageIndex);
    if (vdims.x ~= batchStruct.metadata.voxelSize.x); return; end
    if (vdims.y ~= batchStruct.metadata.voxelSize.y); return; end
    if (vdims.z ~= batchStruct.metadata.voxelSize.z); return; end

    chNames = omeGetChannelNames(omeMeta, imageIndex);
    for i = 1:cCount
        if ~strcmp(chNames{i}, batchChannels(i).metadata.targetName); return; end
    end

    chFluor = omeGetChannelFluors(omeMeta, imageIndex);
    for i = 1:cCount
        if ~strcmp(chFluor{i}, batchChannels(i).metadata.probeName); return; end
    end

    resbool = true;
end

function idims = omeGetStackDims(omeMeta, imageIndex)
    idims = struct('x', 0, 'y', 0, 'z', 0);
    idims.x = omeMeta.getPixelsSizeX(imageIndex).getValue();
    idims.y = omeMeta.getPixelsSizeY(imageIndex).getValue();
    idims.z = omeMeta.getPixelsSizeZ(imageIndex).getValue();
end

function idims = omeGetVoxelSizeNano(omeMeta, imageIndex)
    idims = struct('x', 0, 'y', 0, 'z', 0);
    idims.x = double(omeMeta.getPixelsPhysicalSizeX(imageIndex).value(ome.units.UNITS.NANOMETER));
    idims.y = double(omeMeta.getPixelsPhysicalSizeY(imageIndex).value(ome.units.UNITS.NANOMETER));
    idims.z = double(omeMeta.getPixelsPhysicalSizeZ(imageIndex).value(ome.units.UNITS.NANOMETER));

    %Round
    idims.x = uint32(round(idims.x));
    idims.y = uint32(round(idims.y));
    idims.z = uint32(round(idims.z));
end

function chNames = omeGetChannelNames(omeMeta, imageIndex)
    cCount = omeMeta.getChannelCount(imageIndex);
    chNames = cell(1, cCount);
    for i = 1:cCount
        j = i - 1;
        chNames{i} = omeMeta.getChannelName(imageIndex, j);
    end
end

function chFluor = omeGetChannelFluors(omeMeta, imageIndex)
    cCount = omeMeta.getChannelCount(imageIndex);
    chFluor = cell(1, cCount);
    for i = 1:cCount
        j = i - 1;
        chFluor{i} = omeMeta.getChannelFluor(imageIndex, j);
    end
end