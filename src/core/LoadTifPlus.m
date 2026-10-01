%
%%

function loadResults = LoadTifPlus(tifpath, opStruct)
    if nargin < 2; opStruct = []; end

    %Check for bioformats...
    BFOKAY = false;
    BF_DIR = './thirdparty/bfmatlab';
    if isfolder(BF_DIR)
        addpath(BF_DIR);
        BFOKAY = true;
    end

    loadResults = struct();
    loadResults.inputPath = tifpath;
    loadResults.imageData = [];
    loadResults.imageDims = [];
    loadResults.voxelDims = [];
    loadResults.bioFormatsDetected = BFOKAY;
    loadResults.inputType = [];

    if isempty(opStruct)
        opStruct = LoadTifPlusStructGen();
    end
    
%Verbosity:
%   0 - None
%   1 - Print local warnings, but not tiffread
%   2 - Print everything

    rfstate = RNA_Fisher_State.getStaticState();
    rfstate.tif_verbose = (opStruct.verbosity > 1);
    
    if opStruct.verbosity > 0
        fprintf("LoadTifPlus -- Reading %s\n", tifpath);
    end

    if BFOKAY
        if opStruct.verbosity > 0
            fprintf("LoadTifPlus -- Bio-Formats detected! Using bfopen...\n");
        end

        idatRaw = bfopen(tifpath);
        seriesCount = size(idatRaw, 1);

        if isempty(opStruct.seriesToRead)
            opStruct.seriesToRead = 1:1:seriesCount;
        end

        if ~isempty(opStruct.subimagesToRead)
            slotCount = size(opStruct.subimagesToRead, 2);
            if slotCount < seriesCount
                tempcell = cell(1, seriesCount);
                tempcell(1:slotCount) = opStruct.subimagesToRead;
                opStruct.subimagesToRead = tempcell;
            end
            clear slotCount tempcell
        end

        loadResults.imageData = cell(1, seriesCount);
        loadResults.imageDims = cell(1, seriesCount);
        loadResults.voxelDims = cell(1, seriesCount);
        loadResults.inputType = cell(1, seriesCount);
        for s = 1:seriesCount
            if ~ismember(s, opStruct.seriesToRead)
                continue;
            end
            if opStruct.verbosity > 0
                fprintf("LoadTifPlus -- Processing series %d of %d...\n", s, seriesCount);
            end
            subincl = [];
            if ~isempty(opStruct.subimagesToRead)
                subincl = opStruct.subimagesToRead{s};
            end

            stackDat = idatRaw{s, 1};
            omeMeta = idatRaw{s, 4};
            subImgCount = omeMeta.getImageCount();
            stackPos = 1;

            if isempty(subincl)
                subincl = 1:1:subImgCount;
            end

            s_idat = cell(1, subImgCount);
            s_class = cell(1, subImgCount);
            s_idims(subImgCount) = struct('x', 0, 'y', 0, 'z', 0);
            s_vdims(subImgCount) = struct('x', 0, 'y', 0, 'z', 0);
            
            for i = 1:subImgCount
                if ~ismember(i, subincl)
                    continue;
                end

                if opStruct.verbosity > 0
                    fprintf("LoadTifPlus -- Processing series subimage %d of %d...\n", i, subImgCount);
                end

                idims = omeGetStackDims(omeMeta, i-1);
                vdims = omeGetVoxelSizeNano(omeMeta, i-1);
                chCount = omeMeta.getChannelCount(i-1);

                ichDat = cell(1, chCount);
                chincl = opStruct.channelsToRead;
                if isempty(chincl)
                    chincl = 1:1:chCount;
                end

                itype = class(stackDat{stackPos,1});
                for c = 1:chCount
                    if ismember(c, chincl)
                        myStack = zeros(idims.y, idims.x, idims.z, itype);
                        for z = 1:idims.z
                            myPlane = stackDat{stackPos, 1};
                            myStack(:,:,z) = myPlane(:,:);
                            stackPos = stackPos + 1;
                        end
                        ichDat{c} = myStack;
                        clear z myStack myPlane
                    else
                        stackPos = stackPos + idims.z;
                    end
                end
                clear c

                s_idims(i) = idims;
                s_vdims(i) = vdims;
                s_idat{i} = ichDat;
                s_class{i} = itype;
                clear idims vdims chCount ichDat chincl itype
            end
            clear i

            loadResults.imageData{s} = s_idat;
            loadResults.imageDims{s} = s_idims;
            loadResults.voxelDims{s} = s_vdims;
            loadResults.inputType{s} = s_class;

            clear s_idat s_idims s_vdims subincl stackDat omeMeta subImgCount stackPos s_class
        end
    else
        if opStruct.chCount < 1
            opStruct.chCount = 1;
            if opStruct.verbosity > 0
                fprintf("LoadTifPlus -- Channel count was not provided. Will read as single channel.\n");
            end
        end

        [stack, img_read] = tiffread2(tifpath);
        test_slice = stack(1,1).data;
        Z = img_read/opStruct.chCount;
        Y = size(test_slice, 1);
        X = size(test_slice, 2);

        s_idat = cell(1, 1);
        s_idims = struct('x', X, 'y', Y, 'z', Z);
        s_vdims = struct('x', 0, 'y', 0, 'z', 0);
        s_class = cell(1, 1);
        s_class{1} = class(test_slice);

        channels = cell(1, opStruct.chCount);
        if isempty(opStruct.channelsToRead)
            opStruct.channelsToRead = 1:1:opStruct.chCount;
        end
        sz = size(opStruct.channelsToRead,2);

        if opStruct.verbosity > 0
            fprintf("LoadTifPlus -- %d x %d x %d image loaded!\n", s_idims.x, s_idims.y, s_idims.z);
        end

        for i = 1:sz
            c = opStruct.channelsToRead(i);
            j = c;
            channel = NaN(Y,X,Z);

            for z = 1:Z
                channel(:,:,z) = stack(1,j).data;
                j = j + opStruct.chCount;
            end
            channels{c} = channel;

            if opStruct.verbosity > 0
                fprintf("LoadTifPlus -- Channel %d loaded!\n", c);
            end
        end

        s_idat{1} = channels;

        loadResults.imageData = cell(1,1);
        loadResults.imageDims = cell(1,1);
        loadResults.voxelDims = cell(1,1);
        loadResults.inputType = cell(1,1);

        loadResults.imageData{1} = s_idat;
        loadResults.imageDims{1} = s_idims;
        loadResults.voxelDims{1} = s_vdims;
        loadResults.inputType{1} = s_class;
    end

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