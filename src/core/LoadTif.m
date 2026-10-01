%
%%

function [channels, idims] = LoadTif(tifpath, total_ch, ch_to_read, verbosity)
    if nargin < 2; total_ch = 0; end
    if nargin < 3; ch_to_read = []; end
    if nargin < 4; verbosity = 0; end

    %Check for bioformats...
    BFOKAY = false;
    BF_DIR = '../thirdparty/bfmatlab';
    if isfolder(BF_DIR)
        addpath(BF_DIR);
        BFOKAY = true;
    end
    
%addpath('../thirdparty');

%Verbosity:
%   0 - None
%   1 - Print local warnings, but not tiffread
%   2 - Print everything

    rfstate = RNA_Fisher_State.getStaticState();
    rfstate.tif_verbose = (verbosity > 1);
    
    if verbosity > 0
        fprintf("LoadTif -- Reading %s...\n", tifpath);
    end

    if BFOKAY
        if verbosity > 0
            fprintf("LoadTif -- Bio-Formats detected! Using bfopen...\n");
        end

        idatRaw = bfopen(tifpath);
        seriesCount = size(idatRaw, 1);
        if (seriesCount > 1) & (verbosity > 0)
            fprintf("LoadTif -- More than one series found in input. Reading only first series...\n");
        end

        stackDat = idatRaw{1, 1};
        omeMeta = idatRaw{1, 4};
        subImgCount = omeMeta.getImageCount();
        if (subImgCount > 1) & (verbosity > 0)
            fprintf("LoadTif -- More than subimage found in series. Reading only first image...\n");
        end

        idims = omeGetStackDims(omeMeta, 0);
        total_ch = omeMeta.getChannelCount(0);
        channels = cell(total_ch, 1);
        if isempty(ch_to_read)
            ch_to_read = 1:1:total_ch;
        end

        stackPos = 1;
        for c = 1:total_ch
            if ismember(c, ch_to_read)
                myStack = zeros(idims.y, idims.x, idims.z, class(stackDat{stackPos,1}));
                for z = 1:idims.z
                    myPlane = stackDat{stackPos, 1};
                    myStack(:,:,z) = myPlane(:,:);
                    stackPos = stackPos + 1;
                end
                channels{c, 1} = myStack;
                clear z myStack myPlane
            else
                stackPos = stackPos + idims.z;
            end
        end
    else
        if total_ch < 1
            total_ch = 1;
            if verbosity > 0
                fprintf("LoadTif -- Channel count was not provided. Will read as single channel.\n");
            end
        end

        [stack, img_read] = tiffread2(tifpath);
        Z = img_read/total_ch;
        Y = size(stack(1,1).data, 1);
        X = size(stack(1,1).data, 2);

        idims = struct('x', X, 'y', Y, 'z', Z);

        channels = cell(total_ch, 1);
        if isempty(ch_to_read)
            ch_to_read = 1:1:total_ch;
        end
        sz = size(ch_to_read,2);

        if verbosity > 0
            fprintf("LoadTif -- %d x %d x %d image read!\n", idims.x, idims.y, idims.z);
        end

        for i = 1:sz
            c = ch_to_read(1,i);
            j = c;
            channel = NaN(Y,X,Z);

            for z = 1:Z
                channel(:,:,z) = stack(1,j).data;
                j = j + total_ch;
            end
            channels{c,1} = channel;

            if verbosity > 0
                fprintf("LoadTif -- Channel %d loaded!\n", c);
            end
        end
    end

end

function idims = omeGetStackDims(omeMeta, imageIndex)
    idims = struct('x', 0, 'y', 0, 'z', 0);
    idims.x = omeMeta.getPixelsSizeX(imageIndex).getValue();
    idims.y = omeMeta.getPixelsSizeY(imageIndex).getValue();
    idims.z = omeMeta.getPixelsSizeZ(imageIndex).getValue();
end
