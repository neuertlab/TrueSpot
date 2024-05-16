%
%%

function [channels, idims] = LoadTif(tifpath, total_ch, ch_to_read, verbosity)

%addpath('../thirdparty');

%Verbosity:
%   0 - None
%   1 - Print local warnings, but not tiffread
%   2 - Print everything

    rfstate = RNA_Fisher_State.getStaticState();
    rfstate.tif_verbose = (verbosity > 1);
    
    if verbosity > 0
        fprintf("LoadTif -- Reading %s\n", tifpath);
    end

    [stack, img_read] = tiffread2(tifpath);
    Z = img_read/total_ch;
    Y = size(stack(1,1).data,1);
    X = size(stack(1,1).data,2);
    
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