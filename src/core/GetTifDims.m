%
%%

function [X, Y, Z] = GetTifDims(tif_path, n_channels)

    [stack, img_read] = tiffread2(tif_path);
    Z = img_read/n_channels;
    Y = size(stack(1,1).data,1);
    X = size(stack(1,1).data,2);

end