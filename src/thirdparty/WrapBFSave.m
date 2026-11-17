%

%Requires bio-formats

%%
function WrapBFSave(varargin)
    myImage = varargin{1};
    if iscell(myImage)
        myImage = myImage{1};
    end

    %Assumed to be YX(ZCT)
    %Because of how bfsave iterates, it ends up writing XYZCT, but always
    %gets annotated as XYCZT for some reason?

    %This wrapper rearranges the planes to trick bfsave into writing them
    %in XYCZT order.

    X = size(myImage, 2);
    Y = size(myImage, 1);
    Z = size(myImage, 3);
    C = size(myImage, 4);
    T = size(myImage, 5);

    planeCount = Z * C * T;
    myImage = reshape(myImage, Y, X, planeCount);
    pp = 1:planeCount;
    pz = repmat(1:Z, [1 C]) - 1;
    pc = floor((pp - 1) ./ Z);
    remap = (pz .* C) + pc + 1;

    myNewImage = zeros(Y,X,planeCount,class(myImage));
    myNewImage(:,:,remap) = myImage(:,:,:);
    
    myNewImage = reshape(myNewImage, Y, X, Z, C, T);

    varargin{1} = myNewImage;
    bfsave(varargin{:});
end