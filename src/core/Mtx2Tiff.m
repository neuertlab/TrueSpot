%
%%

%Source: https://www.mathworks.com/matlabcentral/answers/707933-how-to-save-a-series-of-images-tiff-as-a-stack

function success = Mtx2Tiff(my_image, outpath)

    success = true;
    Y = size(my_image,1);
    X = size(my_image,2);

    out_tiff = Tiff(outpath, 'w');
    tagstruct.ImageLength = Y;
    tagstruct.ImageWidth = X;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    dimcount = ndims(my_image);
    
    if dimcount == 2
        setTag(out_tiff, tagstruct);
        write(out_tiff,squeeze(my_image));
    elseif dimcount == 3
        Z = size(my_image,3);
        for z = 1:Z
            if z > 1
                writeDirectory(out_tiff);
            end
            setTag(out_tiff, tagstruct);
            write(out_tiff,squeeze(my_image(:,:,z)));
        end
    elseif dimcount == 4
        C = size(my_image,4);
        for c = 1:C
            if c > 1
                writeDirectory(out_tiff);
            end
            Z = size(my_image,3);
            for z = 1:Z
                if z > 1
                    writeDirectory(out_tiff);
                end
                setTag(out_tiff, tagstruct);
                write(out_tiff,squeeze(my_image(:,:,z)));
            end
        end
    else
        success = false;
    end
    
    close(out_tiff);
end