%
%%

function Bkg_Mask_Core(light_ch, seg_path, out_stem, outputpng)

    bkge = Background_Extractor;
    bkge = bkge.initialize(light_ch, seg_path);
    bkge = bkge.setBkgMaskSavePath([out_stem '.mat']);
    bkge = bkge.applyAndSave();
    
    if outputpng
        bkge.exportMaskedImage([out_stem '.png']);
        bkge.exportHistogramImage([out_stem '_histo.png']);
    end

end