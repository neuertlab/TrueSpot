%
%%
function loadResults = LoadTifPlusSingle(tifpath, opStruct, seriesId, subimageId)
    if nargin < 2; opStruct = []; end
    if nargin < 3; seriesId = 1; end
    if nargin < 4; subimageId = 1; end

    if isempty(opStruct)
        opStruct = LoadTifPlusStructGen();
    end

    opStruct.seriesToRead = seriesId;
    opStruct.subimagesToRead{seriesId} = subimageId;
    
    fullRes = loadTifPlus(tifpath, opStruct);
    loadResults.inputPath = fullRes.inputPath;
    loadResults.bioFormatsDetected = fullRes.bioFormatsDetected;
    s_idat = fullRes.imageData{seriesId};
    loadResults.imageData = s_idat{subimageId};
    s_idims = fullRes.imageDims{seriesId};
    loadResults.imageDims = s_idims(subimageId);
    s_vdims = fullRes.voxelDims{seriesId};
    loadResults.voxelDims = s_vdims(subimageId);
    s_types = fullRes.inputType{seriesId};
    loadResults.inputType = s_types{subimageId};

end