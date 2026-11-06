%
%%

function opStruct = LoadTifPlusStructGen()
    opStruct = struct();
    opStruct.verbosity = 0;
    opStruct.chCount = 0;

    opStruct.channelsToRead = [];
    opStruct.seriesToRead = [];
    opStruct.subimagesToRead = [];

    opStruct.flipZCAxes = false; %For some reason the getPixelsDimensionOrder return string is not consistent???
end