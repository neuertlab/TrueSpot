%
%%

function opStruct = LoadTifPlusStructGen()
    opStruct = struct();
    opStruct.verbosity = 0;
    opStruct.chCount = 0;

    opStruct.channelsToRead = [];
    opStruct.seriesToRead = [];
    opStruct.subimagesToRead = [];
end