%
%%

function Main_ImageStats(varargin)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2026.10.02.00';
VERSION_STRING = 'v1.4.0';

% ========================== Process args ==========================

fprintf('Running Main_ImageStats build %s (TrueSpot %s)!\n', BUILD_STRING, VERSION_STRING);

arg_debug = true; %CONSTANT used for debugging arg parser.
opStruct = genOptionsStruct();

lastkey = [];
for i = 1:nargin
    argval = varargin{i};
    if ischar(argval) & startsWith(argval, "-")
        %Key
        if size(argval,2) >= 2
            lastkey = argval(2:end);
        else
            lastkey = [];
        end
        
        %Account for boolean keys...
        if strcmp(lastkey, "nodpc")
            opStruct.nodpc = true;
            if arg_debug; fprintf("Dead Pixel Clean: Off\n"); end
            lastkey = [];
        end
  
    else
        if isempty(lastkey)
            fprintf("Value without key: %s - Skipping...\n", argval);
            continue;
        end
        
        %Value
        if strcmp(lastkey, "input")
            opStruct.inputPath = argval;
            if arg_debug; fprintf("Input Path Set: %s\n", opStruct.inputPath); end
        elseif strcmp(lastkey, "ostats")
            opStruct.outputPath = argval;
            if arg_debug; fprintf("Output Directory Set: %s\n", opStruct.outputPath); end
        elseif strcmp(lastkey, "oreport")
            opStruct.reportPath = argval;
            if arg_debug; fprintf("Output XML Path Set: %s\n", opStruct.reportPath); end
        elseif strcmp(lastkey, "chtotal")
            opStruct.inputChannelCount = Force2Num(argval);
            if arg_debug; fprintf("Input Channel Count Set: %d\n", opStruct.inputChannelCount); end
        elseif strcmp(lastkey, "iseries")
            opStruct.inputSeries = Force2Num(argval);
            if arg_debug; fprintf("Input Series ID Set: %d\n", opStruct.inputSeries); end
        elseif strcmp(lastkey, "isubi")
            opStruct.inputSeriesImage = Force2Num(argval);
            if arg_debug; fprintf("Input Series Subimage ID Set: %d\n", opStruct.inputSeriesImage); end
        elseif strcmp(lastkey, "gaussrad")
            opStruct.gaussRad = Force2Num(argval);
            if arg_debug; fprintf("LoG Gaussian Radius Set: %d\n", opStruct.gaussRad); end
        elseif strcmp(lastkey, "cellmask")
            opStruct.cellSegPath = argval;
            if arg_debug; fprintf("Cell Mask Path Set: %s\n", opStruct.cellSegPath); end
        elseif strcmp(lastkey, "nucmask")
            opStruct.nucSegPath = argval;
            if arg_debug; fprintf("Nuclear Mask Path Set: %s\n", opStruct.nucSegPath); end
        else
            fprintf("Key not recognized: %s - Skipping...\n", lastkey);
        end
    end
end

%--- Check args (Fill in defaults based on inputs)
if isempty(opStruct.inputPath)
    fprintf("ERROR: Input directory required! Exiting...\n");
    return;
end

if ~isfile(opStruct.inputPath)
    fprintf("ERROR: Input file ""%s"" does not exist!\n", opStruct.inputPath);
    return;
end

[indir, infn, ~] = fileparts(opStruct.inputPath);
if isempty(opStruct.outputPath)
    opStruct.outputPath = [indir filesep infn '_stats.mat'];
    fprintf("Warning: Output file path was not provided. Set to: %s\n", opStruct.outputPath);
end

if isempty(opStruct.reportPath)
    opStruct.reportPath = [indir filesep infn '_report.txt'];
    fprintf("Warning: Output report file path was not provided. Set to: %s\n", opStruct.reportPath);
end

%--- Load image
fprintf("Loading source image...\n");
loadTifOps = LoadTifPlusStructGen();
loadTifOps.chCount = opStruct.inputChannelCount;
loadTifOps.verbosity = 1;
loadTifRes = LoadTifPlusSingle(opStruct.inputPath, loadTifOps, opStruct.inputSeries, opStruct.inputSeriesImage);
clear loadTifOps
idims = loadTifRes.imageDims;
chCount = size(loadTifRes.imageData, 2);

%--- Load cell/nuc segmentation masks (if provided)
nucMask = [];
cellMask = [];
if ~isempty(opStruct.cellSegPath)
    fprintf("Loading cell mask(s)...\n");
    if isfile(opStruct.cellSegPath)
        if endsWith(opStruct.cellSegPath, '.mat')
            %Assumed TS output. Load both cell and nuc masks.
            cellMask = CellSeg.openCellMask(opStruct.cellSegPath);
            nucMask = CellSeg.openNucMask(opStruct.cellSegPath);
        else
            cellMask = CellSeg.openCellMask(opStruct.cellSegPath, idims.y);
        end

        if ~isempty(cellMask)
            %Extend to 3D
            cellMask = repmat(cellMask, [1 1 idims.x]);
        end
    else
        fprintf("Warning: Cell mask ""%s"" does not exist. Skipping...\n", opStruct.cellSegPath);
    end
end

if isempty(nucMask) & ~isempty(opStruct.nucSegPath)
    fprintf("Loading nuclear mask...\n");
    if isfile(opStruct.nucSegPath)
        nucMask = CellSeg.openNucMask(opStruct.nucSegPath, 2, true, idims.y);
    else
        fprintf("Warning: Nuclear mask ""%s"" does not exist. Skipping...\n", opStruct.nucSegPath);
    end
end

%--- Run stats
cellCount = 0;
if ~isempty(cellMask)
    cellCount = max(cellMask, [], 'all', 'omitnan');
end

imageStats = struct();
imageStats.scriptVer = BUILD_STRING;
imageStats.tsVer = VERSION_STRING;
imageStats.modifiedDate = datetime;
imageStats.sourcePath = opStruct.inputPath;
imageStats.channelStats = opStruct.genChannelStatsStruct(cellCount);
for c = 1:chCount
    fprintf("Working on channel %d of %d...\n", c, chCount);

    %Fetch channel
    chDat = loadTifRes.imageData{c};

    %Apply DPC and LoG filter
    fprintf("\tApplying filters...\n");
    [chDat_f, deadpixInfo] = RNAUtils.applyLoGFilter(chDat, opStruct.gaussRad, ~opStruct.nodpc);
    if ~opStruct.nodpc
        chDat = RNAUtils.cleanDeadPixels(chDat, deadpixInfo, true);
    end
    clear deadpixInfo

    myChannelStats = genChannelStatsStruct(cellCount);

    %Whole channel
    myChannelStats.stats_raw = takeStackRegionStats(chDat, []);
    myChannelStats.stats_fLoG = takeStackRegionStats(chDat_f, []);

    %Cell-less bkg (if applicable)
    if cellCount > 0
        myChannelStats.bkgStats.stats_raw = takeStackRegionStats(chDat, (cellMask == 0));
        myChannelStats.bkgStats.stats_fLoG = takeStackRegionStats(chDat_f, (cellMask == 0));
    end

    %Cell by cell (if applicable)
    for i = 1:cellCount
        myCellStats = myChannelStats.cellStats(i);
        myCell = (cellMask == i);

        myCellStats.stats_raw = takeStackRegionStats(chDat, myCell);
        myCellStats.stats_fLoG = takeStackRegionStats(chDat_f, myCell);

        if ~isempty(nucMask)
            myNuc = and(myCell, nucMask);
            myCellStats.nuc.stats_raw = takeStackRegionStats(chDat, myNuc);
            myCellStats.nuc.stats_fLoG = takeStackRegionStats(chDat_f, myNuc);

            myCyto = and(myCell, ~nucMask);
            myCellStats.cyto.stats_raw = takeStackRegionStats(chDat, myCyto);
            myCellStats.cyto.stats_fLoG = takeStackRegionStats(chDat_f, myCyto);
            clear myNuc myCyto
        end

        myChannelStats.cellStats(i) = myCellStats;
        clear myCellStats myCell
    end

    imageStats.channelStats{c} = myChannelStats;
    clear myChannelStats chDat chDat_f i
end
clear c

%--- Clean up a bit 
clear nucMask cellMask idims chCount loadTifRes cellCount

%--- Output primary file
save(opStruct.outputPath, 'imageStats');

%--- Output summary report
%TODO

end

% ========================== Helper functions ======================

function opStruct = genOptionsStruct()
    opStruct = struct();

    opStruct.inputPath = []; %TIF path
    opStruct.outputPath = []; %mat path
    opStruct.reportPath = []; %txt path

    opStruct.inputChannelCount = 0;
    opStruct.inputSeries = 1;
    opStruct.inputSeriesImage = 1;

    opStruct.gaussRad = 7;
    opStruct.nodpc = false;

    opStruct.cellSegPath = []; %cellseg mask path (mask, may or may not be TS output)
    opStruct.nucSegPath = []; %nucseg mask path (mask, may or may not be TS output, though not needed if cellSegPath is TS output)
end

function statsStruct = genChannelStatsStruct(cellCount)
    statsStruct = struct();

    %Overall image
    statsStruct.stats_raw = [];
    statsStruct.stats_fLoG = [];

    if cellCount > 0
        statsStruct.bkgStats = genCellStatsStruct(false);
        statsStruct.cellStats(cellCount) = genCellStatsStruct(true);
    else
        statsStruct.bkgStats = [];
        statsStruct.cellStats = [];
    end
end

function statsStruct = genCellStatsStruct(inclnuc)
    statsStruct = struct();

    statsStruct.stats_raw = [];
    statsStruct.stats_fLoG = [];

    if inclnuc
        statsStruct.nuc = genCellStatsStruct(false);
        statsStruct.cyto = genCellStatsStruct(false);
    else
        statsStruct.nuc = [];
        statsStruct.cyto = [];
    end
end

function statsStruct = genStackRegionStatsStruct()
    statsStruct = struct();

    statsStruct.stackStats = [];
    statsStruct.maxProjStats = [];
    statsStruct.sliceStats = [];
end

function statsStruct = takeStackRegionStats(imgDat, mask)
    statsStruct = genStackRegionStatsStruct();
    if isempty(imgDat); return; end

    itype = class(imgDat);
    Z = size(imgDat, 3);
    %Y = size(imgDat, 1);
    %X = size(imgDat, 2);

    dImgDat = double(imgDat);
    if ~isempty(mask)
        dImgDat(~mask) = NaN;
    end
    intHisto = isinteger(imgDat);

    statsStruct.stackStats = takeStats(dImgDat(:), intHisto);
    statsStruct.stackStats.datType = itype;

    if Z > 1
        statsStruct.sliceStats(Z) = genGeneralStatsStruct();
        for z = 1:Z
            sliceDat = dImgDat(:,:,z);
            sliceStats = takeStats(sliceDat(:), intHisto);
            sliceStats.datType = itype;
            statsStruct.sliceStats(z) = sliceStats;
        end

        maxproj = max(dImgDat, [], 3, 'omitnan');
        statsStruct.maxProjStats = takeStats(maxproj(:), intHisto);
        statsStruct.maxProjStats.datType = itype;
    end
end

function statsStruct = genGeneralStatsStruct()
    statsStruct = struct();

    statsStruct.datType = [];
    statsStruct.volume = 0;

    statsStruct.histo = [];
    statsStruct.min = NaN;
    statsStruct.max = NaN;
    statsStruct.mean = NaN;
    statsStruct.median = NaN;
    statsStruct.stdev = NaN;
    statsStruct.mad = NaN;
    statsStruct.ptiles = NaN(1, 99);
end

function statsStruct = takeStats(dataCollapsed, intHisto)
    statsStruct = genGeneralStatsStruct();

    if isempty(dataCollapsed)
        return;
    end

    notNan = ~isnan(dataCollapsed);
    if nnz(notNan) < 1; return; end
    dataCollapsed = dataCollapsed(notNan);

    statsStruct.volume = size(dataCollapsed, 2);

    statsStruct.min = min(dataCollapsed, [], 'all', 'omitnan');
    statsStruct.max = max(dataCollapsed, [], 'all', 'omitnan');
    statsStruct.mean = mean(dataCollapsed, 'all', 'omitnan');
    statsStruct.median = median(dataCollapsed, 'all', 'omitnan');
    statsStruct.stdev = std(dataCollapsed, 0, 'all', 'omitnan');
    statsStruct.mad = mad(dataCollapsed, 1, 'all');

    statsStruct.ptiles = prctile(dataCollapsed, 1:1:99, 'all');

    if intHisto
        E = 0:1:(statsStruct.max+1);
        [N, E] = histcounts(dataCollapsed, E);
    else
        [N, E] = histcounts(dataCollapsed);
    end

    binCount = size(N, 2);
    statsStruct.histo = NaN(binCount, 2);
    statsStruct.histo(1,:) = E(1:(binCount-1));
    statsStruct.histo(2,:) = N(:);
end