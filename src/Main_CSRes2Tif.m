function Main_CSRes2Tif(matpath)

addpath('./core');
addpath('./thirdparty');

BUILD_STRING = '2024.11.14.00';
VERSION_STRING = 'v1.1.1';

fprintf('Running Main_CSRes2Tif\n');
fprintf('Build: %s\n', BUILD_STRING);
fprintf('TrueSpot Version: %s\n', VERSION_STRING);

if ~isfile(matpath)
    fprintf('ERROR: Input file %s does not exist! Exiting...\n', matpath);
    return;
end

[fdir, fname, ~] = fileparts(matpath);
nuc_out_path = [fdir filesep fname '_nuc.tif'];
cell_out_path = [fdir filesep fname '_cell.tif'];

fprintf('Loading cell mask...\n');
cell_mask = CellSeg.openCellMask(matpath);
if ~isempty(cell_mask)
    cell_mask = uint16(cell_mask);
    tifop = struct();
    tifop.overwrite = true;
    fprintf('Saving cell mask to %s...\n', cell_out_path);
    saveastiff(cell_mask, cell_out_path, tifop);
else
    fprintf('ERROR: Cell mask could not be loaded from input file! Skipping...\n');
end

fprintf('Loading nuc mask...\n');
nuc_mask = CellSeg.openNucMask(matpath);
if ~isempty(nuc_mask)
    nuc_mask = uint16(nuc_mask);

    if ~isempty(cell_mask)
        fprintf('Labeling nuclei...\n');
        Z = size(nuc_mask, 3);
        cm3 = repmat(cell_mask, [1 1 Z]);
        nuc_mask = immultiply(nuc_mask, cm3);
        nmcell_count = max(nuc_mask, [], 'all', 'omitnan');
        fprintf('Nuclei labeled: %d\n', nmcell_count);
    end

    tifop = struct();
    tifop.overwrite = true;
    fprintf('Saving nuc mask to %s...\n', nuc_out_path);
    saveastiff(nuc_mask, nuc_out_path, tifop);
else
    fprintf('ERROR: Nuc mask could not be loaded from input file! Skipping...\n');
end

fprintf('All done!\n');

end