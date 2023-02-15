%
%%
%For just taking the coord table and counting points by cell.

function Main_QuickCellQuant(coord_table_path, cell_seg_path, out_path, thresh_idx, ogdims)

addpath('./core');

%Load things
load(coord_table_path, 'coord_table');
load(cell_seg_path, 'cells');
cell_mask = cells;
clear cells;

%Check
if isempty(coord_table)
    fprintf('Provided coord table was empty! Exiting...\n');
    return;
end
if isempty(cell_mask)
    fprintf('Provided cell mask was empty! Exiting...\n');
    return;
end

th_idx = Force2Num(thresh_idx);
image_dims = parseDims(ogdims);

%Regurgiate input
fprintf('Main_QuickCellQuant (v 23.02.14) --\n');
fprintf('coord_table_path = %s\n', coord_table_path);
fprintf('cell_seg_path = %s\n', cell_seg_path);
fprintf('out_path = %s\n', out_path);
fprintf('image dims = %d x %d x %d\n', image_dims.x, image_dims.y, image_dims.z);
fprintf('thresh_idx = %d\n', th_idx);

%Allocate cell info
cell_count = max(cell_mask(:));
cell_info(cell_count) = SingleCell.newRNACell(cell_count, 0);
for i = 1:cell_count-1
    cell_info(i) = SingleCell.newRNACell(i, 0);
end

%Get requested coordinate table
if iscell(coord_table)
    coord_table = coord_table{th_idx,1};
end

%Per cell...
cell_spot_counts = zeros(cell_count,3);
for c = 1:cell_count
    cell_spot_counts(c,1) = c;
    cell_info(c).dim_z = image_dims.z;
    cell_info(c) = cell_info(c).findBoundaries(cell_mask, []);
    [cell_coord_table, nuc_coord_table] = cell_info(c).getCoordsSubset(coord_table);
    
    if ~isempty(cell_coord_table)
        cell_spot_counts(c,2) = size(cell_coord_table,1);
    end
    if ~isempty(nuc_coord_table)
        cell_spot_counts(c,3) = size(nuc_coord_table,1);
    end
end

%Save
if endsWith(out_path, '.csv')
    tbl_file = fopen(out_path, 'w');
    fprintf(tbl_file, 'CELL_IDX,SPOT_COUNT,NUC_SPOT_COUNT\n');
    for c = 1:cell_count
        fprintf(tbl_file, '%d,%d,%d\n',c,cell_spot_counts(c,2),cell_spot_counts(c,3));
    end
    fclose(tbl_file);
    fprintf('Saved comma delimited table to %s...\n', out_path);
elseif endsWith(out_path, '.tsv')
    tbl_file = fopen(out_path, 'w');
    fprintf(tbl_file, 'CELL_IDX\tSPOT_COUNT,NUC_SPOT_COUNT\n');
    for c = 1:cell_count
        fprintf(tbl_file, '%d\t%d\t%d\n',c,cell_spot_counts(c,2),cell_spot_counts(c,3));
    end
    fclose(tbl_file);
    fprintf('Saved tab delimited table to %s...\n', out_path);
elseif endsWith(out_path, '.mat')
    save(out_path, 'cell_spot_counts');
    fprintf('Saved MAT variable ''cell_spot_counts'' to %s...\n', out_path);
else
    fprintf('Output format not recognized. Exiting...\n');
    return;
end

end

function idims = parseDims(dims_arg)
    idims = struct('x',0,'y',0,'z',0);
	if ischar(dims_arg)
        %Read as string
        %'(x,y,z)'
        repl_str = replace(dims_arg, {'(', ')'}, '');
        split_str = split(repl_str, ',');
        ndims = size(split_str,1);
        idims.x = Force2Num(split_str{1,1});
        if ndims > 1; idims.y = Force2Num(split_str{2,1}); end
        if ndims > 2; idims.z = Force2Num(split_str{3,1}); end
	else
        %Try to read as vector
        if isvector(dims_arg)
            ndims = size(dims_arg,2);
            idims.x = dims_arg(1);
            if ndims > 1; idims.y = dims_arg(2); end
            if ndims > 2; idims.z = dims_arg(3); end
        end
	end
end