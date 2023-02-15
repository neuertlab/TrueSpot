%
%%  !! UPDATE TO YOUR BASE DIR
DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Settings ==========================
addpath('./core');

ImgName = 'scrna_E1R2I4_CTT1';

OutDir = [DataDir filesep 'data' filesep 'deepblink_training'];
ApplyMask = true;
mode = 'maxproj'; %'3d' 'splitz' 'maxproj'

% ========================== Read Table ==========================

InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

rec_row = 0;
rec_count = size(image_table,1);
for r = 1:rec_count
    iname = getTableValue(image_table, r, 'IMGNAME');
    if strcmp(iname, ImgName)
        rec_row = r;
        break;
    end
end

if rec_row < 1
    fprintf('Image with name %s could not be found!\n', ImgName);
    return;
end

% ========================== Load Spots Anno ==========================

iname = getTableValue(image_table, rec_row, 'IMGNAME');
outstem = [DataDir replace(getTableValue(image_table, rec_row, 'OUTSTEM'), '/', filesep)];

spotsrun = RNASpotsRun.loadFrom(outstem);
if isempty(spotsrun)
    fprintf('Run with prefix "%s" was not found!\n', outstem);
    return;
end

%Update spotsrun paths
spotsrun.out_stem = outstem;
spotsrun.saveMe();

if ~RNA_Threshold_SpotSelector.selectorExists(outstem)
    fprintf('Spot selector for "%s" was not found!\n', outstem);
    return;
end

selector = RNA_Threshold_SpotSelector.openSelector(outstem, true);

%Calculate base coords if masked
tifpath = [ImgDir replace(getTableValue(image_table, rec_row, 'IMAGEPATH'), '/', filesep)];
nchan = getTableValue(image_table, rec_row, 'CH_TOTAL');
tchan = getTableValue(image_table, rec_row, 'CHANNEL');

[channels, ~] = LoadTif(tifpath, nchan, [tchan], 1);
my_image = channels{tchan,1};
my_image = uint16(my_image);
clear channels;

x0 = 1; x1 = size(my_image,2);
y0 = 1; y1 = size(my_image,1);
z0 = 1; z1 = size(my_image,3);

if ApplyMask
    z0 = selector.z_min; z1 = selector.z_max;
    if ~isempty(selector.selmcoords)
        x0 = selector.selmcoords(1,1);
        x1 = selector.selmcoords(2,1);
        y0 = selector.selmcoords(3,1);
        y1 = selector.selmcoords(4,1);
    end
end

%Okay, now make sure dimensions are powers of 2 I guess.
x_dim = dimto2pow(x1 - x0 + 1); x1 = x0 + x_dim - 1;
y_dim = dimto2pow(y1 - y0 + 1); y1 = y0 + y_dim - 1;
z_dim = z1 - z0 + 1;

if strcmp(mode, '3d')
    z_dim = dimto2pow(z1 - z0 + 1);
    flip = false;
    while (z1 - z0 + 1) > z_dim
        if ~flip
            z0 = z0 + 1;
        else
            z1 = z1 - 1;
        end
        flip = ~flip;
    end
    z_dim = z1 - z0 + 1;
end

%Output the ref set tables
if isempty(selector.ref_coords)
    fprintf('Spot selector for "%s" does not have reference set!\n', outstem);
    return;
end

%Determine which z slices to use by counting true spots on each slice
refspots = size(selector.ref_coords, 1);
fprintf('Trim dimensions: %d x %d x %d\n', x_dim, y_dim, z_dim);

if strcmp(mode, '3d')
    RefDumpPath = [OutDir filesep iname '.csv'];
    coords_csv = fopen(RefDumpPath, 'w');
    fprintf(coords_csv, ' ,X,Y,Slice\n');
    j = 1;
    for i = 1:refspots
        %Make sure coords fall within mask...
        if ApplyMask
            if selector.ref_coords(i,1) < x0; continue; end
            if selector.ref_coords(i,1) > x1; continue; end
            if selector.ref_coords(i,2) < y0; continue; end
            if selector.ref_coords(i,2) > y1; continue; end
            if selector.ref_coords(i,3) < z0; continue; end
            if selector.ref_coords(i,3) > z1; continue; end
        end
        
        fprintf(coords_csv, '%d,', j);
        fprintf(coords_csv, '%f,', selector.ref_coords(i,1) - x0); %x and y are zero based?? But not z??
        fprintf(coords_csv, '%f,', selector.ref_coords(i,2) - y0);
        fprintf(coords_csv, '%d\n', selector.ref_coords(i,3) - z0 + 1);
        
        j = j+1;
    end
    fclose(coords_csv);
    
    if ApplyMask
        export_tif_path = [OutDir filesep iname '.tif'];
        
        fprintf('Exporting masked image: (%d,%d,%d)(%d,%d,%d)\n', x0, y0, z0, x1, y1, z1);
        my_image_crop = my_image(y0:y1,x0:x1,z0:z1);
        
        tifops.overwrite = true;
        tifops.color = false;
        saveastiff(my_image_crop, export_tif_path, tifops);
    end
    clear my_image;
    
elseif strcmp(mode, 'splitz')
    z_counts = zeros(1,size(my_image,3));
    for i = 1:refspots
        zval = selector.ref_coords(i,3);
        if ApplyMask
            if zval < z0; continue; end
            if zval > z1; continue; end
        end
        z_counts(zval) = z_counts(zval) + 1;
    end
    
    for z = z0:z1
        if z_counts(z) < 1; continue; end
        
        RefDumpPath = [OutDir filesep iname '_z' num2str(z) '.csv'];
        coords_csv = fopen(RefDumpPath, 'w');
        fprintf(coords_csv, ' ,X,Y,Slice\n');
        j = 1;
        for i = 1:refspots
            %Make sure coords fall within mask...
            if ApplyMask
                if selector.ref_coords(i,1) < x0; continue; end
                if selector.ref_coords(i,1) > x1; continue; end
                if selector.ref_coords(i,2) < y0; continue; end
                if selector.ref_coords(i,2) > y1; continue; end
                if selector.ref_coords(i,3) ~= z; continue; end
            end
            
            fprintf(coords_csv, '%d,', j);
            fprintf(coords_csv, '%f,', selector.ref_coords(i,1) - x0); %x and y are zero based?? But not z??
            fprintf(coords_csv, '%f,', selector.ref_coords(i,2) - y0);
            %fprintf(coords_csv, '%d\n', selector.ref_coords(i,3) - z0 + 1);
            fprintf(coords_csv, '1\n');
            
            j = j+1;
        end
        fclose(coords_csv);
        
        %Dump mask trimmed TIFs, if requested
        if ApplyMask
            export_tif_path = [OutDir filesep iname '_z' num2str(z) '.tif'];
            
            %fprintf('Exporting masked image: (%d,%d,%d)(%d,%d,%d)\n', x0, y0, z0, x1, y1, z1);
            my_image_z = my_image(y0:y1,x0:x1,z);
            
            tifops.overwrite = true;
            tifops.color = false;
            saveastiff(my_image_z, export_tif_path, tifops);
        end
    end
    clear my_image;
elseif strcmp(mode, 'maxproj')
    RefDumpPath = [OutDir filesep iname '_maxproj.csv'];
    coords_csv = fopen(RefDumpPath, 'w');
    fprintf(coords_csv, ' ,X,Y,Slice\n');
    j = 1;
    for i = 1:refspots
        %Make sure coords fall within mask...
        if ApplyMask
            if selector.ref_coords(i,1) < x0; continue; end
            if selector.ref_coords(i,1) > x1; continue; end
            if selector.ref_coords(i,2) < y0; continue; end
            if selector.ref_coords(i,2) > y1; continue; end
            if selector.ref_coords(i,3) < z0; continue; end
            if selector.ref_coords(i,3) > z1; continue; end
        end
        
        fprintf(coords_csv, '%d,', j);
        fprintf(coords_csv, '%f,', selector.ref_coords(i,1) - x0); %x and y are zero based?? But not z??
        fprintf(coords_csv, '%f,', selector.ref_coords(i,2) - y0);
        fprintf(coords_csv, '1\n');
        
        j = j+1;
    end
    fclose(coords_csv);
    
    if ApplyMask
        export_tif_path = [OutDir filesep iname '_maxproj.tif'];
        
        %fprintf('Exporting masked image: (%d,%d,%d)(%d,%d,%d)\n', x0, y0, z0, x1, y1, z1);
        my_image_crop = my_image(:,:,z0:z1);
        my_image_crop = max(my_image_crop,[],3);
        my_image_crop = my_image_crop(y0:y1,x0:x1);
        
        tifops.overwrite = true;
        tifops.color = false;
        saveastiff(my_image_crop, export_tif_path, tifops);
    end
end


% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end

function output_dim = dimto2pow(input_dim)
    output_dim = 1;
    nextpow = 2;

    while nextpow <= input_dim
        output_dim = nextpow;
        nextpow = nextpow * 2;
    end

end
