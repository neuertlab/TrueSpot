%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Load csv Table ==========================
addpath('./core');
InputTablePath = [DataDir filesep 'test_images.csv'];
image_table = testutil_opentable(InputTablePath);

out_table_path = [DataDir filesep 'idims.tsv'];

% ========================== Load csv Table ==========================

out_table = fopen(out_table_path, 'w');

rec_count = size(image_table,1);
for r = 124:243
    iname = getTableValue(image_table, r, 'IMGNAME');
    
    if endsWith(iname, 'rnasum')
        fprintf(out_table, '%s\t', iname);
        fprintf(out_table, '0\t0\t0\n');
        continue;
    end
    
    tifpath = [ImgDir replace(getTableValue(image_table, r, 'IMAGEPATH'), '/', filesep)];
    
    fprintf(out_table, '%s\t', iname);
    if isfile(tifpath)
        fprintf('%d of %d: "%s" found!\n', r, rec_count, tifpath);
        if endsWith(tifpath, '.tif')
            ch_tot = getTableValue(image_table, r, 'CH_TOTAL');
            ch_trg = getTableValue(image_table, r, 'CHANNEL');
        
            [channels, idims] = LoadTif(tifpath, ch_tot, [ch_trg], 1);
            clear channels;
            fprintf(out_table, '%d\t', idims.x);
            fprintf(out_table, '%d\t', idims.y);
            fprintf(out_table, '%d\n', idims.z);
        elseif endsWith(tifpath, '.mat')
            load(tifpath, 'imgdat');
            fprintf(out_table, '%d\t', size(imgdat,2));
            fprintf(out_table, '%d\t', size(imgdat,1));
            fprintf(out_table, '%d\n', size(imgdat,3));
            clear imgdat;
        else
            fprintf('%d of %d: "%s" not a valid image file. Skipping...\n', r, rec_count, tifpath);
            fprintf(out_table, '0\t0\t0\n');
        end
        
    else
        fprintf('%d of %d: image "%s" not found. Skipping...\n', r, rec_count, tifpath);
        fprintf(out_table, '0\t0\t0\n');
    end
    
end

fclose(out_table);

% ========================== Helper Functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end