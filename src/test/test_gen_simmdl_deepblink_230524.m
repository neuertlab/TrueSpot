%
%%  !! UPDATE TO YOUR BASE DIR
ImgDir = 'C:\Users\hospelb\labdata\imgproc';
%ImgDir = 'D:\usr\bghos\labdat\imgproc';

DataDir = 'D:\Users\hospelb\labdata\imgproc\imgproc';
%DataDir = 'D:\usr\bghos\labdat\imgproc';

% ========================== Constants ==========================

ModelName = 'SimRandoModel';
DeepBlinkModelDir = [DataDir filesep 'data' filesep 'deepblink_training'];
ResultsDir = [DataDir filesep 'data' filesep 'results'];

ZEROPROP_MAX = 0.7;
COORDS_ONE_BASED = false;

% ========================== Prep ==========================
InputTablePath_A = [DataDir filesep 'test_images_simytc.csv'];
InputTablePath_B = [DataDir filesep 'test_images_simvarmass.csv'];

image_table_A = testutil_opentable(InputTablePath_A);
image_table_B = testutil_opentable(InputTablePath_B);

ModelDir = [DeepBlinkModelDir filesep ModelName];
LabelDir = [ModelDir filesep 'labels'];
PredDir = [ModelDir filesep 'predict'];
if ~isfolder(ModelDir)
    mkdir(ModelDir);
    mkdir(LabelDir);
    mkdir(PredDir);
end

ListPathTraining = [DeepBlinkModelDir filesep 'simrando_train.txt'];
ListPathPredict = [DeepBlinkModelDir filesep 'simrando_pred.txt'];

ListFileTraining = fopen(ListPathTraining, 'w');
ListFilePredict = fopen(ListPathPredict, 'w');

% ========================== Sim YTC Table ==========================

train_count_cy5l = 0;
train_count_tmrl = 0;
pred_count_cy5l = 0;
pred_count_tmrl = 0;

entry_count = size(image_table_A, 1);
for r = 1:entry_count

    myname = getTableValue(image_table_A, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

    pass_zeroprop = true;
    if ZEROPROP_MAX > 0.0
        pass_zeroprop = false;
        ResFilePath = [ResultsDir filesep 'simytc' filesep myname '_summary.mat'];

        if isfile(ResFilePath)
            load(ResFilePath, 'analysis');
            if isfield(analysis, 'results_hb')
                nzprop = analysis.results_hb.fprop_nz;
                zprop = 1.0 - nzprop;
                if zprop <= ZEROPROP_MAX
                    pass_zeroprop = true;
                end
            end
            clear analysis;
        end
    end

    if pass_zeroprop
        %decide if train or predict...
        astrain = false;
        if contains(myname, '_CY5L_')
            if train_count_cy5l < 100
                if train_count_cy5l < pred_count_cy5l
                    astrain = true;
                end
            end

            if astrain
                train_count_cy5l = train_count_cy5l + 1; 
            else
                pred_count_cy5l = pred_count_cy5l + 1;
            end
        else
            if train_count_tmrl < 100
                if train_count_tmrl < pred_count_tmrl
                    astrain = true;
                end
            end

            if astrain
                train_count_tmrl = train_count_tmrl + 1; 
            else
                pred_count_tmrl = pred_count_tmrl + 1;
            end
        end

        tif_path_raw = getTableValue(image_table_A, r, 'IMAGEPATH');
        tif_path = [ImgDir replace(tif_path_raw, '/', filesep)];
        if astrain
            %Find and load the source mat file.
            [tifdir, fn, ~] = fileparts(tif_path);
            [matdir, ~, ~] = fileparts(tifdir);
            matpath = [matdir filesep fn '.mat'];
            load(matpath, 'imgdat', 'key');

            key_len = size(key,2);
            key_mtx = NaN(key_len,3);
            key_mtx(:,1) = [key.x];
            key_mtx(:,2) = [key.y];
            key_mtx(:,3) = [key.z];
            
            %We will have to split by slice...
            Z = size(imgdat, 3);
            for z = 1:Z
                znum = sprintf('%02d', z);
                slice_tif_path = [ModelDir filesep myname '_' znum '.tif'];
                slice_key_path = [LabelDir filesep myname '_' znum '.csv'];

                %Filter key hits that are on this slice.
                keyz_rows = find(key_mtx(:,3) == z);
                if ~isempty(keyz_rows)
                    key_x = key_mtx(keyz_rows,1);
                    key_y = key_mtx(keyz_rows,2);

                    if ~COORDS_ONE_BASED
                        key_x = key_x - 1;
                        key_y = key_y - 1;
                    end

                    slice_key_file = fopen(slice_key_path, 'w');
                    fprintf(slice_key_file, ' ,X,Y,Slice\n');

                    key_ct = size(key_x,1);
                    for j = 1:key_ct
                        fprintf(slice_key_file, '%d,%f,%f,', j, key_x(j,1), key_y(j,1));
                        if COORDS_ONE_BASED
                            fprintf(slice_key_file, '1\n');
                        else
                            fprintf(slice_key_file, '0\n');
                        end
                    end

                    fclose(slice_key_file);
                else
                    %No hits on this slice, so skip it.
                    continue;
                end

                my_slice = imgdat(:,:,z);
                saveastiff(my_slice, slice_tif_path);
            end

            clear matpath imgdat key
            fprintf(ListFileTraining, '%s\n', myname);
        else
            %Just copy the tif to the pred dir
            copyfile(tif_path, [PredDir filesep myname '.tif']);
            fprintf(ListFilePredict, '%s\n', myname);
        end
    end
end

% ========================== Sim Mass Table ==========================

train_count_mass = 0;
pred_count_mass = 0;

entry_count = size(image_table_B, 1);
for r = 1:entry_count

    myname = getTableValue(image_table_B, r, 'IMGNAME');
    fprintf('> Now processing %s (%d of %d)...\n', myname, r, entry_count);

    pass_zeroprop = true;
    if ZEROPROP_MAX > 0.0
        pass_zeroprop = false;
        ResFilePath = [ResultsDir filesep 'simvarmass' filesep myname '_summary.mat'];

        if isfile(ResFilePath)
            load(ResFilePath, 'analysis');
            if isfield(analysis, 'results_hb')
                nzprop = analysis.results_hb.fprop_nz;
                zprop = 1.0 - nzprop;
                if zprop <= ZEROPROP_MAX
                    pass_zeroprop = true;
                end
            end
            clear analysis;
        end
    end

    if pass_zeroprop
        %decide if train or predict...
        astrain = false;
        if train_count_mass < 100
            if train_count_mass < pred_count_mass
                astrain = true;
            end
        end

        if astrain
            train_count_mass = train_count_mass + 1;
        else
            pred_count_mass = pred_count_mass + 1;
        end

        tif_path_raw = getTableValue(image_table_B, r, 'IMAGEPATH');
        tif_path = [ImgDir replace(tif_path_raw, '/', filesep)];
        if astrain
            %Find and load the source mat file.
            [tifdir, fn, ~] = fileparts(tif_path);
            [matdir, ~, ~] = fileparts(tifdir);
            matpath = [matdir filesep fn '.mat'];
            load(matpath, 'imgdat', 'key');

            key_len = size(key,2);
            key_mtx = NaN(key_len,3);
            key_mtx(:,1) = [key.x];
            key_mtx(:,2) = [key.y];
            key_mtx(:,3) = [key.z];
            
            %We will have to split by slice...
            Z = size(imgdat, 3);
            for z = 1:Z
                znum = sprintf('%02d', z);
                slice_tif_path = [ModelDir filesep myname '_' znum '.tif'];
                slice_key_path = [LabelDir filesep myname '_' znum '.csv'];

                %Filter key hits that are on this slice.
                keyz_rows = find(key_mtx(:,3) == z);
                if ~isempty(keyz_rows)
                    key_x = key_mtx(keyz_rows,1);
                    key_y = key_mtx(keyz_rows,2);

                    if ~COORDS_ONE_BASED
                        key_x = key_x - 1;
                        key_y = key_y - 1;
                    end

                    slice_key_file = fopen(slice_key_path, 'w');
                    fprintf(slice_key_file, ' ,X,Y,Slice\n');

                    key_ct = size(key_x,1);
                    for j = 1:key_ct
                        fprintf(slice_key_file, '%d,%f,%f,', j, key_x(j,1), key_y(j,1));
                        if COORDS_ONE_BASED
                            fprintf(slice_key_file, '1\n');
                        else
                            fprintf(slice_key_file, '0\n');
                        end
                    end

                    fclose(slice_key_file);
                else
                    %No hits on this slice, so skip it.
                    continue;
                end

                my_slice = imgdat(:,:,z);
                saveastiff(my_slice, slice_tif_path);
            end

            clear matpath imgdat key
            fprintf(ListFileTraining, '%s\n', myname);
        else
            %Just copy the tif to the pred dir
            copyfile(tif_path, [PredDir filesep myname '.tif']);
            fprintf(ListFilePredict, '%s\n', myname);
        end
    end
end

% ========================== Cleanup ==========================

fclose(ListFileTraining);
fclose(ListFilePredict);

% ========================== Helper functions ==========================

function val = getTableValue(mytable, row_index, field)
    val = mytable{row_index,field};
    if iscell(val)
        val = val{1,1};
    end
end