%
%%
classdef CellPoseTS
    methods(Static)

        %%
        function param_struct = genCellposeSubParamStruct()
            param_struct = struct();
            param_struct.model_name = [];
            param_struct.avg_dia = NaN; %In pixels
            param_struct.cell_threshold = 0;
            param_struct.flow_threshold = 0.4;
            param_struct.min_size = 15;
            param_struct.max_size = NaN;
            param_struct.normalize_bool = false;
            param_struct.ensemble_bool = false;
            param_struct.do3D = true;
        end

        %%
        function param_struct = genCellposeParamStruct()
            param_struct = struct();
            param_struct.cyto_params = CellPoseTS.genCellposeSubParamStruct();
            param_struct.nuc_params = CellPoseTS.genCellposeSubParamStruct();

            param_struct.cyto_params.model_name = 'cyto'; %'cyto' or 'cyto2'
            param_struct.nuc_params.model_name = 'nuclei';
            param_struct.cyto_params.do3D = false;

            param_struct.z2xy = 1.0; %Z to XY ratio (For 3D)

            %Optional trimming
            param_struct.nzmin = NaN;
            param_struct.nzmax = NaN;
            param_struct.czmin = NaN;
            param_struct.czmax = NaN;
            param_struct.xtrim = 4;
            param_struct.ytrim = 4;
        end

        %%
        function [nuc_lbl, nuc_img] = runNuclearSegmentation(img_channel, param_struct, verbose)
            if nargin < 3; verbose = false; end

            %Z Trim (if applicable)
            if verbose; fprintf("[%s][CellPoseTS.runNuclearSegmentation] Applying z trim...\n", datetime); end
            Z = size(img_channel, 3);
            nzmin = 1;
            nzmax = Z;
            if ~isnan(param_struct.nzmin)
                if param_struct.nzmin > 1
                    nzmin = param_struct.nzmin;
                end
            end
            if ~isnan(param_struct.nzmax)
                if param_struct.nzmax < Z
                    nzmax = param_struct.nzmax;
                end
            end

            ZZ = Z;
            if Z > 1
                if param_struct.nuc_params.do3D
                    img_channel = img_channel(:,:,nzmin:nzmax);
                    ZZ = size(img_channel, 3);
                else
                    img_channel = max(img_channel(:,:,nzmin:nzmax), [], 3, 'omitnan');
                    ZZ = 1;
                end
            end

            nuc_img = img_channel;

            if verbose; fprintf("[%s][CellPoseTS.runNuclearSegmentation] Rescaling channel intensity...\n", datetime); end
            imin = min(img_channel, [], 'all', 'omitnan');
            img_channel = img_channel - imin;
            imax = max(img_channel, [], 'all', 'omitnan');
            img_channel = img_channel ./ imax;

            %Prep model
            if verbose; fprintf("[%s][CellPoseTS.runNuclearSegmentation] Prepping model (%s)...\n", datetime, param_struct.nuc_params.model_name); end
            cpObj = cellpose(Model=param_struct.nuc_params.model_name, ...
                UseEnsemble=param_struct.nuc_params.ensemble_bool);
            if ZZ > 1
                if verbose; fprintf("[%s][CellPoseTS.runNuclearSegmentation] Running 3D segmentation...\n", datetime); end
                raw_lbl = segmentCells3D(cpObj, img_channel, ...
                    ImageCellDiameter=param_struct.nuc_params.avg_dia, ...
                    CellThreshold=param_struct.nuc_params.cell_threshold, ...
                    NormalizeInput=param_struct.nuc_params.normalize_bool, ...
                    VoxelSpacing=[1 1 param_struct.z2xy], ...
                    MinVolume=param_struct.nuc_params.min_size);

               % FlowErrorThreshold=param_struct.nuc_params.flow_threshold, ...
            else
                if verbose; fprintf("[%s][CellPoseTS.runNuclearSegmentation] Running 2D segmentation...\n", datetime); end
                raw_lbl = segmentCells2D(cpObj, img_channel, ...
                    ImageCellDiameter=param_struct.nuc_params.avg_dia, ...
                    CellThreshold=param_struct.nuc_params.cell_threshold, ...
                    FlowErrorThreshold=param_struct.nuc_params.flow_threshold, ...
                    NormalizeInput=param_struct.nuc_params.normalize_bool);
            end

            %Additional cleanup
            if ~isnan(param_struct.nuc_params.max_size)
                fprintf("[%s][CellPoseTS.runNuclearSegmentation] Filtering out large objects...\n", datetime); 
                raw_lbl = CellSeg.removeObjectsAboveSize(raw_lbl, param_struct.nuc_params.max_size);
            end
            fprintf("[%s][CellPoseTS.runNuclearSegmentation] Trimming edges...\n", datetime); 
            if param_struct.xtrim > 0
                xtrim = param_struct.xtrim;
                X = size(raw_lbl, 2);
                raw_lbl(:, 1:xtrim, :) = 0;
                raw_lbl(:, (X-xtrim+1):X, :) = 0;
            end
            if param_struct.ytrim > 0
                ytrim = param_struct.ytrim;
                Y = size(raw_lbl, 1);
                raw_lbl(1:ytrim, :, :) = 0;
                raw_lbl((Y-ytrim+1):Y, :, :) = 0;
            end

            if verbose 
                fprintf("[%s][CellPoseTS.runNuclearSegmentation] Nuclei found: %d\n", datetime, max(raw_lbl, [], 'all', 'omitnan')); 
                fprintf("[%s][CellPoseTS.runNuclearSegmentation] Cleaning up...\n", datetime); 
            end
            nuc_lbl = uint16(raw_lbl);
        end

        %%
        function [cell_lbl, cell_img] = runCytoSegmentation(light_channel, nuc_channel, param_struct, verbose)
            if nargin < 4; verbose = false; end

            if verbose; fprintf("[%s][CellPoseTS.runCytoSegmentation] Applying z trim...\n", datetime); end

            %Z Trim (if applicable)
            Z = size(light_channel, 3);
            czmin = 1;
            czmax = Z;
            if ~isnan(param_struct.czmin)
                if param_struct.czmin > 1
                    czmin = param_struct.czmin;
                end
            end
            if ~isnan(param_struct.czmax)
                if param_struct.czmax < Z
                    czmax = param_struct.czmax;
                end
            end

            ZZ = Z;
            if Z > 1
                if param_struct.cyto_params.do3D
                    light_channel = light_channel(:,:,czmin:czmax);
                    if ~isempty(nuc_channel)
                        nuc_channel = nuc_channel(:,:,czmin:czmax);
                    end
                    ZZ = size(light_channel, 3);
                else
                    light_channel = max(light_channel(:,:,czmin:czmax), [], 3, 'omitnan');
                    if ~isempty(nuc_channel)
                        nuc_channel = max(nuc_channel(:,:,czmin:czmax), [], 3, 'omitnan');
                    end
                    ZZ = 1;
                end
            end

            cell_img = light_channel;

            if verbose; fprintf("[%s][CellPoseTS.runCytoSegmentation] Rescaling channel intensities...\n", datetime); end
            imin = min(light_channel, [], 'all', 'omitnan');
            light_channel = light_channel - imin;
            imax = max(light_channel, [], 'all', 'omitnan');
            light_channel = light_channel ./ imax;

            if isempty(nuc_channel)
                nuc_channel = zeros(size(light_channel));
            else
                imin = min(nuc_channel, [], 'all', 'omitnan');
                nuc_channel = nuc_channel - imin;
                imax = max(nuc_channel, [], 'all', 'omitnan');
                nuc_channel = nuc_channel ./ imax;
            end
            clear imin imax
            
            %Prep model
            if verbose; fprintf("[%s][CellPoseTS.runCytoSegmentation] Prepping model (%s)...\n", datetime, param_struct.cyto_params.model_name); end
            cpObj = cellpose(Model=param_struct.cyto_params.model_name, ...
                UseEnsemble=param_struct.cyto_params.ensemble_bool);
            if ZZ > 1
                if verbose; fprintf("[%s][CellPoseTS.runCytoSegmentation] Running 3D segmentation...\n", datetime); end
                raw_lbl = segmentCells3D(cpObj, light_channel, ...
                    AuxiliaryChannel=nuc_channel, ...
                    ImageCellDiameter=param_struct.cyto_params.avg_dia, ...
                    CellThreshold=param_struct.cyto_params.cell_threshold, ...
                    FlowErrorThreshold=param_struct.cyto_params.flow_threshold, ...
                    NormalizeInput=param_struct.cyto_params.normalize_bool, ...
                    VoxelSpacing=[1 1 param_struct.z2xy], ...
                    MinVolume=param_struct.cyto_params.min_size);
            else
                if verbose; fprintf("[%s][CellPoseTS.runCytoSegmentation] Running 2D segmentation...\n", datetime); end
                raw_lbl = segmentCells2D(cpObj, light_channel, ...
                    AuxiliaryChannel=nuc_channel, ...
                    ImageCellDiameter=param_struct.cyto_params.avg_dia, ...
                    CellThreshold=param_struct.cyto_params.cell_threshold, ...
                    FlowErrorThreshold=param_struct.cyto_params.flow_threshold, ...
                    NormalizeInput=param_struct.cyto_params.normalize_bool);
            end

            %Additional cleanup
            if ~isnan(param_struct.cyto_params.max_size)
                fprintf("[%s][CellPoseTS.runCytoSegmentation] Filtering out large objects...\n", datetime); 
                raw_lbl = CellSeg.removeObjectsAboveSize(raw_lbl, param_struct.cyto_params.max_size);
            end
            fprintf("[%s][CellPoseTS.runCytoSegmentation] Trimming edges...\n", datetime); 
            if param_struct.xtrim > 0
                xtrim = param_struct.xtrim;
                X = size(raw_lbl, 2);
                raw_lbl(:, 1:xtrim, :) = 0;
                raw_lbl(:, (X-xtrim+1):X, :) = 0;
            end
            if param_struct.ytrim > 0
                ytrim = param_struct.ytrim;
                Y = size(raw_lbl, 1);
                raw_lbl(1:ytrim, :, :) = 0;
                raw_lbl((Y-ytrim+1):Y, :, :) = 0;
            end

            if verbose 
                fprintf("[%s][CellPoseTS.runCytoSegmentation] Cells found: %d\n", datetime, max(raw_lbl, [], 'all', 'omitnan')); 
                fprintf("[%s][CellPoseTS.runCytoSegmentation] Cleaning up...\n", datetime); 
            end
            cell_lbl = uint16(raw_lbl);
        end


        %%
        function [nuc_lbl, nuc_stats, bkg_stats] = runHybrid3DNuclearSegmentation(nuc_channel, param_struct, verbose)
            %Use cellpose to segment max projection, then use this info to
            %expand to 3D
            if verbose; fprintf("[%s][CellPoseTS.runHybrid3DNuclearSegmentation] Detecting nuclear labels from max proj...\n", datetime); end
            param_struct.nuc_params.do3D = false;
            [nuc_lbl2, nuc_img] = CellPoseTS.runNuclearSegmentation(nuc_channel, param_struct, verbose);

            Z = size(nuc_channel, 3);
            if Z < 2; return; end
            if verbose; fprintf("[%s][CellPoseTS.runHybrid3DNuclearSegmentation] Refining nuclear mask...\n", datetime); end
            nuc_channel = double(nuc_channel);
            nucCount = max(nuc_lbl2, [], 'all', 'omitnan');
            nuc_lbl = zeros(size(nuc_channel));
            if nucCount < 1; return; end

            %DEBUG
            figure(615);
            imshow(nuc_lbl2, []);

            if verbose; fprintf("[%s][CellPoseTS.runHybrid3DNuclearSegmentation] Calculating background stats...\n", datetime); end
            bkg_stats = struct();
            nuc_mask2 = (nuc_lbl2 == 0);
            idat = nuc_img;
            idat(~nuc_mask2) = NaN;
            bkg_stats.mean = mean(idat, 'all', 'omitnan');
            bkg_stats.stdev = std(idat, 0, 'all', 'omitnan');
            bkg_stats.median = median(idat, 'all', 'omitnan');
            bkg_stats.init_th = (bkg_stats.median + bkg_stats.stdev);

            nuc_stats(nucCount) = bkg_stats;

            for n = 1:nucCount
                if verbose; fprintf("[%s][CellPoseTS.runHybrid3DNuclearSegmentation] Working on nucleus %d...\n", datetime, n); end
                nstats = nuc_stats(n);
                nuc_mask2 = (nuc_lbl2 == n);
                idat = nuc_img;
                idat(~nuc_mask2) = NaN;
                nuc_mask2_dbl = double(nuc_mask2);

                nstats.mean = mean(idat, 'all', 'omitnan');
                nstats.stdev = std(idat, 0, 'all', 'omitnan');
                nstats.median = median(idat, 'all', 'omitnan');
                nstats.init_th = (nstats.median - (nstats.stdev * 1.5));

                for z = 1:Z
                    slice_filt = immultiply(nuc_mask2_dbl, nuc_channel(:,:,z) >= nstats.init_th);
                    slice_filt = imclose(slice_filt, strel('disk', 5));
                    slice_filt = bwareaopen(slice_filt, 15); %Filter out little blips
                    slice_filt = immultiply(nuc_mask2_dbl, slice_filt);
                    nuc_lbl(:,:,z) = nuc_lbl(:,:,z) + (slice_filt .* double(n));
                end
                clear slice_filt

                %TODO Fill in to get same max proj?

                nuc_stats(n) = nstats;
            end

            nuc_lbl = uint16(nuc_lbl);

            %DEBUG
            % for z = 1:Z
            %     figure(z);
            %     imshow(nuc_lbl(:,:,z), []);
            % end
            figure(100);
            imshow(max(nuc_lbl, [], 3, 'omitnan'), []);
        end

        %%
        function [cell_lbl, nuc_lbl] = matchCellNucLabels(cell_lbl, nuc_lbl, nuc_bkg_stats)
            %Make sure nuc label is consistant with cell label.
            Znuc = size(nuc_lbl, 3);
            Zcell = size(cell_lbl, 3);

            nuc_lbl = double(nuc_lbl ~= 0);
            if Znuc == Zcell
                cell_lbl_temp = double(cell_lbl);
                nuc_lbl = immultiply(nuc_lbl, cell_lbl_temp);
                nuc_lbl = uint16(nuc_lbl);
            else
                %Collapse cell mask if not already 2D and do nuc mask by
                %slice
                if Zcell > 1
                    cell_lbl_2 = double(max(cell_lbl, [], 3, 'omitnan'));
                else
                    cell_lbl_2 = double(cell_lbl);
                end

                for z = 1:Znuc
                    nuc_lbl(:,:,z) = immultiply(cell_lbl_2, nuc_lbl(:,:,z));
                end
                nuc_lbl = uint16(nuc_lbl);
            end

            %TODO Retry finding nuc for cells without nuclei

            %DEBUG
            figure(200);
            imshow(max(nuc_lbl, [],3, 'omitnan'), []);

            figure(201);
            imshow(cell_lbl, []);

        end

    end
end