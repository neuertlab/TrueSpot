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
        end

        %%
        function nuc_lbl = runNuclearSegmentation(img_channel, param_struct)
            %Z Trim (if applicable)
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

            %Prep model
            cpObj = cellpose(Model=param_struct.nuc_params.model_name, ...
                UseEnsemble=param_struct.nuc_params.ensemble_bool);
            if ZZ > 1
                raw_lbl = segmentCells3D(cpObj, img_channel, ...
                    ImageCellDiameter=param_struct.nuc_params.avg_dia, ...
                    CellThreshold=param_struct.nuc_params.cell_threshold, ...
                    FlowErrorThreshold=param_struct.nuc_params.flow_threshold, ...
                    NormalizeInput=param_struct.nuc_params.normalize_bool, ...
                    VoxelSpacing=[1 1 param_struct.z2xy], ...
                    MinVolume=param_struct.nuc_params.min_size);
            else
                raw_lbl = segmentCells2D(cpObj, img_channel, ...
                    ImageCellDiameter=param_struct.nuc_params.avg_dia, ...
                    CellThreshold=param_struct.nuc_params.cell_threshold, ...
                    FlowErrorThreshold=param_struct.nuc_params.flow_threshold, ...
                    NormalizeInput=param_struct.nuc_params.normalize_bool);
            end

            nuc_lbl = uint16(raw_lbl);
        end

        %%
        function cell_lbl = runCytoSegmentation(light_channel, nuc_channel, param_struct)
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

            if isempty(nuc_channel)
                nuc_channel = zeros(size(light_channel));
            end

            %Prep model
            cpObj = cellpose(Model=param_struct.cyto_params.model_name, ...
                UseEnsemble=param_struct.cyto_params.ensemble_bool);
            if ZZ > 1
                raw_lbl = segmentCells3D(cpObj, light_channel, ...
                    AuxiliaryChannel=nuc_channel, ...
                    ImageCellDiameter=param_struct.cyto_params.avg_dia, ...
                    CellThreshold=param_struct.cyto_params.cell_threshold, ...
                    FlowErrorThreshold=param_struct.cyto_params.flow_threshold, ...
                    NormalizeInput=param_struct.cyto_params.normalize_bool, ...
                    VoxelSpacing=[1 1 param_struct.z2xy], ...
                    MinVolume=param_struct.cyto_params.min_size);
            else
                raw_lbl = segmentCells2D(cpObj, light_channel, ...
                    AuxiliaryChannel=nuc_channel, ...
                    ImageCellDiameter=param_struct.cyto_params.avg_dia, ...
                    CellThreshold=param_struct.cyto_params.cell_threshold, ...
                    FlowErrorThreshold=param_struct.cyto_params.flow_threshold, ...
                    NormalizeInput=param_struct.cyto_params.normalize_bool);
            end

            cell_lbl = uint16(raw_lbl);
        end

    end
end