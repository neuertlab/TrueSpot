%
%%

classdef ImageChannelSet
    
    properties
        dat_rna_sample
        dat_trans_sample
        dat_rna_control
        dat_trans_control
    end

    methods(Static)

        function img_ch_set = loadMAT(mat_stem)
            img_ch_set = ImageChannelSet;

            mat_path = [mat_stem '_datSample.mat'];
            if isfile(mat_path)
                load(mat_path, 'imgdat');
                img_ch_set.dat_rna_sample = imgdat;
            end

            mat_path = [mat_stem '_transSample.mat'];
            if isfile(mat_path)
                load(mat_path, 'imgdat');
                img_ch_set.dat_trans_sample = imgdat;
            end

            mat_path = [mat_stem '_datControl.mat'];
            if isfile(mat_path)
                load(mat_path, 'imgdat');
                img_ch_set.dat_rna_control = imgdat;
            end

            mat_path = [mat_stem '_transControl.mat'];
            if isfile(mat_path)
                load(mat_path, 'imgdat');
                img_ch_set.dat_trans_control = imgdat;
            end

        end

        function saveMAT(img_ch_set, mat_stem)

            if isempty(img_ch_set); return; end

            if ~isempty(img_ch_set.dat_rna_sample)
                mat_path = [mat_stem '_datSample.mat'];
                imgdat = img_ch_set.dat_rna_sample;
                save(mat_path, 'imgdat');
            end

            if ~isempty(img_ch_set.dat_trans_sample)
                mat_path = [mat_stem '_transSample.mat'];
                imgdat = img_ch_set.dat_trans_sample;
                save(mat_path, 'imgdat');
            end

            if ~isempty(img_ch_set.dat_rna_control)
                mat_path = [mat_stem '_datControl.mat'];
                imgdat = img_ch_set.dat_rna_control;
                save(mat_path, 'imgdat');
            end

            if ~isempty(img_ch_set.dat_trans_control)
                mat_path = [mat_stem '_transControl.mat'];
                imgdat = img_ch_set.dat_trans_control;
                save(mat_path, 'imgdat');
            end

        end

    end
    
end