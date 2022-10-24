%
%%

classdef MatImages
    
    methods (Static)

        %------------ Read ------------

        function image_channel_set = loadImageChannels(filepath, varname)
            %Checks to see if it's just a raw matrix of ImageChannelSet.

            if nargin < 2
                varname = 'imgdat';
            end

            if isfile(filepath)
                %Assumed single mat
                image_channel_set = ImageChannelSet;
                load(filepath, varname);
                image_channel_set.dat_rna_sample = eval(varname);
            else
                %Assumed channel set
                image_channel_set = ImageChannelSet.loadMAT(filepath);
            end

        end

        %------------ Write ------------

        %------------ View ------------

        function fig_handle = viewMaxProjection(channel_data, bool_apply_filter, figno)

            if nargin < 3
                figno = round(rand() * 10000);
            end

            channel_data = uint16(channel_data);

            if bool_apply_filter
                channel_data = RNA_Threshold_Common.applyGaussianFilter(channel_data, 7, 2);
                channel_data = RNA_Threshold_Common.applyEdgeDetectFilter(channel_data);
                channel_data = RNA_Threshold_Common.blackoutBorders(channel_data, 7, 0);
                
                %Clean up borders
                chdbl = double(channel_data);
                border_mask = RNAUtils.genBorderMask(size(channel_data), [7 7 0]);
                borderidx = find(border_mask);
                chdbl(borderidx) = NaN;
                chmed = nanmedian(chdbl, 'all');
                chdbl(borderidx) = chmed;
                channel_data = uint16(chdbl);
            end
            
            max_proj = double(max(channel_data,[],3));
            Lmin = min(max_proj, [], 'all', 'omitnan');
            Lmax = median(max_proj(:)) + round(10 * std(max_proj(:)));

            fig_handle = figure(figno);
            imshow(max_proj, [Lmin Lmax]);
        end

    end

end