%
%%
classdef CellsegDrawer
    
    %%
    properties
        cell_color;
        cell_alpha;
        nuc_color;
        nuc_alpha;
        
        cell_mask;
        nuc_mask;

        useMaxProj;
    end
    
    %%
    methods

        %%
        function obj = initializeMe(obj)
            obj.cell_color = [0.8 0.0 0.8];
            obj.cell_alpha = 0.25;
            obj.nuc_color = [0.369 0.510 0.788];
            obj.nuc_alpha = 0.25;
            obj.cell_mask = [];
            obj.nuc_mask = [];
            obj.useMaxProj = true;
        end

        %%
        function imgOut = preprocessInputImage(obj, imgIn, z)

            %If stack, collapse. But don't collapse if RGB.
            %Also need to adjust range.

            imgOut = imgIn;
            if ndims(imgIn) > 2
                %See if it's already RGB
                if isa(imgIn, 'uint8') & (size(imgIn, 3) == 3)
                    return;
                end

                if obj.useMaxProj | (z < 1)
                    imgOut = max(imgIn, [], 3);
                else
                    imgOut = imgIn(:,:,z);
                end
            end

%             if isa(imgIn, 'double')
%                 inMax = max(imgIn, [], 'all', 'omitnan');
%                 if inMax <= 1.0
%                     return;
%                 end
%             end
% 
%             imgOut = uint16(imgOut);
%             imgOut = imadjust(imgOut);
        end

        %%
        function imgOut = applyCellMask(obj, imgIn, z)
            if nargin < 3; z = 0; end
            if isempty(imgIn)
                imgOut = [];
                return;
            end

            imgIn = obj.preprocessInputImage(imgIn, z);

            if isempty(obj.cell_mask)
                imgOut = imgIn;
                return;
            end

            mm = obj.cell_mask;
            if ndims(obj.cell_mask) > 2
                if obj.useMaxProj | (z < 1)
                    mm = max(obj.cell_mask, [], 3);
                else
                    mm = obj.cell_mask(:,:,z);
                end
            end

            imgOut = labeloverlay(imgIn, mm, 'Colormap', obj.cell_color, 'Transparency', 1.0 - obj.cell_alpha);
        end

        %%
        function imgOut = applyNucMask(obj, imgIn, z)
            if nargin < 3; z = 0; end
            if isempty(imgIn)
                imgOut = [];
                return;
            end

            imgIn = obj.preprocessInputImage(imgIn, z);

            if isempty(obj.nuc_mask)
                imgOut = imgIn;
                return;
            end

            mm = obj.nuc_mask;
            if ndims(obj.nuc_mask) > 2
                if obj.useMaxProj | (z < 1)
                    mm = max(obj.nuc_mask, [], 3);
                else
                    mm = obj.nuc_mask(:,:,z);
                end
            end

            imgOut = labeloverlay(imgIn, mm, 'Colormap', obj.nuc_color, 'Transparency', 1.0 - obj.nuc_alpha);
        end

    end
    
    %%
    methods (Static)
        
        %% ========================== Structs ==========================
        
        
        %% ========================== Rendering ==========================
        
        
        %% ========================== Drawing ==========================
        
        
        %% ========================== Interface ==========================
        
        
        
    end
end