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
        drawCellMaskAsOutline = true;
        drawNucMaskAsOutline = true;
        cellOutlineDiskSize = 3;
        nucOutlineDiskSize = 3;

        greyscaleLUT = [];
    end
    
    %%
    methods

        %%
        function obj = initializeMe(obj)
            obj.cell_color = [0.8 0.0 0.8];
            obj.cell_alpha = 0.5;
            obj.nuc_color = [0.369 0.510 0.788];
            obj.nuc_alpha = 0.5;
            obj.cell_mask = [];
            obj.nuc_mask = [];
            obj.useMaxProj = true;
            obj.greyscaleLUT = VisCommon.genGreyscaleLUT();
        end

        %%
        function [imgOut, xyView] = preprocessInputImage(obj, imgIn, z, xyView)
            if nargin < 4; xyView = []; end

            %If stack, collapse. But don't collapse if RGB.
            %Also need to adjust range.

            imgOut = imgIn;
            %XYTrim
            if ~isempty(xyView)
                x0 = 1; x1 = size(imgIn, 2);
                y0 = 1; y1 = size(imgIn, 1);
                if ~isnan(xyView.x0); x0 = xyView.x0; end
                if ~isnan(xyView.x1); x1 = xyView.x1; end
                if ~isnan(xyView.y0); y0 = xyView.y0; end
                if ~isnan(xyView.y1); y1 = xyView.y1; end

                imgOut = imgOut(y0:y1,x0:x1,:);
                xyView = struct('x0', x0, 'x1', x1, 'y0', y0, 'y1', y1);
            else
                xyView = struct('x0', 1, 'x1', size(imgIn, 2), 'y0', 1, 'y1', size(imgIn, 1));
            end

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

            imgOut = VisCommon.bw2rgb(imgOut, obj.greyscaleLUT, false);

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
        function imgOut = applyCellMask(obj, imgIn, z, xyView)
            if nargin < 3; z = 0; end
            if nargin < 4; xyView = []; end
            if isempty(imgIn)
                imgOut = [];
                return;
            end

            [imgIn, xyView] = obj.preprocessInputImage(imgIn, z, xyView);

            if isempty(obj.cell_mask)
                imgOut = imgIn;
                return;
            end

            x0 = xyView.x0;
            x1 = xyView.x1;
            y0 = xyView.y0;
            y1 = xyView.y1;

            mm = obj.cell_mask;
            if ndims(obj.cell_mask) > 2
                if obj.useMaxProj | (z < 1)
                    mm = max(obj.cell_mask(y0:y1,x0:x1,:), [], 3);
                else
                    mm = obj.cell_mask(y0:y1,x0:x1,z);
                end
            end

            %imgOut = labeloverlay(imgIn, mm, 'Colormap', obj.cell_color, 'Transparency', 1.0 - obj.cell_alpha);
            imgOut = VisCommon.compositeMaskOverlay(imgIn, mm, obj.cell_color, obj.cell_alpha, obj.drawCellMaskAsOutline, obj.cellOutlineDiskSize);
        end

        %%
        function imgOut = applyNucMask(obj, imgIn, z, xyView)
            if nargin < 3; z = 0; end
            if nargin < 4; xyView = []; end
            if isempty(imgIn)
                imgOut = [];
                return;
            end

            [imgIn, xyView] = obj.preprocessInputImage(imgIn, z, xyView);

            if isempty(obj.nuc_mask)
                imgOut = imgIn;
                return;
            end

            x0 = xyView.x0;
            x1 = xyView.x1;
            y0 = xyView.y0;
            y1 = xyView.y1;

            mm = obj.nuc_mask;
            if ndims(obj.nuc_mask) > 2
                if obj.useMaxProj | (z < 1)
                    mm = max(obj.nuc_mask(y0:y1,x0:x1,:), [], 3);
                else
                    mm = obj.nuc_mask(y0:y1,x0:x1,z);
                end
            end

            %imgOut = labeloverlay(imgIn, mm, 'Colormap', obj.nuc_color, 'Transparency', 1.0 - obj.nuc_alpha);
            imgOut = VisCommon.compositeMaskOverlay(imgIn, mm, obj.nuc_color, obj.nuc_alpha, obj.drawNucMaskAsOutline, obj.nucOutlineDiskSize);
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