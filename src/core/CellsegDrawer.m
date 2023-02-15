%
%%
classdef CellsegDrawer
    
    %%
    properties
        ch_light_proj;
        ch_nuc_proj;
        
        cell_mask;
        nuc_mask;
    end
    
    %%
    methods
    end
    
    %%
    methods (Static)
        
        %% ========================== Structs ==========================
        
        %%
        function render_struct = generateRenderer(X, Y, Z)
            render_struct = struct('dimX', X);
            render_struct.dimY = Y;
            render_struct.dimZ = Z;
            render_struct.r = [];
            render_struct.g = [];
            render_struct.b = [];
            render_struct.r_draw = [];
            render_struct.g_draw = [];
            render_struct.b_draw = [];
            render_struct.cell_pos = []; %For labeling cells. Cols: center_x, center_y, z_min, z_max
        end
        
        %% ========================== Rendering ==========================
        
        %%
        function render_struct = renderCellMask(cell_mask, cell_color, bkg_color)
            if nargin < 2; cell_color = [0.5 0.5 0.5]; end
            if nargin < 3; bkg_color = [0 0 0]; end
            if isempty(bkg_color); bkg_color = [0 0 0]; end
            if isempty(cell_color); cell_color = [0.5 0.5 0.5]; end
            if isempty(cell_mask); return; end
            
            cell_count = max(cell_mask(:));
            mask_3d = ndims(cell_mask) >= 3;
            X = size(cell_mask,2);
            Y = size(cell_mask,1);
            if mask_3d
                Z = size(cell_mask,3);
            else
                Z = 1;
            end
            
            %Render cell mask image
            render_struct = CellsegDrawer.generateRenderer(X, Y, Z);
            render_struct.r = NaN(Y,X,Z); render_struct.r(:,:,:) = bkg_color(1);
            render_struct.g = NaN(Y,X,Z); render_struct.g(:,:,:) = bkg_color(2);
            render_struct.b = NaN(Y,X,Z); render_struct.b(:,:,:) = bkg_color(3);
            cell_colored = find(cell_mask);
            render_struct.r(cell_colored) = cell_color(1);
            render_struct.g(cell_colored) = cell_color(2);
            render_struct.b(cell_colored) = cell_color(3);
            render_struct.r_draw = uint8(render_struct.r * 255.0); render_struct.r = [];
            render_struct.g_draw = uint8(render_struct.g * 255.0); render_struct.g = [];
            render_struct.b_draw = uint8(render_struct.b * 255.0); render_struct.b = [];
            
            %Determine cell positions for labels
            render_struct.cell_pos = zeros(cell_count, 4);
            if mask_3d
                for c = 1:cell_count
                    cell_filter = logical(cell_mask == c);
                    if nnz(cell_filter) == 0
                        continue;
                    end
                    
                    rp = regionprops3(cell_filter,'BoundingBox','Area');
                    cell_box = rp.BoundingBox;
                    render_struct.cell_pos(c,1) = cell_box(1) + round(cell_box(4)/2);
                    render_struct.cell_pos(c,2) = cell_box(2) + round(cell_box(5)/2);
                    render_struct.cell_pos(c,3) = cell_box(3);
                    render_struct.cell_pos(c,4) = cell_box(3) + cell_box(6);
                end
            else
                for c = 1:cell_count
                    cell_filter = logical(cell_mask == c);
                    if nnz(cell_filter) == 0
                        continue;
                    end
                    
                    rp = regionprops(cell_filter,'BoundingBox','Area');
                    cell_box = rp.BoundingBox;
                    render_struct.cell_pos(c,1) = cell_box(1) + round(cell_box(3)/2);
                    render_struct.cell_pos(c,2) = cell_box(2) + round(cell_box(4)/2);
                    render_struct.cell_pos(c,3) = 1;
                    render_struct.cell_pos(c,4) = 1;
                end
            end
        end
        
        %% ========================== Drawing ==========================
        
        %%
        function drawCellsToFigure(render_struct, z, fighandle, label_color)
            if nargin < 4; label_color = []; end %Don't label.
            if isempty(render_struct); return; end
            if z < 1; z = 1; end
            if z > render_struct.dimZ; z = render_struct.dimZ; end
            
            figure(fighandle);
            hold on;
            
            draw_image = uint8(zeros(render_struct.dimY, render_struct.dimX, 3));
            draw_image(:,:,1) = render_struct.r_draw(:,:,z);
            draw_image(:,:,2) = render_struct.g_draw(:,:,z);
            draw_image(:,:,3) = render_struct.b_draw(:,:,z);
            
            imshow(draw_image);
            
            if ~isempty(label_color)
                cell_count = size(render_struct.cell_pos,1);
                for c = 1:cell_count
                    if z < render_struct.cell_pos(c,3); continue; end
                    if z > render_struct.cell_pos(c,4); continue; end
                    tx = render_struct.cell_pos(c,1);
                    ty = render_struct.cell_pos(c,2);
                    text(tx, ty, num2str(c), 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'Color', label_color, 'FontSize', 18, ...
                        'FontWeight', 'bold');
                end
            end
        end
        
        %% ========================== Interface ==========================
        
        
        
    end
end