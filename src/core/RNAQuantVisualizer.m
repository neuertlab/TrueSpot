%%
classdef RNAQuantVisualizer
    
    %   space (or any key) - Ready listener
	%   x - Exit
	%   [ - Down 1 z-slice
    %   ] - Up 1 z-slice
    
    %%
    properties
        img_raw; %Original image channel
        render_struct; %Renderer
        
        cell_rna_data; %Quant output
        
        fig_handle;
        crosshair_color = [0 0 0];
        cell_color = [0.5 0.5 0.5];
        nuc_color = [0 0 1];
        sig1_color = [1 0 0];
        sig2_color = [0 1 0];
        
        imgraw_lmin = 0;
        imgraw_lmax = 65535;
        
        current_slice; %Current z-slice in view (either for orig image or 3D view) 
        gui_open = false; %true if GUI is open. Affects if things are redrawn
        loop_breaker; %For async control
    end
    
    %%
    methods
        
        %%
        function obj = initialize(obj, raw_image_channel, quant_cells)
            %Prepares everything and renders the image.
            %Does not draw figure, though.
            if ~isempty(obj.fig_handle)
                close(obj.fig_handle);
            end
            
            obj.gui_open = false;
            obj.img_raw = uint16(raw_image_channel); %Scaled to max proj contrast when drawn.
            obj.cell_rna_data = quant_cells;
            obj.loop_breaker = 0;
            
            i3d = ndims(obj.img_raw) >= 3;
            
            X = size(obj.img_raw,2);
            Y = size(obj.img_raw,1);
            Z = 1;
            obj.current_slice = 1;
            
            if i3d
                Z = size(obj.img_raw,3);
                obj.current_slice = max(round(Z/2), 1);
            end
            obj.render_struct = RNAQuantVisualizer.generateRenderer(X,Y,Z);
            
            obj.render_struct.cell_color = obj.cell_color;
            obj.render_struct.nuc_color = obj.nuc_color;
            obj.render_struct.rna1_color = obj.sig1_color;
            obj.render_struct.rna2_color = obj.sig2_color;
            
            %Determine spot/cloud render ranges...
            ch_min = min(obj.img_raw(:));
            ch_max = max(obj.img_raw(:));
            obj.render_struct.range_alpha = [ch_min ch_max];
            
            %Determine original image contrast
            max_proj = double(max(obj.img_raw,[],3));
            obj.imgraw_lmin = min(max_proj(:));
            obj.imgraw_lmax = median(max_proj(:)) + round(10 * std(max_proj(:)));
            
            %Render.
            obj.render_struct = RNAQuantVisualizer.renderQuantResults(obj.render_struct, quant_cells);
        end
        
        %%
        function obj = launchGUI(obj, figno)
            %Draws figure using existing rendered image.
            if obj.gui_open; return; end
            
            obj.fig_handle = figure(figno);
            obj = obj.drawImages();
            
            obj.gui_open = true;
            obj.loop_breaker = 0;
            while obj.loop_breaker == 0
                w = waitforbuttonpress;
                if w
                    obj = obj.onReadyKey();
                end
            end
        end
        
        %%
        function obj = drawImages(obj)
            figure(obj.fig_handle);
            
            %Raw image
            subplot(1,2,1);
            imshow(obj.img_raw(:,:,obj.current_slice), [obj.imgraw_lmin obj.imgraw_lmax]);
            title('Original Image');
            
            %Quant render
            subplot(1,2,2);
            RNAQuantVisualizer.drawToCurrentPlot(obj.render_struct, obj.current_slice);
            title('Quantifier Fit');
            
            if obj.render_struct.dimZ > 1
                sgtitle(['RNA quantification results (z = ' num2str(obj.current_slice) ')']);
            else
                sgtitle('RNA quantification results');
            end
        end
        
        %%
        function obj = onReadyKey(obj)
            if ~obj.gui_open; return; end
            [~,~,btn] = ginput_color(1, obj.crosshair_color);
            
            if btn == 'x' %Exit
                obj = obj.exitGUI();
            elseif btn == '[' %Z down
                newZ = obj.current_slice - 1;
                if newZ < 1
                    fprintf('Already at bottom!\n');
                else
                    obj.current_slice = newZ;
                    obj = obj.drawImages();
                end
            elseif btn == ']' %Z up
                newZ = obj.current_slice + 1;
                if newZ > obj.render_struct.dimZ
                    fprintf('Already at top!\n');
                else
                    obj.current_slice = newZ;
                    obj = obj.drawImages();
                end
            end
        end
        
        %%
        function obj = exitGUI(obj)
            obj.loop_breaker = 1;
            if ~isempty(obj.fig_handle)
                close(obj.fig_handle);
            end
            obj.gui_open = false;
        end
        
    end
    
    %%
    methods (Static)
        
        %% ========================== Structs ==========================
        
        %%
        function render_struct = generateRenderer(X, Y, Z)
            render_struct = struct('dimX', X);
            render_struct.dimY = Y;
            render_struct.dimZ = Z;
            render_struct.range_alpha = [0 65536];
            render_struct.r = [];
            render_struct.g = [];
            render_struct.b = [];
            render_struct.r_draw = [];
            render_struct.g_draw = [];
            render_struct.b_draw = [];
            render_struct.cell_color = [0.500 0.500 0.500];
            render_struct.nuc_color =  [0.000 0.000 1.000];
            render_struct.rna1_color = [1.000 0.000 0.000];
            render_struct.rna2_color = [0.000 1.000 0.000];
        end
        
        %% ========================== Rendering Common ==========================
        
        %%
        function img_out = rescaleImageIntensities(img_in, scale_range)
            base_val = scale_range(1);
            stretch_val = scale_range(2) - scale_range(1);
            img_out = (img_in - base_val) ./ stretch_val;
        end
        
        %% ========================== Rendering Cells ==========================
        
        %%
        function render_struct = renderCellCloud(render_struct, abs_coords, my_cloud, color)
            if isempty(render_struct); return; end
            if isempty(my_cloud); return; end
            if isempty(abs_coords); return; end
            
            if isempty(color)
                color = [1 0 0];
            end
            
            %For readability
            x0 = abs_coords.x0;
            x1 = abs_coords.x1;
            y0 = abs_coords.y0;
            y1 = abs_coords.y1;
            z0 = abs_coords.z0;
            z1 = abs_coords.z1;
            
            %Resulting mask is treated as alpha, not direct value
            alpha_cloud = RNAQuantVisualizer.rescaleImageIntensities(my_cloud.cloud_data, render_struct.range_alpha);
            invert_alpha = 1.0 - alpha_cloud;
            
            %Red
            mix_data = render_struct.r(y0:y1,x0:x1,z0:z1);
            render_struct.r(y0:y1,x0:x1,z0:z1) = (mix_data .* invert_alpha) + (alpha_cloud .* color(1));
                
            %Green
            mix_data = render_struct.g(y0:y1,x0:x1,z0:z1);
            render_struct.g(y0:y1,x0:x1,z0:z1) = (mix_data .* invert_alpha) + (alpha_cloud .* color(2));
                
            %Blue
            mix_data = render_struct.b(y0:y1,x0:x1,z0:z1);
            render_struct.b(y0:y1,x0:x1,z0:z1) = (mix_data .* invert_alpha) + (alpha_cloud .* color(3));
            
        end
        
        %%
        function render_struct = renderCellSpot(render_struct, abs_coords, my_spot, color)
            if isempty(render_struct); return; end
            if isempty(my_spot); return; end
            if isempty(abs_coords); return; end
            
            if isempty(color)
                color = [1 0 0];
            end

            %Simulate
            xy_rad = 4;
            sim_spot = my_spot.generateSimSpotFromFit(xy_rad);
            
            %Shave off x,y,z slices that exceed image boundaries
            %Trim X
            xdim = size(sim_spot,2);
            x0 = abs_coords.x - xy_rad;
            if x0 < 1
                amt = 1 - x0;
                x0 = 1;
                sim_spot = sim_spot(:,amt+1:xdim,:);
                xdim = size(sim_spot,2);
            end
            
            x1 = abs_coords.x + xy_rad;
            if x1 > render_struct.dimX
                amt = x1 - render_struct.dimX;
                x1 = render_struct.dimX;
                sim_spot = sim_spot(:,1:xdim-amt,:);
            end
            
            %Trim Y
            ydim = size(sim_spot,1);
            y0 = abs_coords.y - xy_rad;
            if y0 < 1
                amt = 1 - y0;
                y0 = 1;
                sim_spot = sim_spot(amt+1:ydim,:,:);
                ydim = size(sim_spot,1);
            end
            
            y1 = abs_coords.y + xy_rad;
            if y1 > render_struct.dimY
                amt = y1 - render_struct.dimY;
                y1 = render_struct.dimY;
                sim_spot = sim_spot(1:ydim-amt,:,:);
            end
            
            %Trim Z
            zdim = size(sim_spot,3);
            z_rad = round((zdim - 1) ./ 2);
            z0 = abs_coords.z - z_rad;
            if z0 < 1
                amt = 1 - z0;
                z0 = 1;
                sim_spot = sim_spot(:,:,amt+1:zdim);
                zdim = size(sim_spot,3);
            end
            
            z1 = abs_coords.z + z_rad;
            if z1 > render_struct.dimZ
                amt = z1 - render_struct.dimZ;
                z1 = render_struct.dimZ;
                sim_spot = sim_spot(:,:,1:zdim-amt);
            end
            
            %Rescale to alpha
            alpha_cloud = RNAQuantVisualizer.rescaleImageIntensities(sim_spot, render_struct.range_alpha);
            invert_alpha = 1.0 - alpha_cloud;
            
            %Red
            mix_data = render_struct.r(y0:y1,x0:x1,z0:z1);
            render_struct.r(y0:y1,x0:x1,z0:z1) = (mix_data .* invert_alpha) + (alpha_cloud .* color(1));
                
            %Green
            mix_data = render_struct.g(y0:y1,x0:x1,z0:z1);
            render_struct.g(y0:y1,x0:x1,z0:z1) = (mix_data .* invert_alpha) + (alpha_cloud .* color(2));
                
            %Blue
            mix_data = render_struct.b(y0:y1,x0:x1,z0:z1);
            render_struct.b(y0:y1,x0:x1,z0:z1) = (mix_data .* invert_alpha) + (alpha_cloud .* color(3));
            
        end
        
        %%
        function render_struct = renderCellClouds(render_struct, my_cell)
            if isempty(render_struct); return; end
            if isempty(my_cell); return; end
            if isempty(my_cell.clouds); return; end
            
            x0 = my_cell.cell_loc.left;
            x1 = my_cell.cell_loc.right;
            y0 = my_cell.cell_loc.top;
            y1 = my_cell.cell_loc.bottom;
            z0 = 1; z1 = 1;
            if render_struct.dimZ > 1
                z0 = my_cell.cell_loc.z_bottom;
                z1 = my_cell.cell_loc.z_top;
            end
            
            cloud_count = size(my_cell.clouds,2);
            abspos = struct('x0', 1, 'x1', 1, 'y0', 1, 'y1', 1, 'z0', 1, 'z1', 1);
            for c = 1:cloud_count
                my_cloud = my_cell.clouds(c);
                
                %Absolute coords
                abspos.x0 = x0 + my_cloud.cloud_box.left - 1;
                abspos.x1 = x1 + my_cloud.cloud_box.left - 1;
                abspos.y0 = y0 + my_cloud.cloud_box.top - 1;
                abspos.y1 = y1 + my_cloud.cloud_box.top - 1;
                abspos.z0 = z0 + my_cloud.cloud_box.z_bottom - 1;
                abspos.z1 = z1 + my_cloud.cloud_box.z_bottom - 1;
                
                render_struct = RNAQuantVisualizer.renderCellCloud(render_struct, abspos, my_cloud, render_struct.rna1_color);
            end
        end
        
        %%
        function render_struct = renderCellSpots(render_struct, my_cell)
            if isempty(render_struct); return; end
            if isempty(my_cell); return; end
            if isempty(my_cell.spots); return; end
            
            x0 = my_cell.cell_loc.left;
            y0 = my_cell.cell_loc.top;
            z0 = 1;
            if render_struct.dimZ > 1
                z0 = my_cell.cell_loc.z_bottom;
            end
            
            spot_count = size(my_cell.spots,2);
            abspos = struct('x', 1, 'y', 1, 'z', 1);
            for s = 1:spot_count
                my_spot = my_cell.spots(s);
                
                %Absolute coords (check these...)
                abspos.x = my_spot.gauss_fit.xfit + x0 - 1;
                abspos.y = my_spot.gauss_fit.yfit + y0 - 1;
                abspos.z = my_spot.gauss_fit.zabs + z0 - 1;
                
                render_struct = RNAQuantVisualizer.renderCellCloud(render_struct, abspos, my_cloud, render_struct.rna1_color);
            end
        end
        
        %%
        function render_struct = renderCellNucleus(render_struct, my_cell)
            if isempty(render_struct); return; end
            if isempty(my_cell); return; end
            
            nucmask_3d = ndims(my_cell.mask_nuc) >= 3;
            x0 = my_cell.cell_loc.left;
            x1 = my_cell.cell_loc.right;
            y0 = my_cell.cell_loc.top;
            y1 = my_cell.cell_loc.bottom;
            if nucmask_3d
                mask_dbl = double(my_cell.mask_nuc(:,:,:));
                not_mask_dbl = double(~my_cell.mask_nuc(:,:,:));
                z0 = my_cell.cell_loc.z_bottom;
                z1 = my_cell.cell_loc.z_top;
                
                bkg = render_struct.r(y0:y1,x0:x1,z0:z1) .* not_mask_dbl;
                render_struct.r(y0:y1,x0:x1,z0:z1) =...
                    bkg + (mask_dbl(:,:,:) .* render_struct.nuc_color(1));
                
                bkg = render_struct.g(y0:y1,x0:x1,z0:z1) .* not_mask_dbl;
                render_struct.g(y0:y1,x0:x1,z0:z1) =...
                    bkg + (mask_dbl(:,:,:) .* render_struct.nuc_color(2));
                
                bkg = render_struct.b(y0:y1,x0:x1,z0:z1) .* not_mask_dbl;
                render_struct.b(y0:y1,x0:x1,z0:z1) =...
                    bkg + (mask_dbl(:,:,:) .* render_struct.nuc_color(3));
                clear mask_dbl;
                clear not_mask_dbl;
            else
                mask_dbl = double(my_cell.mask_nuc(:,:));
                not_mask_dbl = double(~my_cell.mask_nuc(:,:));
                Z = render_struct.dimZ;
                
                for z = 1:Z
                    bkg = render_struct.r(y0:y1,x0:x1,z) .* not_mask_dbl;
                    render_struct.r(y0:y1,x0:x1,z) =...
                        bkg + (mask_dbl(:,:) .* render_struct.nuc_color(1));
                    
                    bkg = render_struct.g(y0:y1,x0:x1,z) .* not_mask_dbl;
                    render_struct.g(y0:y1,x0:x1,z) =...
                        bkg + (mask_dbl(:,:) .* render_struct.nuc_color(2));
                    
                    bkg = render_struct.b(y0:y1,x0:x1,z) .* not_mask_dbl;
                    render_struct.b(y0:y1,x0:x1,z) =...
                        bkg + (mask_dbl(:,:) .* render_struct.nuc_color(3));
                end
                clear mask_dbl;
                clear not_mask_dbl;
            end
        end
        
        %%
        function render_struct = renderCell(render_struct, my_cell)
            if isempty(render_struct); return; end
            if isempty(my_cell); return; end
            
            %Cell boundaries
            cellmask_3d = ndims(my_cell.mask_cell) >= 3;
            x0 = my_cell.cell_loc.left;
            x1 = my_cell.cell_loc.right;
            y0 = my_cell.cell_loc.top;
            y1 = my_cell.cell_loc.bottom;
            if cellmask_3d
                mask_dbl = double(my_cell.mask_cell(:,:,:));
                not_mask_dbl = double(~my_cell.mask_cell(:,:,:));
                z0 = my_cell.cell_loc.z_bottom;
                z1 = my_cell.cell_loc.z_top;
                bkg = render_struct.r(y0:y1,x0:x1,z0:z1) .* not_mask_dbl;
                render_struct.r(y0:y1,x0:x1,z0:z1) =...
                    bkg + (mask_dbl(:,:,:) .* render_struct.cell_color(1));
                
                bkg = render_struct.g(y0:y1,x0:x1,z0:z1) .* not_mask_dbl;
                render_struct.g(y0:y1,x0:x1,z0:z1) =...
                    bkg + (mask_dbl(:,:,:) .* render_struct.cell_color(2));
                
                bkg = render_struct.b(y0:y1,x0:x1,z0:z1) .* not_mask_dbl;
                render_struct.b(y0:y1,x0:x1,z0:z1) =...
                    bkg + (mask_dbl(:,:,:) .* render_struct.cell_color(3));
                clear mask_dbl;
                clear not_mask_dbl;
            else
                mask_dbl = double(my_cell.mask_cell(:,:));
                not_mask_dbl = double(~my_cell.mask_cell(:,:));
                Z = render_struct.dimZ;
                for z = 1:Z
                    bkg = render_struct.r(y0:y1,x0:x1,z) .* not_mask_dbl;
                    render_struct.r(y0:y1,x0:x1,z) =...
                        bkg + (mask_dbl(:,:) .* render_struct.cell_color(1));
                    
                    bkg = render_struct.g(y0:y1,x0:x1,z) .* not_mask_dbl;
                    render_struct.g(y0:y1,x0:x1,z) =...
                        bkg + (mask_dbl(:,:) .* render_struct.cell_color(2));
                    
                    bkg = render_struct.b(y0:y1,x0:x1,z) .* not_mask_dbl;
                    render_struct.b(y0:y1,x0:x1,z) =...
                        bkg + (mask_dbl(:,:) .* render_struct.cell_color(3));
                end
                clear mask_dbl;
                clear not_mask_dbl;
            end
            
            %Nucleus
            render_struct = RNAQuantVisualizer.renderCellNucleus(render_struct, cell_rna_data(c));
            
            %Clouds
            render_struct = RNAQuantVisualizer.renderCellClouds(render_struct, cell_rna_data(c));
            
            %Spots
            render_struct = RNAQuantVisualizer.renderCellSpots(render_struct, cell_rna_data(c));
        end
        
        %%
        function render_struct = renderQuantResults(render_struct, cell_rna_data)
            %Reset rendered data to 0
            render_struct.r = double(zeroes(render_struct.dimY, render_struct.dimX, render_struct.dimZ));
            render_struct.g = double(zeroes(render_struct.dimY, render_struct.dimX, render_struct.dimZ));
            render_struct.b = double(zeroes(render_struct.dimY, render_struct.dimX, render_struct.dimZ));
            
            if isempty(cell_rna_data); return; end
            cell_count = size(cell_rna_data,2);
            for c = 1:cell_count
                render_struct = RNAQuantVisualizer.renderCell(render_struct, cell_rna_data(c));
            end
            
            %Scale down to 8-bit channels to save memory.
            render_struct.r_draw = uint8(render_struct.r .* 255);
            render_struct.g_draw = uint8(render_struct.g .* 255);
            render_struct.b_draw = uint8(render_struct.b .* 255);
            render_struct.r = [];
            render_struct.g = [];
            render_struct.b = [];
        end
        
        %% ========================== Drawing ==========================
        
        %%
        function drawToCurrentPlot(render_struct, z)
            rgb_slice = uint8(zeros(render_struct.dimY, render_struct.dimX, 3));
            rgb_slice(:,:,1) = render_struct.r_draw(:,:,z);
            rgb_slice(:,:,2) = render_struct.g_draw(:,:,z);
            rgb_slice(:,:,3) = render_struct.b_draw(:,:,z);
            imshow(rgb_slice);
        end
        
        %%
        function drawToTIF(render_struct, tifpath)
            tifop.color = true;
            tifop.overwrite = true;
            
            Z = render_struct.dimZ;
            tif_dat = uint8(zeros(render_struct.dimY, render_struct.dimX, 3, Z));
            for z = 1:Z
                tif_dat(:,:,1,z) = render_struct.r_draw(:,:,z);
                tif_dat(:,:,2,z) = render_struct.g_draw(:,:,z);
                tif_dat(:,:,3,z) = render_struct.b_draw(:,:,z);
            end
            saveastiff(tif_dat, tifpath, tifop);
        end
        
    end
    
end