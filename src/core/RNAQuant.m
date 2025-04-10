%
%Methods for RNA quantification
%Blythe Hospelhorn
%Modified from code written by Ben Kesler & Gregor Neuert
%Version 1.1.1
%Updated Dec 12, 2022

%Update Log:
%   1.0.0 | 22.10.14
%       Init Doc
%   1.1.0 | 22.12.10
%       Round 1 debugging. Appears to work thru Gaussian fit.
%   1.1.1 | 22.12.12
%       Round 2 debugging. Fixed parallel, tweaked some cloud stuff.

classdef RNAQuant
    
    %%
    properties (Constant)
        GPARAM_MU1 = 1; %X gaussian radius
        GPARAM_MU2 = 2; %Y gaussian radius
        GPARAM_S1 = 3; %XY (or X) gaussian strength
            
        GPARAM_S2 = 4; %Y gaussian strength
        GPARAM_A_1DIM = 4; %Gaussian peak intensity
        GPARAM_A_2DIM = 5; %Gaussian peak intensity
            
        DEFOC_PARAM_WO = 1;
        DEFOC_PARAM_C = 2;
        DEFOC_PARAM_D = 3;
        DEFOC_PARAM_A = 4;
        DEFOC_PARAM_B = 5;
        
        MIN_CLOUD_VOL_DEFO = 7*7*3;
    end
    
    %%
    methods (Static)
        
        %% ========================== Structs ==========================
        
        %%
        function quant_info_struct = genRNAQuantInfoStruct()
            %Params
            quant_info_struct = struct('img_raw', []);
            quant_info_struct.threshold = 0; %This is now used for final counts only
            quant_info_struct.thresholds = [];
            quant_info_struct.t_coord_table = [];
            quant_info_struct.small_obj_size = 3;
            quant_info_struct.gaussian_radius = 7;
            quant_info_struct.connect_size = 8;
            quant_info_struct.spotzoom_r_xy = 4;
            quant_info_struct.spotzoom_r_z = 2;
            quant_info_struct.workdir = [];
            quant_info_struct.cell_mask = [];
            quant_info_struct.nuc_mask = [];
            quant_info_struct.do_clouds = false;
            quant_info_struct.do_gauss_fit = true;
            quant_info_struct.z_adj = 1.0;
            quant_info_struct.workers = 1;
            quant_info_struct.dbgcell = 0;
            quant_info_struct.incl_cell_zero = false;
            
            %Results
            quant_info_struct.plane_stats = [];
            quant_info_struct.cell_rna_data = [];
            quant_info_struct.cell_zero = [];
            quant_info_struct.clouds_mask = [];
        end
        
        %%
        function param_struct = genEmptyGaussFitParamStruct(init_me)
            if init_me
                param_struct.mu1 = 0.0;
                param_struct.mu2 = 0.0;
                param_struct.s1 = 0.0;
                param_struct.s2 = 0.0;
                param_struct.A = 0.0;
                param_struct.xFWHM = 0.0;
                param_struct.yFWHM = 0.0;
            else
                param_struct = struct('mu1', {}, 'mu2', {}, 's1', {}, 's2', {}, ...
                    'A', {}, 'xFWHM', {}, 'yFWHM', {});
            end
        end
        
        %%
        function param_struct = genGaussFitParamStruct(mu1, mu2, s1, s2, A)
            param_struct.mu1 = mu1;
            param_struct.mu2 = mu2;
            param_struct.s1 = s1;
            param_struct.s2 = s2;
            param_struct.A = A;
            param_struct.xFWHM = 2*sqrt(2*log(2)).*s1;
            param_struct.yFWHM = 2*sqrt(2*log(2)).*s2;
        end
        
        %%
        function param_vec = genGaussFitParamVector(mu1, mu2, s1, A, s2)
            if nargin < 5
                param_vec = [mu1 mu2 s1 A];
            else
                param_vec = [mu1 mu2 s1 s2 A];
            end
        end
        
        %%
        function param_vec = genGauss3FitParamVector(mu1, mu2, mu3, s1, s2, s3, A)
            param_vec = [mu1 mu2 mu3 s1 s2 s3 A];
        end
        
        %%
        function param_struct = genGaussFitResultStruct()
            param_struct.xfit = NaN;        %x postion from gaussian fit.
            param_struct.yfit = NaN;        %y postion from gaussian fit.
            param_struct.xinit = NaN;       %x initial position from nongaussian fit.
            param_struct.yinit = NaN;       %y initial position from nongaussian fit.
            param_struct.xgw = NaN;         %x gaussian width
            param_struct.ygw = NaN;         %y gaussian width
            param_struct.xFWHM = NaN;       %x Full Width at Half Max of gaussian fit
            param_struct.yFWHM = NaN;       %y Full Width at Half Max of gaussian fit
            param_struct.expMInt = NaN;     %Maximum RNA pixel intensity experiment
            param_struct.fitMInt = NaN;     %Maximum RNA spot intensity from fit
            param_struct.TotExpInt = NaN;   %Integrated Pixel Intensity for experiment (with background subtracted out)
            param_struct.TotFitInt = NaN;   %Integrated gaussian intensity fit (with background subtraction)
            param_struct.r = NaN;           %Ratio of long to short axis of experimental spot
            param_struct.rFit = NaN;        %Ratio of long to short axis of Fitted gaussian spot
            param_struct.back = NaN;        %Background intensity of the RNA spot
            param_struct.xsem = NaN;        %x fit Standard Error of the Mean in pixels found from the number of photons at each pixel
            param_struct.ysem = NaN;        %y fit Standard Error of the Mean in pixels found from the number of photons at each pixel
            param_struct.zxfit = NaN;       %z position from x fit
            param_struct.zyfit = NaN;       %z position from y fit
            param_struct.zinit = NaN;       %z initial postion from the nongaussian fit.  Integer number
            param_struct.zint = NaN;        %z position from intensity variance through stacks
            param_struct.zrel = NaN;        %Relative z position from the mean of the three z fits
            param_struct.zabs = NaN;        %absolute position in the z stack from the fits
            param_struct.zstd = NaN;        %Nan standard deviation of the z position from the three estimates
            param_struct.zqFit = NaN;       %Quality of zfit.  If <0 pass if >0 curvature is positive at some point and fit is fail
            param_struct.nucRNA = false;    %Is in nucleus?
            param_struct.distRNA = NaN;     %distance in pixels of the nuclear RNA to the nuclear membrane
            param_struct.szNUC = NaN;       %size in square pixels of the nucleus in the z plane of the RNA spot
            param_struct.cytoRNA = false;   %Is in cytoplasm?
            param_struct.xabsloc = NaN;     %x Absolute location in the image.  Column!
            param_struct.yabsloc = NaN;     %y aboslute location in image. Row! Plot x y on RNA image is in FunPlotGen.m function.  Must plot(yabsloc,xabsloc) on image!
            param_struct.normdist = NaN;    %Distance of the nuclear RNA to the nuclear membrane normalized by the radius of the nucleus
            param_struct.nucR = NaN;        %Ratio of long axis to short axis in nucleus.
        end
        
        %% ========================== Fitter Functions ==========================
        
         %%
        function g = Gauss3DFit(param_struct, xydim)
            mu1 = param_struct.mu1;
            mu2 = param_struct.mu2;
            s1 = param_struct.s1;
            A = param_struct.A;
            [XX,YY] = meshgrid(1:1:xydim);
            
            x_factor = (XX - mu1 - 1).^2;
            y_factor = (YY - mu2 - 1).^2;
            w_factor = 2 * (s1^2);
            
            g = A * exp(-((x_factor + y_factor) ./ w_factor));
        end
        
        %%
        function g = Gauss3DFit_Vec(param_vec, xydim)
            mu1 = param_vec(1); %X radius
            mu2 = param_vec(2); %Y radius
            s1 = param_vec(3); %XY width
            A = param_vec(4); %Peak
            [XX,YY] = meshgrid(1:1:xydim);
            
            x_factor = (XX - mu1 - 1).^2;
            y_factor = (YY - mu2 - 1).^2;
            w_factor = 2 * (s1^2);
            
            g = A * exp(-((x_factor + y_factor) ./ w_factor));
        end
        
        %%
        function g = Gauss3DFitS1S2(param_struct, xydim)
            mu1 = param_struct.mu1;
            mu2 = param_struct.mu2;
            s1 = param_struct.s1;
            s2 = param_struct.s2;
            A = param_struct.A;
            [XX,YY] = meshgrid(1:1:xydim);
            
            x_factor = (XX - mu1 - 1).^2;
            y_factor = (YY - mu2 - 1).^2;
            xw_factor = 2 * (s1^2);
            yw_factor = 2 * (s2^2);
            
            g = A * exp(-((x_factor ./ xw_factor) + (y_factor ./ yw_factor)));
        end
        
        %%
        function g = Gauss3DFitS1S2_Vec(param_vec, xydim)
            mu1 = param_vec(1);
            mu2 = param_vec(2);
            s1 = param_vec(3);
            s2 = param_vec(4);
            A = param_vec(5);
            [XX,YY] = meshgrid(1:1:xydim);
            
            x_factor = (XX - mu1 - 1).^2;
            y_factor = (YY - mu2 - 1).^2;
            xw_factor = 2 * (s1^2);
            yw_factor = 2 * (s2^2);
            
            g = A * exp(-((x_factor ./ xw_factor) + (y_factor ./ yw_factor)));
        end
         
        %%
        function df = Defoc4Order_Vec(param_vec, z)
            wo = param_vec(1);
            c = param_vec(2);
            d = param_vec(3);
            A = param_vec(4);
            B = param_vec(5);
            df = wo.*sqrt(1 + ((z-c)./d).^2 + A.*((z-c)./d).^3 + B.*((z-c)./d).^4 );
        end
        
        %% ========================== Helper Functions ==========================
        
        %%
        function spot_params = checkNucRNA(spot_params, this_cell, cx, cy, cz)
            c_width = this_cell.cell_loc.width_i;
            c_height = this_cell.cell_loc.height_i;
            
            spot_params.nucRNA = this_cell.mask_nuc(cy, cx, cz);
            if spot_params.nucRNA
                %Look for nuc edge in x
                for i = 1:c_width
                    x_right = round(cx + i);
                    x_left = round(cx - i);
                            
                    if (x_right > c_width) | (x_left < 1)
                        %If it's outside the cell box, assumed
                        %outside the nucleus?
                        xloc = i;
                        break;
                    end
                            
                    if ~this_cell.mask_nuc(cy, x_right, cz)
                        xloc = i;       %Number of pixels to nuclear membrane
                        break;
                    elseif ~this_cell.mask_nuc(cy, x_left, cz)
                        xloc = i;
                        break;
                    end
                end
                
                %Look for nuc edge in y
                for i = 1:c_height
                    y_down = round(cy + i);
                    y_up = round(cy - i);
                            
                    if (y_down > c_height) | (y_up < 1)
                        yloc = i;
                        break;
                    end
                            
                    if ~this_cell.mask_nuc(y_down, cx, cz)
                        yloc = i;
                        break;
                    elseif ~this_cell.mask_nuc(y_up, cx, cz)
                        yloc = i;
                        break;
                    end
                end
                
                %Take the smaller distance.
                if xloc < yloc
                    spot_params.distRNA = xloc;  %the result of the distance in voxels of all nuclear RNA
                elseif yloc < xloc
                    spot_params.distRNA = yloc;
                end
                spot_params.szNUC = nnz(this_cell.mask_nuc(:,:,cz)).^2;
                spot_params.normdist = spot_params.distRNA ./ sqrt(spot_params.szNUC ./ pi()); %normalized to radius
            end
        end
        
        %%
        function nuc_ellipticity = findNucEllipticityAtZ(this_cell, cz)
            %walk in four directions from center of image until the
            %edge is reached to find the approximate ellipticity
                    
            c_width = this_cell.cell_loc.width_i;
            c_height = this_cell.cell_loc.height_i;
            
            x_mid = round(c_width./2);
            y_mid = round(c_height./2);
            
            for i = 1:c_width
                xx = x_mid+i;
                if (xx > c_width) | (~this_cell.mask_nuc(y_mid, xx, cz))
                    mxloc = abs(i);
                    break;
                end
            end
                        
            for i = 1:c_width
                xx = x_mid-i;
                if (xx < 1) | (~this_cell.mask_nuc(y_mid, xx, cz))
                    pxloc = abs(i);
                    break;
                end
            end
                        
            
            for i = 1:c_height
                yy = y_mid+i;
                if (yy > c_height) | ~this_cell.mask_nuc(yy, x_mid, cz)
                    myloc = abs(i);
                    break;
                end
            end
                        
            for i = 1:c_height
                yy = y_mid-i;
                if (yy < 1) | ~this_cell.mask_nuc(yy, x_mid, cz)
                    pyloc = abs(i);
                    break;
                end
            end
            
            nuc_ellipticity = (mxloc+pxloc)./(myloc+pyloc);
        end
        
        %%
        function min_cloud_vol = suggestMinCloudVolume(cell_rna_data)
            % > 1 stdev above mean spot volume
            
            %Count spots for vector allocation
            total_spots = 0;
            cell_count = size(cell_rna_data,2);
            for c = 1:cell_count
                this_cell = cell_rna_data(c); %READ ONLY
                spot_count = this_cell.getSpotCount();
                total_spots = total_spots + spot_count;
            end
            
            %Nab spot volumes
            spot_vols = NaN(1,total_spots);
            pos = 1;
            for c = 1:cell_count
                this_cell = cell_rna_data(c); %READ ONLY
                spot_count = this_cell.getSpotCount();
                if spot_count <= 0; continue; end
                for s = 1:spot_count
                    spot_vols(1,pos) = this_cell.spots(s).fit_volume;
                    pos = pos + 1;
                end
            end
            
            %Do the stats.
            vol_mean = nanmean(spot_vols);
            vol_std = nanstd(spot_vols);
            
            min_cloud_vol = vol_mean + vol_std;
        end
        
        %%
        function b_filters = generateGaussianBFilters(xydim)
            b_filters.A = double(ones(xydim,xydim));
            b_filters.A(2:xydim-1,2:xydim-1) = 0;
            b_filters.B = (1.0 - b_filters.A) * 0.5;
            b_filters.C = b_filters.B;
            b_filters.C(3:xydim-2,3:xydim-2) = 2.0;
            b_filters.C = b_filters.C ./ 2.0;
        end
        
        %%
        function cloud_boxes = mergeOverlappingClouds(cloud_boxes, X, Y, Z)
            box_mask = false(Y,X,Z);
            box_count = size(cloud_boxes,1);
            
            for i = 1:box_count
                x0 = cloud_boxes(i,1);
                y0 = cloud_boxes(i,2);
                z0 = cloud_boxes(i,3);
                x1 = cloud_boxes(i,4) + x0;
                y1 = cloud_boxes(i,5) + y0;
                z1 = cloud_boxes(i,6) + z0;
                box_mask(y0:y1,x0:x1,z0:z1) = true;
            end
            
            imrprop = regionprops3(box_mask, 'BoundingBox');
            cloud_boxes = imrprop.BoundingBox;
            
            for i = 1:3
                d0 = cloud_boxes(:,i);
                dd = cloud_boxes(:,i+3);
                d1 = d0 + dd;
                d0 = floor(d0);
                d1 = ceil(d1);
                dd = d1 - d0;
                cloud_boxes(:,i) = d0;
                cloud_boxes(:,i+3) = dd;

                clear d0 d1 dd
            end

        end
        
        %%
        function cloud_mask = fillClouds(cloud_mask)
            %TODO may play with connectivity to see if that works better...
            cloud_mask = imfill(cloud_mask, 'holes');
        end
        
        %% ========================== Fit Primary Functions ==========================

        %%
        function [img_processed, imrmax_mask, plane_stats] = FilterCleanedImage(img_dpcleaned, threshold, small_obj_size, connect_size, gaussian_radius)
            if nargin < 5; gaussian_radius = 7; end
            if nargin < 3; small_obj_size = 3; end
            if nargin < 4; connect_size = 8; end
            
            X = size(img_dpcleaned,2);
            Y = size(img_dpcleaned,1);
            Z = size(img_dpcleaned,3);
            plane_stats(Z) = struct('stdev', 0.0, 'mean', 0.0, 'variance', 0.0, 'median', 0.0);

            img_processed = uint16(zeros(Y,X,Z));
            for z = 1:Z
                z_raw = img_dpcleaned(:,:,z);
                z_filtered = RNA_Threshold_Common.applyGaussianFilter(z_raw, gaussian_radius, 2);
                z_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(z_filtered);
                
                z_thresh = z_filtered > threshold;
                z_thresh = bwareaopen(z_thresh, small_obj_size, connect_size);
                img_processed(:,:,z) = immultiply(z_thresh, z_raw);
                
                plane_stats(z).stdev = std2(z_filtered(:));
                plane_stats(z).mean = mean2(z_filtered(:));
                plane_stats(z).variance = var(double(z_filtered(:)));
                plane_stats(z).median = median(double(z_filtered(:)));
            end
            
            img_processed = RNA_Threshold_Common.blackoutBorders(img_processed, gaussian_radius, 0);
            imrmax_mask = imregionalmax(img_processed,26);
        end

        %%
        function [img_dpcleaned, img_bkg] = RNAProcess3DimFast(img_raw, workdir)
            if nargin < 2; workdir = pwd; end
            
            img_dpcleaned = [];
            img_bkg = [];
        
            if isempty(img_raw) 
                return; 
            end

            dead_pix_path = [workdir filesep 'deadpix.mat'];

            %Pre-filtering...
            img_dpcleaned = uint16(img_raw);
            clear img_raw;
            RNA_Threshold_Common.saveDeadPixels(img_dpcleaned, dead_pix_path);
            img_dpcleaned = RNA_Threshold_Common.cleanDeadPixels(img_dpcleaned, dead_pix_path);
            
            if isfile(dead_pix_path)
                delete(dead_pix_path);
            end

            X = size(img_dpcleaned,2);
            Y = size(img_dpcleaned,1);
            Z = size(img_dpcleaned,3);

            img_bkg = uint16(zeros(Y,X,Z));
            for z = 1:Z
                z_raw = img_dpcleaned(:,:,z);
                img_bkg(:,:,z) = medfilt2(z_raw, [20 20]); 
            end
        end
        
        %%
        function [img_dpcleaned, img_processed, imrmax_mask, img_bkg, plane_stats] = RNAProcess3Dim(img_raw, threshold, small_obj_size, connect_size, gaussian_radius, workdir)
            %B2_RNAprocess3Dim.m
            [img_dpcleaned, img_bkg] = RNAQuant.RNAProcess3DimFast(img_raw, workdir);
            [img_processed, imrmax_mask, plane_stats] = RNAQuant.FilterCleanedImage(img_dpcleaned, threshold, small_obj_size, connect_size, gaussian_radius);
        end
        
        %%
        function [cell_rna_data, cell_zero] = RNAPosFromTable(img_raw, coords, cell_mask, nuc_mask, include_cell_zero)
            if nargin < 4; include_cell_zero = false; end

            Z = 1;
            if ndims(img_raw) >= 3
                Z = size(img_raw, 3);
            end
            
            cell_count = max(cell_mask(:));

            cell_rna_data(cell_count) = SingleCell;
            for i = cell_count:-1:1
                cell_rna_data(i) = SingleCell.newRNACell(i, 0);
                cell_rna_data(i).dim_z = Z;
            end

            if include_cell_zero
                cell_zero = SingleCell.newRNACell(0, 0);
                cell_zero.dim_z = Z;
            else
                cell_zero = [];
            end

            cstart = 1;
            if include_cell_zero; cstart = 0; end
            for c = cstart:cell_count
                if c > 0
                    myCell = cell_rna_data(c);
                else
                    myCell = cell_zero;
                end
                myCell = myCell.findBoundaries(cell_mask, nuc_mask);

                [rna_in_cell, rna_in_nuc] = myCell.getCoordsSubset(coords);
                rna_count_total = 0;
                rna_count_nuc = 0;
                if ~isempty(rna_in_cell)
                    rna_count_total = size(rna_in_cell,1);
                end

                if ~isempty(rna_in_nuc)
                    rna_count_nuc = size(rna_in_nuc,1);
                end

                myCell.spotcount_total = rna_count_total;
                myCell.spotcount_nuc = rna_count_nuc;
                myCell.spotcount_cyto = rna_count_total - rna_count_nuc;

                if rna_count_total > 0
                    myCell = myCell.preallocSpots(rna_count_total);
                    for i = 1:rna_count_total
                        myCell.spots(i).x = rna_in_cell(i,1);
                        myCell.spots(i).y = rna_in_cell(i,2);
                        myCell.spots(i).z = rna_in_cell(i,3);
                        if size(rna_in_cell, 2) > 3
                            myCell.spots(i).dropout_thresh = rna_in_cell(i,4);
                        end
                    end
                end

                if c > 0
                    cell_rna_data(c) = myCell;
                else
                    cell_zero = myCell;
                end
            end

        end

        %%
        function [cell_rna_data, cell_zero] = NonGaussRNAPos(img_f, imrmax_mask, cell_mask, nuc_mask, include_cell_zero)
            if nargin < 5; include_cell_zero = false; end
            
            %B3_nongaussRNApos
            Z = 1;
            if ndims(img_f) >= 3
                Z = size(img_f, 3);
            end
            
            cell_count = max(cell_mask(:));

            cell_rna_data(cell_count) = SingleCell;
            for i = cell_count:-1:1
                cell_rna_data(i) = SingleCell.newRNACell(i, 0);
                cell_rna_data(i).dim_z = Z;
            end
            if include_cell_zero
                cell_zero = SingleCell.newRNACell(0, 0);
                cell_zero.dim_z = Z;
            else
                cell_zero = [];
            end
            
            cstart = 1;
            if include_cell_zero; cstart = 0; end
            img_imrmax = immultiply(img_f, imrmax_mask);
            for c = cstart:cell_count
                if c > 0
                    myCell = cell_rna_data(c);
                else
                    myCell = cell_zero;
                end
                myCell = myCell.findBoundaries(cell_mask, nuc_mask);
                
                X0 = myCell.cell_loc.left;
                X1 = myCell.cell_loc.right;
                Y0 = myCell.cell_loc.top;
                Y1 = myCell.cell_loc.bottom;
                
                cell_mask_dbl = double(myCell.get3DCellMask());
                cell_rna_imrmax = immultiply(double(img_imrmax(Y0:Y1,X0:X1,:)), cell_mask_dbl);
                cell_nuc = immultiply(cell_rna_imrmax, double(myCell.mask_nuc));
                cell_cyto = immultiply(cell_rna_imrmax, double(myCell.mask_cyto));
                
                [lbl_cell, numr_cell] = bwlabeln(cell_rna_imrmax > 0,26); 
                [~, numr_cyto] = bwlabeln(cell_cyto > 0,26); 
                [~, numr_nuc] = bwlabeln(cell_nuc > 0,26); 
                
                if numr_cell > 0
                    myCell = myCell.preallocSpots(numr_cell);
                    for i = 1:numr_cell
                        [row, col, vec] = ind2sub(size(lbl_cell), find(lbl_cell == i));
                        if size(row,1) == 1
                            y = row;
                            x = col;
                            z = vec;
                        elseif size(row,1) > 1
                            len = (size(row,1));
                            tempvec(len) = 0;
                            for j = 1:len
                                tempvec(j) = cell_rna_imrmax(row(j),col(j),vec(j));
                            end
                            [~, indx] = max(tempvec);
                            y = row(indx);
                            x = col(indx);
                            z = vec(indx);
                            clear tempvec;
                        end
                        
                        myCell.spots(i).x = x;
                        myCell.spots(i).y = y;
                        myCell.spots(i).z = z;
                        myCell.spots(i).dropout_thresh = img_f((y+Y0),(x+X0),z);
                    end
                end
                myCell.spotcount_total = numr_cell;
                myCell.spotcount_nuc = numr_nuc;
                myCell.spotcount_cyto = numr_cyto;

                if c > 0
                    cell_rna_data(c) = myCell;
                else
                    cell_zero = myCell;
                end
            end
        end
        
        %%
        function my_cell = FitGaussians3Cell(cell_img_nobkg, my_cell, xy_rad, z_rad, b_filters)

            my_cell = my_cell.calculateNuclearEllipticity(); %Precalculate to reduce redundancy
            my_cell = my_cell.getBasicStats(cell_img_nobkg);
            spot_count = my_cell.getSpotCount();
            if spot_count <= 0; return; end %No spots in this cell. No need to do anything else.
            
            %Calculate some values...
            xydim = (xy_rad * 2) + 1;
            zdim = (z_rad * 2) + 1;
            z_mid = z_rad + 1; %Index of middle slice in spot box.
            
            %For convenience/code tidyness
            X0 = my_cell.cell_loc.left;
            Y0 = my_cell.cell_loc.top;
            
            GPARAM_MU1 = RNAQuant.GPARAM_MU1;
            GPARAM_MU2 = RNAQuant.GPARAM_MU2;
            GPARAM_S1 = RNAQuant.GPARAM_S1;
            
            GPARAM_S2 = RNAQuant.GPARAM_S2;
            GPARAM_A_1DIM = RNAQuant.GPARAM_A_1DIM;
            GPARAM_A_2DIM = RNAQuant.GPARAM_A_2DIM;
            
            DEFOC_PARAM_WO = RNAQuant.DEFOC_PARAM_WO;
            DEFOC_PARAM_C = RNAQuant.DEFOC_PARAM_C;
            DEFOC_PARAM_D = RNAQuant.DEFOC_PARAM_D;
            DEFOC_PARAM_A = RNAQuant.DEFOC_PARAM_A;
            DEFOC_PARAM_B = RNAQuant.DEFOC_PARAM_B;
            
            spot_bkg(zdim) = 0.0;
            axis_ratios(zdim) =  0.0;
            for s = 1:spot_count
                %fprintf("DEBUG || RNAQuant.FitGaussians3Cell -- Cell %d: Fitting spot %d of %d...\n", my_cell.cell_number, s, spot_count);
                my_spot = my_cell.spots(s); %READ ONLY
                        
                %Isolate the pixels immediately around the spot
                spot_data = RNAUtils.isolateSpotData(cell_img_nobkg, my_spot.x, my_spot.y, my_spot.z, xy_rad, z_rad);
                
                spot_filtered = NaN(xydim, xydim, zdim);
                for z = 1:zdim
                    spot_bkg(z) = sum(sum(immultiply(spot_data(:,:,z),double(b_filters.A)),1))./sum(b_filters.A(:));
                    spot_filtered(:,:,z) = immultiply(spot_data(:,:,z),double(b_filters.B)); %Not used? Is this a bug?
                    spot_filtered(:,:,z) = round(immultiply(spot_data(:,:,z),double(b_filters.C))./1);
                end
                
                %Try to fit Gaussian to spot
                gparams_vec = RNAQuant.genGaussFitParamVector(2, 3, 1, double(max(spot_filtered(:))));
                %[mxX, mxY] = meshgrid(1:1:xydim);
                options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9);
                
                gauss_sim = NaN(xydim, xydim, zdim);
                fitted_params = NaN(5,zdim);
                for z = 1:zdim
                    slice_filtered = double(spot_filtered(:,:,z));
                    slice_scaled = slice_filtered ./ max(slice_filtered(:));
                    slice_scaled = im2bw(slice_scaled,0.5);
                    iprops = regionprops(slice_scaled,'MinorAxisLength','MajorAxisLength');
                    try 
                        axis_ratios(z) = iprops.MajorAxisLength ./ iprops.MinorAxisLength; %try/catch?
                    catch
                        axis_ratios(z) = NaN;
                    end
                    
                    lb = [2 2 gparams_vec(GPARAM_S1) gparams_vec(GPARAM_A_1DIM)];
                    ub = [xydim-2 xydim-2  gparams_vec(GPARAM_S1) 5*gparams_vec(GPARAM_A_1DIM)] + 0.001;
                    [resparams,~,~,~,~,~,~] = lsqcurvefit(@RNAQuant.Gauss3DFit_Vec, gparams_vec, xydim, slice_filtered, lb, ub, options);
                    
                    lb = [resparams(GPARAM_MU1) resparams(GPARAM_MU2) 0 0];
                    ub = [resparams(GPARAM_MU1) resparams(GPARAM_MU2) 5 2*resparams(GPARAM_A_1DIM)] + 0.001;
                    [resparams,~,~,~,~,~,~] = lsqcurvefit(@RNAQuant.Gauss3DFit_Vec, gparams_vec, xydim, slice_filtered, lb, ub, options);
                    
                    lb = [resparams(GPARAM_MU1) resparams(GPARAM_MU2) 0 0 0];
                    ub = [resparams(GPARAM_MU1) resparams(GPARAM_MU2) resparams(GPARAM_S1)*3 resparams(GPARAM_S1)*3 resparams(GPARAM_A_1DIM)*3] + 0.001;
                    resparams([GPARAM_MU1 GPARAM_MU2 GPARAM_S1 GPARAM_A_2DIM]) = resparams;
                    resparams(GPARAM_S2) = resparams(GPARAM_S1);
                    [resparams,~,~,~,~,~,~] = lsqcurvefit(@RNAQuant.Gauss3DFitS1S2_Vec, resparams, xydim, slice_filtered, lb, ub, options);
                    
                    gauss_sim(:,:,z) = RNAUtils.generateGaussian2D(...
                        xydim, xydim, resparams(GPARAM_MU1), resparams(GPARAM_MU2), ...
                        resparams(GPARAM_S1), resparams(GPARAM_S2), resparams(GPARAM_A_2DIM));
                    fitted_params(:,z) = resparams;
                    fitted_params(GPARAM_MU1,z) = fitted_params(GPARAM_MU1,z) + my_spot.x - xy_rad;
                    fitted_params(GPARAM_MU2,z) = fitted_params(GPARAM_MU2,z) + my_spot.y - xy_rad;
                end
                
                %Save some data to the return struct
                my_cell.spots(s) = my_cell.spots(s).saveFittedParamsFromVector(fitted_params);
                spot_params = my_cell.spots(s).gauss_fit; %No idea if it copies if saved back. But wth writing out the whole thing is so ugly.
                spot_params.xfit = fitted_params(GPARAM_MU1,z_mid);
                spot_params.yfit = fitted_params(GPARAM_MU2,z_mid);
                spot_params.xinit = my_spot.x;
                spot_params.yinit = my_spot.y;
                spot_params.xgw = fitted_params(GPARAM_S1,z_mid);
                spot_params.ygw = fitted_params(GPARAM_S2,z_mid);
                spot_params.xFWHM = 2*sqrt(2*log(2)).*fitted_params(GPARAM_S1,z_mid);
                spot_params.yFWHM = 2*sqrt(2*log(2)).*fitted_params(GPARAM_S2,z_mid);
                spot_params.expMInt = double(max(spot_filtered(:)));
                spot_params.fitMInt = fitted_params(GPARAM_A_2DIM,z_mid);
                spot_params.TotExpInt = nansum(nansum(spot_filtered(:,:,z_mid)));
                spot_params.TotFitInt = nansum(nansum(gauss_sim(:,:,z_mid)));
                spot_params.r = axis_ratios(z_mid);
                if spot_params.xFWHM > spot_params.yFWHM
                    spot_params.rFit = spot_params.xFWHM./spot_params.yFWHM;
                else
                    spot_params.rFit = spot_params.yFWHM./spot_params.xFWHM;
                end
                spot_params.back = spot_bkg(z_mid);
                
                %NOTE: In Gregor & Ben's code this is scaled to nm
                %using the voxel size.
                %Because I don't want dependence on that here, this
                %value is measured in units of VOXELS.
                spot_params.xsem = sqrt( (spot_params.xgw^2)/spot_params.TotExpInt + ...
                    (1/(12*spot_params.TotExpInt) + (8*pi()*spot_params.xgw^4) * spot_params.back.^2) ...
                    ./(spot_params.TotExpInt.^2));
                spot_params.ysem = sqrt( (spot_params.ygw^2)/spot_params.TotExpInt + ...
                    (1/(12*spot_params.TotExpInt) + (8*pi()*spot_params.ygw^4) * spot_params.back.^2) ...
                    ./(spot_params.TotExpInt.^2));
                
                %Determine a z position for the Gaussian
                %1. Using x width
                [~,minidx] = min(fitted_params(GPARAM_S1,:));
                if (minidx > 1 & minidx < 5)
                    df_params = double([1 3 2 0.67 0.3]);
                    lb = [0 0 0 0 0]; % lower bounds for finding the gaussian width
                    ub = [10 10 10 10 10] + 0.001; % upper bounds for finding the gaussian width
                    z_vec = 1:5;
                    options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9); % fit conditions
                    [df_params] = lsqcurvefit(@RNAQuant.Defoc4Order_Vec, df_params, z_vec, fitted_params(GPARAM_S1,:), lb, ub, options);
                    z_vec = [1:0.01:5];
                    tempvec = df_params(DEFOC_PARAM_WO) .* sqrt(1 + ((z_vec - df_params(DEFOC_PARAM_C)) ./ df_params(DEFOC_PARAM_D)).^2 + ...
                        df_params(DEFOC_PARAM_A) .* ((z_vec - df_params(DEFOC_PARAM_C)) ./ df_params(DEFOC_PARAM_D)).^3 + ...
                        df_params(DEFOC_PARAM_B) .* ((z_vec - df_params(DEFOC_PARAM_C)) ./ df_params(DEFOC_PARAM_D)).^4 );
                    [~,minidx] = min(tempvec);
                    spot_params.zxfit = z_vec(minidx); %gm1
                else
                    spot_params.zxfit = NaN; %gm1
                end
                
                %2. Using y width
                [~,minidx] = min(fitted_params(GPARAM_S2,:));
                if (minidx > 1 & minidx < 5)
                    df_params = double([1 3 2 0.67 0.3]);
                    lb = [0 0 0 0 0]; % lower bounds for finding the gaussian width
                    ub = [10 10 10 10 10] + 0.001; % upper bounds for finding the gaussian width
                    z_vec = 1:5;
                    options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9); % fit conditions
                    [df_params] = lsqcurvefit(@RNAQuant.Defoc4Order_Vec, df_params, z_vec, fitted_params(GPARAM_S2,:), lb, ub, options);
                    z_vec = 1:0.01:5;
                    tempvec = df_params(DEFOC_PARAM_WO) .* sqrt(1 + ((z_vec - df_params(DEFOC_PARAM_C)) ./ df_params(DEFOC_PARAM_D)).^2 + ...
                        df_params(DEFOC_PARAM_A) .* ((z_vec - df_params(DEFOC_PARAM_C)) ./ df_params(DEFOC_PARAM_D)).^3 + ...
                        df_params(DEFOC_PARAM_B) .* ((z_vec - df_params(DEFOC_PARAM_C)) ./ df_params(DEFOC_PARAM_D)).^4 );
                    [~,minidx] = min(tempvec);
                    spot_params.zyfit = z_vec(minidx); %gm2
                else
                    spot_params.zyfit = NaN; %gm2
                end
                
                %3. Using intensity
                [~,maxidx] = max(fitted_params(GPARAM_A_2DIM,:));
                if (maxidx > 1 & maxidx < 5)
                    p3 = polyfit([1:5], fitted_params(GPARAM_A_2DIM,:), 4);
                    tempvec = 1:0.01:5;
                    [~,maxidx] = max(polyval(p3,tempvec));
                    %spot_params.zqFit = max(diff(diff(polyval(p3,tempvec),2)));
                    spot_params.zqFit = max(diff(polyval(p3,tempvec),2)); %I think this was a typo
                    spot_params.zint = tempvec(maxidx);
                else
                    spot_params.zint = NaN; %gm3
                    spot_params.zqFit = NaN; %r1
                end
                
                spot_params.zrel = mean([spot_params.zxfit spot_params.zyfit spot_params.zint], 'all', 'omitnan');
                spot_params.zstd = std([spot_params.zxfit spot_params.zyfit spot_params.zint], 0, 'all', 'omitnan');
                spot_params.zinit = my_spot.z;
                spot_params.zabs = spot_params.zrel + my_spot.z - z_mid;

                if isnan(spot_params.zabs)
                    spot_params.zabs = spot_params.zinit;
                end
                
                %Check if nuclear RNA, and calculate related
                %properties.
                mu1_mid = round(fitted_params(GPARAM_MU1,z_mid));
                mu2_mid = round(fitted_params(GPARAM_MU2,z_mid));
                spot_params = RNAQuant.checkNucRNA(spot_params, my_cell, mu1_mid, mu2_mid, my_spot.z);
                
                spot_params.nucR = my_cell.nuc_ellip(my_spot.z);
                
                spot_params.cytoRNA = my_cell.mask_cyto(mu2_mid, mu1_mid, my_spot.z); % Cytoplasmic
                spot_params.xabsloc = X0 + spot_params.xfit - 1;
                spot_params.yabsloc = Y0 + spot_params.yfit - 1;
                
                my_cell.spots(s).gauss_fit = spot_params;
                my_cell.spots(s) = my_cell.spots(s).sumFromGaussFit();
            end
        end
        
        %%
        function my_cell = RNAGaussThreshCell(my_cell)
            spot_count = my_cell.getSpotCount(); %nn in og code
            if spot_count <= 0; return; end
            
            %smNuc appears to be the nuclear mask inside cell box
            for s = 1:spot_count
                adj_param = my_cell.spots(s).gauss_fit;
                if ~isnan(adj_param.xfit)
                    %See if adjusted, fit position is inside nucleus
                    fit_x = round(adj_param.xfit);
                    fit_y = round(adj_param.yfit);
                    fit_z = round(adj_param.zinit);
                    
                    adj_param = RNAQuant.checkNucRNA(adj_param, my_cell, fit_x, fit_y, fit_z);
                    adj_param.cytoRNA = my_cell.mask_cyto(fit_x, fit_y, fit_z); % Cytoplasmic
                    
                    %I have no idea what you are trying to do here...
                    adj_param.nucR = my_cell.nuc_ellip(adj_param.zinit);
                end
                my_cell.spots(s).adj_gauss_fit = adj_param;
            end
        end
        
        %%
        function [cell_cloud_mask, cloud_boxes] = DetectRNACloudsCell(cell_img_raw, start_coords, z_adj)
            %This just runs the DETECT function, does not mess with Spots.
            %That way, can be run at same time as FitGaussians3Cell
            clouddet = RNA_Clouds;
            [cell_cloud_mask, cloud_boxes] = clouddet.detectClouds(cell_img_raw, start_coords, z_adj);
        end
        
        %%
        function my_cell = FinishRNACloudsCell(my_cell, cell_img_nobkg, cell_cloud_mask, cloud_boxes, min_cloud_vol)
            if isempty(cloud_boxes); return; end %Nothing there...
            
            cw = my_cell.cell_loc.width_i;
            ch = my_cell.cell_loc.height_i;
            Z = my_cell.dim_z;

            cloud_boxes = RNAQuant.mergeOverlappingClouds(cloud_boxes, cw, ch, Z);
            cloud_count = size(cloud_boxes,1);
            temp_mem(cloud_count) = RNACloud;
            temp_pos = 1;

            %cell_cloud_mask = RNAQuant.fillClouds(cell_cloud_mask);
            kept_cloud_mask = false(ch,cw,Z);
            for i = 1:cloud_count
                x0 = max(cloud_boxes(i,1) - 4, 1);
                y0 = max(cloud_boxes(i,2) - 4, 1);
                z0 = max(cloud_boxes(i,3) - 4, 1);
                x1 = min(cloud_boxes(i,4) + x0 + 5, cw);
                y1 = min(cloud_boxes(i,5) + y0 + 5, ch);
                z1 = min(cloud_boxes(i,6) + z0 + 5, Z);
                
                my_cloud = RNACloud.newRNACloud();
                my_cloud.cloud_box = ... %Update this if ever need to add 2D compat.
                    SingleCell.generateRecPrismStruct(x0, x1, y0, y1, z0, z1);
                my_cloud.cloud_mask = cell_cloud_mask(y0:y1,x0:x1,z0:z1);
                
                my_cloud = my_cloud.recalculateVolume();
                if my_cloud.cloud_volume < min_cloud_vol; continue; end
                
                %Determine total intensity
                cloud_img = cell_img_nobkg(y0:y1,x0:x1,z0:z1);
                cloud_img = immultiply(cloud_img, my_cloud.cloud_mask);
                my_cloud.total_intensity = sum(cloud_img, 'all');
                my_cloud.cloud_data = cloud_img;
                ccmask = false(ch,cw,Z);
                ccmask(y0:y1,x0:x1,z0:z1) = my_cloud.cloud_mask;
                kept_cloud_mask = or(kept_cloud_mask, ccmask);
                
                %Determine if in nucleus
                cnuc = immultiply(ccmask, my_cell.mask_nuc);
                nuc_pix = sum(cnuc, 'all');
                if nuc_pix > (my_cloud.cloud_volume / 2)
                    my_cloud.is_nuc = true;
                else
                    my_cloud.is_cyto = true;
                end
                
                %Save to temp array
                temp_mem(temp_pos) = my_cloud;
                temp_pos = temp_pos + 1;
            end
            
            %Move passed clouds from temp to cell obj
            if temp_pos > 1
                my_cell.clouds = temp_mem(1:temp_pos-1);
            end
            clear temp_mem;
            
            %Go through previously fitted spots and see if they fall
            %   inside clouds...
            if temp_pos > 1
                if isempty(my_cell.spots); return; end
                cell_spots = size(my_cell.spots, 2);
                for s = 1:cell_spots
                    my_spot = my_cell.spots(s); %READ ONLY
%                     x = my_spot.gauss_fit.xfit;
%                     y = my_spot.gauss_fit.yfit;
%                     z = max(my_spot.gauss_fit.zabs, 1); %TODO This is sometimes 0. That needs to be fixed in the fitter.
%                     z = min(z, Z); %Also going outside max...
                    % Use seed snap instead of fit
                    y = my_spot.y;
                    x = my_spot.x;
                    z = my_spot.z;
                    if kept_cloud_mask(y,x,z)
                        my_cell.spots(s).in_cloud = true;
                    end
                end
            end
        end
        
        %% ========================== Parallelization ==========================
        
        %%
        function rmres = removeParInternalDir(pardir)
            parinternaldir = [pardir filesep 'mlinternal'];
            rmres = rmdir(parinternaldir, 's');
        end
        
        %%
        function workdir = initParallel(threads, workers, workdir)
            rn = rand();
            rn = round(rn * 2000000000);
            subdir_name = sprintf("quantjob_%08x", rn);
            workdir = strjoin([workdir filesep subdir_name],'');
            mkdir(workdir);
            
            pardir = strjoin([workdir filesep 'mlinternal'], '');
            mkdir(pardir);
            
            local_cluster = parcluster();
            local_cluster.NumThreads = threads;
            local_cluster.NumWorkers = workers;
            local_cluster.JobStorageLocation = pardir;
            parpool(local_cluster);
        end
        
        %%
        function shutdownParallel()
            delete(gcp('nocreate'));
        end
        
        %% ========================== Fit Job Distribution ==========================
        
        %%
        function FitGaussians3CellPar(my_cell, xy_rad, z_rad, b_filters, pardir)
            my_cell = RNAQuant.FitGaussians3Cell(my_cell.img_nobkg, my_cell, xy_rad, z_rad, b_filters);
            savepath = [pardir filesep sprintf('gaussians_%04d.mat', my_cell.cell_number)];
            spots = my_cell.spots;
            nuc_ellip = my_cell.nuc_ellip;
            save(savepath, 'spots', 'nuc_ellip', '-v7.3');
        end
        
        %%
        function DetectRNACloudsCellPar(my_cell, z_adj, pardir)
            [cell_cloud_mask, cloud_boxes] = RNAQuant.DetectRNACloudsCell(my_cell.img_raw, my_cell.coord_list, z_adj);
            savepath = [pardir filesep sprintf('clouds_%04d.mat', my_cell.cell_number)];
            save(savepath, 'cell_cloud_mask', 'cloud_boxes', '-v7.3');
        end
        
        %%
        function quant_struct = FitRNA_P(quant_struct, img_clean, img_bkg)
            
            %This looks stupid, but it's a trick to get MATLAB to
            %more efficiently distribute data to workers.
            cell_list_a = quant_struct.cell_rna_data;
            %cell_list_b = [];
            
            %Pre-divide image data into cells
            cell_count = size(cell_list_a,2);
            img_raw_dbl = double(img_clean);
            if quant_struct.no_bkg_subtract
                img_use_dbl = img_raw_dbl;
            else
                img_bkg_dbl = double(img_bkg);
                img_use_dbl = img_raw_dbl - img_bkg_dbl;
            end
            clear img_bkg_dbl;
            for c = 1:cell_count
                this_cell = cell_list_a(c);
                X0 = this_cell.cell_loc.left;
                X1 = this_cell.cell_loc.right;
                Y0 = this_cell.cell_loc.top;
                Y1 = this_cell.cell_loc.bottom;
                
                cell_dbl = double(this_cell.get3DCellMask());
                cell_list_a(c).img_nobkg = immultiply(img_use_dbl(Y0:Y1,X0:X1,:), cell_dbl);
                cell_list_a(c).img_raw = immultiply(img_raw_dbl(Y0:Y1,X0:X1,:), cell_dbl);
                [cell_list_a(c).coord_list, ~] = this_cell.getCoordsSubset(quant_struct.t_coord_table);
            end
            clear img_raw_dbl;
            
            %Handle some common read-only values...
            xy_rad = quant_struct.spotzoom_r_xy;
            z_rad = quant_struct.spotzoom_r_z;
            z_adj = quant_struct.z_adj;
            xydim = (xy_rad * 2) + 1;
            do_gauss = quant_struct.do_gauss_fit;
            b_filters = RNAQuant.generateGaussianBFilters(xydim);
            
            task_count = cell_count;
            do_both = false;
            if quant_struct.do_gauss_fit & quant_struct.do_clouds
                %task_count = task_count * 2;
                do_both = true;
                %cell_list_b = cell_list_a;
            end
            %cell_list_b = cell_list_a;
            
            %Distribute Cloud detect and Gauss fit jobs to parallel workers
            pardir = char(RNAQuant.initParallel(1, quant_struct.workers, quant_struct.workdir));
            %pardir = char(quant_struct.workdir);
            parfor t = 1:task_count
                %tt = max(t - cell_count,1);
                if do_both
%                     if t > cell_count
%                         %This appears to be causing problems (maybe with the pre-compiling). Fix at some point.
%                         %Something wrong with indexing
%                         
%                         %Cloud
%                         my_cell = cell_list_b(tt);
%                         RNAQuant.DetectRNACloudsCellPar(my_cell, z_adj, pardir);
%                     else
%                         %Gauss fit
%                         my_cell = cell_list_a(t);
%                         RNAQuant.FitGaussians3CellPar(my_cell, xy_rad, z_rad, b_filters, pardir);
%                         RNAQuant.DetectRNACloudsCellPar(my_cell, z_adj, pardir);
%                     end

                    my_cell = cell_list_a(t);
                    RNAQuant.FitGaussians3CellPar(my_cell, xy_rad, z_rad, b_filters, pardir);
                    RNAQuant.DetectRNACloudsCellPar(my_cell, z_adj, pardir);
                    cell_list_a(t).img_nobkg = [];
                elseif do_gauss
                    %Just gaussian fit.
                    my_cell = cell_list_a(t);
                    RNAQuant.FitGaussians3CellPar(my_cell, xy_rad, z_rad, b_filters, pardir);
                    cell_list_a(t).img_nobkg = []; %Frees some memory
                else
                    %Just cloud
                    my_cell = cell_list_a(t);
                    RNAQuant.DetectRNACloudsCellPar(my_cell, z_adj, pardir);
                    cell_list_a(t).img_nobkg = [];
                    cell_list_a(t).coord_list = [];
                end
            end
            RNAQuant.shutdownParallel();
            
            %Cleanup par variables to free some memory.
            clear cell_list_a;
            %clear cell_list_b;
            clear b_filters;
            
            %Load back in gaussian fit data
            if quant_struct.do_gauss_fit
                for c = 1:cell_count
                    gfilepath = [pardir filesep sprintf('gaussians_%04d.mat', c)];
                    if isfile(gfilepath)
                        %Load spots
                        load(gfilepath, 'spots', 'nuc_ellip');
                        quant_struct.cell_rna_data(c).spots = spots;
                        quant_struct.cell_rna_data(c).nuc_ellip = nuc_ellip;
                        delete(gfilepath);
                    end
                end
            end
            
            %Load back in cloud data, and finish cloud processing
            min_cloud_vol = 0;
            if quant_struct.do_gauss_fit
                min_cloud_vol = RNAQuant.suggestMinCloudVolume(quant_struct.cell_rna_data);
            end
            if quant_struct.do_clouds
                for c = 1:cell_count
                    cfilepath = [pardir filesep sprintf('clouds_%04d.mat', c)];
                    
                    if isfile(cfilepath)
                        %Reload cell image.
                        my_cell = quant_struct.cell_rna_data(c);
                        cell_nobkg_dbl = my_cell.isolateCellBox(img_use_dbl);
                        
                        %Load and finish processing clouds
                        load(cfilepath, 'cell_cloud_mask', 'cloud_boxes');
                        quant_struct.cell_rna_data(c) = ...
                            RNAQuant.FinishRNACloudsCell(my_cell, cell_nobkg_dbl, cell_cloud_mask, cloud_boxes, min_cloud_vol);
                        delete(cfilepath);
                    end
                end
            end
            RNAQuant.removeParInternalDir(pardir);
            rmres = rmdir(pardir, 's');
            
            for c = 1:cell_count
                quant_struct.cell_rna_data(c) = quant_struct.cell_rna_data(c).updateSpotAndSignalValues();
            end
        end
        
        %%
        function quant_struct = FitRNA_S(quant_struct, img_clean, img_bkg)
            %Get some common stuff.
            xydim = (quant_struct.spotzoom_r_xy * 2) + 1;
            img_raw_dbl = double(img_clean);
            if quant_struct.no_bkg_subtract
                img_use_dbl = img_raw_dbl;
            else
                img_bkg_dbl = double(img_bkg);
                img_use_dbl = img_raw_dbl - img_bkg_dbl;
            end
            clear img_bkg_dbl;
            b_filters = RNAQuant.generateGaussianBFilters(xydim);
            
            %We can just do these cell by cell.
            cell_count = size(quant_struct.cell_rna_data,2);
            min_cloud_vol = 0;
            fprintf("[%s] RNAQuant.FitRNA_S -- Cells found: %d\n", datetime, cell_count);
            if quant_struct.do_gauss_fit
                fprintf("[%s] RNAQuant.FitRNA_S -- Fitting Gaussians to spots...\n", datetime);
                if ~isempty(quant_struct.cell_zero)
                    fprintf("[%s] RNAQuant.FitRNA_S -- Working on cell zero...\n", datetime);
                    my_cell = quant_struct.cell_zero;
                    cell_nobkg = my_cell.isolateCellBox(img_use_dbl);
                    
                    my_cell = RNAQuant.FitGaussians3Cell(cell_nobkg, my_cell,...
                        quant_struct.spotzoom_r_xy, quant_struct.spotzoom_r_z, b_filters);
                    quant_struct.cell_zero = my_cell;
                end

                for c = 1:cell_count
                    if (quant_struct.dbgcell > 0) & (c ~= quant_struct.dbgcell); continue; end
                    fprintf("[%s] RNAQuant.FitRNA_S -- Working on cell %d...\n", datetime, c);
                    my_cell = quant_struct.cell_rna_data(c);
                    cell_nobkg = my_cell.isolateCellBox(img_use_dbl);
                    
                    my_cell = RNAQuant.FitGaussians3Cell(cell_nobkg, my_cell,...
                        quant_struct.spotzoom_r_xy, quant_struct.spotzoom_r_z, b_filters);
                    quant_struct.cell_rna_data(c) = my_cell;
                end
                if quant_struct.do_clouds
                    min_cloud_vol = RNAQuant.suggestMinCloudVolume(quant_struct.cell_rna_data);
                end
            end
            
            if quant_struct.do_clouds
                fprintf("[%s] RNAQuant.FitRNA_S -- Finding clouds...\n", datetime);
                for c = 1:cell_count
                    if (quant_struct.dbgcell > 0) & (c ~= quant_struct.dbgcell); continue; end
                    my_cell = quant_struct.cell_rna_data(c);
                    cell_raw = my_cell.isolateCellBox(img_raw_dbl);
                    
                    [start_coords, ~] = my_cell.getCoordsSubset(quant_struct.t_coord_table);
                    if ~isempty(start_coords)
                        [cell_cloud_mask, cloud_boxes] = ...
                            RNAQuant.DetectRNACloudsCell(cell_raw, start_coords, quant_struct.z_adj);
                        clear cell_raw;
                        cell_nobkg = my_cell.isolateCellBox(img_use_dbl);
                        my_cell = RNAQuant.FinishRNACloudsCell(my_cell, cell_nobkg, cell_cloud_mask, cloud_boxes, min_cloud_vol);
                    else
                        clear cell_raw;
                    end

                    quant_struct.cell_rna_data(c) = my_cell;
                end
            end

            quant_struct = RNAQuant.updateQuantCounts(quant_struct, true, quant_struct.threshold);
        end
        
        %% ========================== Interface Internal ==========================
        
        %%
        function quant_output = FitRNA_WithRefilter(quant_input)
            quant_output = quant_input;
            clear quant_input;

            thList = quant_output.thresholds;
            if isempty(thList); return; end

            thList = sort(thList);
            thMin = thList(1);
            %Use this minimum to generate a coord table

            %1. Preprocess
            fprintf("[%s] RNAQuant.FitRNA_WithRefilter -- Running prefilter...\n", datetime);
            [img_clean, img_processed, imrmax_mask, img_bkg, quant_output.plane_stats] =...
                RNAQuant.RNAProcess3Dim(quant_output.img_raw, thMin, quant_output.small_obj_size,...
                quant_output.connect_size, quant_output.gaussian_radius, quant_output.workdir);

            %2. NonGauss RNA Pos
            fprintf("[%s] RNAQuant.FitRNA_WithRefilter -- Initial cellseg load and spot identification...\n", datetime);
            [quant_output.cell_rna_data, quant_output.cell_zero] = RNAQuant.NonGaussRNAPos(img_processed,...
                imrmax_mask, quant_output.cell_mask, quant_output.nuc_mask, quant_output.incl_cell_zero);
            clear imrmax_mask;

            %3. Gaussian Fit & Cloud Detection (where applicable)
            fprintf("[%s] RNAQuant.FitRNA_WithRefilter -- Now starting gaussian & cloud fitting...\n", datetime);
            if quant_output.workers > 1
                quant_output = RNAQuant.FitRNA_P(quant_output, img_clean, img_bkg);
            else
                quant_output = RNAQuant.FitRNA_S(quant_output, img_clean, img_bkg);
            end

            % Determine dropout thresholds (this is quite slow right now)
            cellCount = size(quant_output.cell_rna_data);
            if ~isempty(quant_output.cell_zero)
                myCell = quant_output.cell_zero;
                if isempty(obj.spotTable) & ~isempty(obj.spots)
                    myCell = myCell.convertSpotStorage();
                end
                myCell.spotTable{:, 'dropout_thresh'} = uint16(thMin);
                quant_output.cell_zero = myCell;
            end
            for c = 1:cellCount
                myCell = quant_output.cell_rna_data(c);
                if isempty(obj.spotTable) & ~isempty(obj.spots)
                    myCell = myCell.convertSpotStorage();
                end
                myCell.spotTable{:, 'dropout_thresh'} = uint16(thMin);
                quant_output.cell_rna_data(c) = myCell;
            end

            T = size(thList, 2);
            for t = 2:T
                thVal = thList(t);

                fprintf("[%s] RNAQuant.FitRNA_WithRefilter -- Finding dropouts at threshold: %d...\n", datetime, thVal);
                [~, imrmax_mask, ~] = ...
                    RNAQuant.FilterCleanedImage(img_clean, thVal, ...
                        quant_output.small_obj_size, quant_output.connect_size, quant_output.gaussian_radius);
                for c = 0:cellCount
                    if c > 0
                        myCell = quant_output.cell_rna_data(c);
                    else
                        myCell = quant_output.cell_zero;
                    end

                    imrmax_mask_cell = myCell.isolateCellBox(imrmax_mask);
                    ismax = find(imrmax_mask_cell);

                    boxSize = size(imrmax_mask_cell);
                    xx = myCell.spotTable{:, 'xinit'};
                    yy = myCell.spotTable{:, 'yinit'};
                    zz = myCell.spotTable{:, 'zinit'};
                    cellSpots1D = sub2ind(boxSize, yy, xx, zz);
                    foundAtTh = ismember(cellSpots1D, ismax);
                    myCell.spotTable{foundAtTh, 'dropout_thresh'} = thVal;

                    if c > 0
                        quant_output.cell_rna_data(c) = myCell;
                    else
                        quant_output.cell_zero = myCell;
                    end
                end
            end

        end

        %%
        function quant_output = FitRNA_WithoutRefilter(quant_input)
            quant_output = quant_input;
            clear quant_input;

            fprintf("[%s] RNAQuant.FitRNA_WithoutRefilter -- No refilter requested. Will try to use previous coordinates.\n", datetime);
            if ~isempty(quant_output.t_coord_table)
                fprintf('\tSpots in coord table: %d\n', size(quant_output.t_coord_table, 1));
                fprintf('\tCoord table includes dropout thresholds: %d\n', (size(quant_output.t_coord_table, 2) > 3));
            end

            %1. Preprocess
            fprintf("[%s] RNAQuant.FitRNA_WithoutRefilter -- Calculating background/cleaning dead pixels...\n", datetime);
            [img_clean, img_bkg] =...
                RNAQuant.RNAProcess3DimFast(quant_output.img_raw, quant_output.workdir);

            %2. Pos from table
            fprintf("[%s] RNAQuant.FitRNA_WithoutRefilter -- Initial cellseg load and spot identification...\n", datetime);
            [quant_output.cell_rna_data, quant_output.cell_zero] = RNAQuant.RNAPosFromTable(img_clean,...
                quant_output.t_coord_table, quant_output.cell_mask, quant_output.nuc_mask, quant_output.incl_cell_zero);

            %3. Gaussian Fit & Cloud Detection (where applicable)
            fprintf("[%s] RNAQuant.FitRNA_WithoutRefilter -- Now starting gaussian & cloud fitting...\n", datetime);
            if quant_output.workers > 1
                quant_output = RNAQuant.FitRNA_P(quant_output, img_clean, img_bkg);
            else
                quant_output = RNAQuant.FitRNA_S(quant_output, img_clean, img_bkg);
            end
        end

        %% ========================== Interface ==========================

        %%
        function quant_struct = updateQuantCounts(quant_struct, ignoreLikelyDups, useThreshold)
            if isempty(quant_struct); return; end

            cell_count = size(quant_struct.cell_rna_data, 2);
            sCount = 0;
            for c = 1:cell_count
                myCell = quant_struct.cell_rna_data(c);
                if isempty(myCell.spotTable) & ~isempty(myCell.spots)
                    myCell = myCell.convertSpotStorage();
                end

                if ~isempty(myCell.spotTable)
                    notCloud = ~myCell.spotTable{:, 'in_cloud'};
                    passTh = true(size(myCell.spotTable, 1), 1);
                    if useThreshold > 0
                        if ismember('dropout_thresh', myCell.spotTable.Properties.VariableNames)
                            passTh = myCell.spotTable{:, 'dropout_thresh'} >= useThreshold;
                        end
                    end
                    sCount = sCount + nnz(and(notCloud, passTh));
                end
                quant_struct.cell_rna_data(c) = myCell;
            end

            allInt = NaN(1, sCount);
            nowRow = 1;
            for c = 1:cell_count
                myCell = quant_struct.cell_rna_data(c);

                if ~isempty(myCell.spotTable)
                    notCloud = ~myCell.spotTable{:, 'in_cloud'};
                    passTh = true(size(myCell.spotTable, 1), 1);
                    if useThreshold > 0
                        if ismember('dropout_thresh', myCell.spotTable.Properties.VariableNames)
                            passTh = myCell.spotTable{:, 'dropout_thresh'} >= useThreshold;
                        end
                    end
                    cellInts = myCell.spotTable{and(notCloud, passTh), 'TotFitInt'};
                    edRow = nowRow + size(cellInts, 1) - 1;

                    allInt(nowRow:edRow) = cellInts';
                    nowRow = edRow + 1;
                end
            end
            clear myCell notCloud cellInts nowRow edRow

            iMean = mean(allInt, 'all', 'omitnan');
            iStd = std(allInt, 0, 'all', 'omitnan');
            globalBrightTh = iMean + (iStd * 2.0);
            globalSingleInt = median(allInt(allInt < globalBrightTh), 'all', 'omitnan');

            for c = 1:cell_count
                quant_struct.cell_rna_data(c) = ...
                    quant_struct.cell_rna_data(c).updateSpotAndSignalValues(globalBrightTh, globalSingleInt, ignoreLikelyDups, useThreshold);
            end

            if ~isempty(quant_struct.cell_zero)
                quant_struct.cell_zero = ...
                    quant_struct.cell_zero.updateSpotAndSignalValues(globalBrightTh, globalSingleInt, ignoreLikelyDups, useThreshold);
            end

            quant_struct.globalBrightTh = globalBrightTh;
            quant_struct.globalSingleInt = globalSingleInt;

            quant_struct.lastGlobalTestThreshold = useThreshold;
        end
        
        %%
        function quant_output = FitRNA(quant_input)
            if quant_input.do_refilter
                quant_output = RNAQuant.FitRNA_WithRefilter(quant_input);
            else
                quant_output = RNAQuant.FitRNA_WithoutRefilter(quant_input);
            end
        end
        
        %%
        function cell_data_table = cellData2Table(cell_rna_data)
            cell_data_table = table.empty();
            if isempty(cell_rna_data); return; end
            
            var_types = {'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' ...
                'uint32' 'uint32' 'uint32' 'double' 'double' 'double'};
            var_names = {'CellIndex' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'...
                'spots_total' 'spots_nuc' 'spots_cyto' 'signal_total' 'signal_nuc' 'signal_cyto'};
            field_count = size(var_names,2);
            
            cell_count = size(cell_rna_data,2);
            cell_data_table = table('Size', [cell_count field_count], 'VariableTypes', var_types, 'VariableNames', var_names);
            for c = 1:cell_count
                cell_rna_data(c) = cell_rna_data(c).updateSpotAndSignalValues();
                
                this_cell = cell_rna_data(c); %READ ONLY
                cell_data_table{c, 'CellIndex'} = this_cell.cell_number;
                cell_data_table{c, 'x0'} = this_cell.cell_loc.left;
                cell_data_table{c, 'x1'} = this_cell.cell_loc.right;
                cell_data_table{c, 'y0'} = this_cell.cell_loc.top;
                cell_data_table{c, 'y1'} = this_cell.cell_loc.bottom;
                cell_data_table{c, 'z0'} = this_cell.cell_loc.z_bottom;
                cell_data_table{c, 'z1'} = this_cell.cell_loc.z_top;
                
                %Save
                cell_data_table{c, 'spots_nuc'} = this_cell.spotcount_nuc;
                cell_data_table{c, 'spots_cyto'} = this_cell.spotcount_cyto;
                cell_data_table{c, 'spots_total'} = this_cell.spotcount_total;
                cell_data_table{c, 'signal_nuc'} = this_cell.signal_nuc;
                cell_data_table{c, 'signal_cyto'} = this_cell.signal_cyto;
                cell_data_table{c, 'signal_total'} = this_cell.signal_total;
            end
        end
        
        %%
        function results_pkg = results2SavePackage(quant_results)
            results_pkg = quant_results;
            %For now, only cell_rna_data is converted.
            if ~isempty(quant_results.cell_rna_data)
                myCells = arrayfun(@(mycell) mycell.packageForSave(), quant_results.cell_rna_data);
                results_pkg.cellData = myCells;
            else
                results_pkg.cellData = [];
            end

            if isfield(results_pkg, 'cell_zero') & ~isempty(quant_results.cell_zero)
                results_pkg.cellZeroData = quant_results.cell_zero.packageForSave();
            else
                results_pkg.cellZeroData = [];
            end

            results_pkg = rmfield(results_pkg, 'cell_rna_data');
            if isfield(results_pkg, 'cell_zero')
                results_pkg = rmfield(results_pkg, 'cell_zero');
            end
        end

        %%
        function quant_results = readResultsSavePackage(results_pkg)
            quant_results = results_pkg;

            if ~isempty(results_pkg.cellData)
                myCells = arrayfun(@(mycell) SingleCell.readFromSavePackage(mycell), results_pkg.cellData);
                quant_results.cell_rna_data = myCells;
            else
                quant_results.cell_rna_data = [];
            end

            if isfield(results_pkg, 'cellZeroData')
                if ~isempty(results_pkg.cellZeroData)
                    quant_results.cell_zero = SingleCell.readFromSavePackage(results_pkg.cellZeroData);
                else
                    quant_results.cell_zero = [];
                end
            else
                quant_results.cell_zero = [];
            end

            quant_results = rmfield(quant_results, 'cellData');
        end

    end
    
end