%
%Methods for RNA quantification
%Blythe Hospelhorn
%Modified from code written by Ben Kesler & Gregor Neuert
%Version 1.0.0
%Updated Oct 14, 2022

%Update Log:
%   1.0.0 | 22.10.14
%       Init Doc

classdef RNAQuant
    
    methods (Static)
        
        %%
        function param_struct = genGaussFitParamStruct()
            param_struct.mu1 = 0.0;
            param_struct.mu2 = 0.0;
            param_struct.s1 = 0.0;
            param_struct.s2 = 0.0;
            param_struct.A = 0.0;
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
        function g = Gauss3DFit(param_struct, xydim)
            mu1 = param_struct.mu1;
            mu2 = param_struct.mu2;
            s1 = param_struct.s1;
            A = param_struct.A;
            [X1,Y1] = meshgrid(1:1:xydim);
            g = A*exp(-((X1-mu1-1).^2)./(2*s1^2)-((Y1-mu2-1).^2)./(2*s1^2));
        end
        
        %%
        function g = Gauss3DFit_Vec(param_vec, xydim)
            mu1 = param_vec(1);
            mu2 = param_vec(2);
            s1 = param_vec(3);
            A = param_vec(4);
            [X1,Y1] = meshgrid(1:1:xydim);
            g = A*exp(-((X1-mu1-1).^2)./(2*s1^2)-((Y1-mu2-1).^2)./(2*s1^2));
        end
        
        %%
        function g = Gauss3DFitS1S2(param_struct, xydim)
            mu1 = param_struct.mu1;
            mu2 = param_struct.mu2;
            s1 = param_struct.s1;
            s2 = param_struct.s2;
            A = param_struct.A;
            [X1,Y1] = meshgrid(1:1:xydim);
            g = A*exp(-((X1-mu1-1).^2)./(2*s1^2)-((Y1-mu2-1).^2)./(2*s2^2));
        end
        
        %%
        function g = Gauss3DFitS1S2_Vec(param_vec, xydim)
            mu1 = param_vec(1);
            mu2 = param_vec(2);
            s1 = param_vec(3);
            s2 = param_vec(4);
            A = param_vec(5);
            [X1,Y1] = meshgrid(1:1:xydim);
            g = A*exp(-((X1-mu1-1).^2)./(2*s1^2)-((Y1-mu2-1).^2)./(2*s2^2));
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
        
        %%
        function [img_processed, imrmax_mask, img_bkg, plane_stats] = RNAProcess3Dim(img_raw, threshold, small_obj_size, connect_size, gaussian_radius, workdir)
            %B2_RNAprocess3Dim.m
            img_processed = [];
            
            if isempty(img_raw) 
                plane_stats = [];
                return; 
            end
            
            if nargin < 6; workdir = pwd; end
            
            dead_pix_path = [workdir filesep 'deadpix.mat'];

            %Pre-filtering...
            IMG3D = uint16(img_raw);
            clear img_raw;
            RNA_Threshold_Common.saveDeadPixels(IMG3D, dead_pix_path);
            IMG3D = RNA_Threshold_Common.cleanDeadPixels(IMG3D, dead_pix_path);
            
            if isfile(dead_pix_path)
                delete(dead_pix_path);
            end

            if nargin < 5; gaussian_radius = 7; end
            if nargin < 3; small_obj_size = 3; end
            if nargin < 4; connect_size = 8; end
            
            X = size(IMG3D,2);
            Y = size(IMG3D,1);
            Z = size(IMG3D,3);
            plane_stats(Z) = struct('stdev', 0.0, 'mean', 0.0, 'variance', 0.0, 'median', 0.0);
            
            img_processed = uint16(zeros(Y,X,Z));
            img_bkg = uint16(zeros(Y,X,Z));
            for z = 1:Z
                z_raw = IMG3D(:,:,z);
                z_filtered = RNA_Threshold_Common.applyGaussianFilter(z_raw, gaussian_radius, 2);
                z_filtered = RNA_Threshold_Common.applyEdgeDetectFilter(z_filtered);
                
                z_thresh = z_filtered > threshold;
                z_thresh = bwareaopen(z_thresh, small_obj_size, connect_size);
                img_processed(:,:,z) = immultiply(z_thresh, z_raw);
                
                img_bkg(:,:,z) = medfilt2(z_raw, [20 20]); 
                
                plane_stats(z).stdev = std2(z_filtered(:));
                plane_stats(z).mean = mean2(z_filtered(:));
                plane_stats(z).variance = var(double(z_filtered(:)));
                plane_stats(z).median = median(double(z_filtered(:)));
            end
            
            img_processed = RNA_Threshold_Common.blackoutBorders(img_processed, gaussian_radius, 0);
            imrmax_mask = imregionalmax(img_processed,26);
        end
        
        %%
        function [cell_max_rna, rna_pos] = NonGaussRNAPos(img_raw, imrmax_mask, cell_mask, nuc_mask)
            %B3_nongaussRNApos

            Y = size(img_raw,1);
            X = size(img_raw,2);
            Z = size(img_raw,3);
            cell_count = max(cell_mask(:));
            
            rna_pos = cell(1, cell_count);
            cell_max_rna(cell_count) = struct('cell_rna', 0, 'nuc_rna', 0, 'cyto_rna', 0);
            
            img_imrmax = immultiply(img_raw, imrmax_mask);
            for c = 1:cell_count
                cell_filter = uint16(cell_mask == c);
                cell_box = regionprops(cell_filter,'BoundingBox','Area').BoundingBox;
                
                X0 = round(cell_box(1)) - 4;
                Y0 = round(cell_box(2)) - 4;
                X1 = round(cell_box(1)+ cell_box(3)) + 4;
                Y1 = round(cell_box(2)+ cell_box(4)) + 4;
                
                if X0 < 1; X0 = 1; end
                if Y0 < 1; Y0 = 1; end
                if X1 > X; X1 = X; end
                if Y1 > Y; Y1 = Y; end
                box_filter = cell_filter(Y0:Y1,X0:X1);
                
                CY = Y1-Y0+1;
                CX = X1-X0+1;
                
                cell_rna_imrmax = NaN(CY,CX,Z);
                cell_nuc = NaN(CY,CX,Z);
                cell_cyto = NaN(CY,CX,Z);
                cell_dbl = NaN(CY,CX,Z);
                for z = 1:Z
                    cell_rna_imrmax(:,:,z) = immultiply(img_imrmax(Y0:Y1,X0:X1,z), box_filter);
                    cell_nuc(:,:,z) = immultiply(double(nuc_mask(Y0:Y1,X0:X1,z)), double(box_filter));
                    cell_cyto(:,:,z) = imsubtract(double(cell_filter(Y0:Y1,X0:X1)),double(cell_nuc(:,:,z)));
                    cell_dbl(:,:,z) = double(cell_filter(Y0:Y1,X0:X1)); 
                end
                
                cell_cyto = immultiply(double(cell_rna_imrmax),double(cell_cyto)); 
                cell_nuc = immultiply(double(cell_rna_imrmax),double(cell_nuc)); 
                cell_rna_imrmax = immultiply(double(cell_rna_imrmax),cell_dbl); 
                
                [lbl_cell, numr_cell] = bwlabeln(cell_rna_imrmax > 0,26); 
                [~, numr_cyto] = bwlabeln(cell_cyto > 0,26); 
                [~, numr_nuc] = bwlabeln(cell_nuc > 0,26); 
                
                if numr_cell > 0
                    pos_cell(numr_cell) = struct('x', 0, 'y', 0, 'z', 0);
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
                        
                        pos_cell(i).x = x;
                        pos_cell(i).y = y;
                        pos_cell(i).z = z;
                    end
                    rna_pos{c} = pos_cell;
                end
                cell_max_rna(c).cell_rna = numr_cell;
                cell_max_rna(c).nuc_rna = numr_nuc;
                cell_max_rna(c).cyto_rna = numr_cyto;
            end
        end
        
        %%
        function param_struct = genGaussFitResultStruct()
            param_struct.xfit = NaN;
            param_struct.yfit = NaN;
            param_struct.xinit = NaN;
            param_struct.yinit = NaN;
            param_struct.xgw = NaN;
            param_struct.ygw = NaN;
            param_struct.xFWHM = NaN;
            param_struct.yFWHM = NaN;
            param_struct.expMInt = NaN;
            param_struct.fitMInt = NaN;
            param_struct.TotExpInt = NaN;
            param_struct.TotFitInt = NaN;
            param_struct.r = NaN;
            param_struct.rFit = NaN;
            param_struct.back = NaN;
            param_struct.xsem = NaN;
            param_struct.ysem = NaN;
            param_struct.zxfit = NaN;
            param_struct.zyfit = NaN;
            param_struct.zinit = NaN;
            param_struct.zint = NaN;
            param_struct.zrel = NaN;
            param_struct.zabs = NaN;
            param_struct.zstd = NaN;
            param_struct.zqFit = NaN;
            param_struct.nucRNA = false;
            param_struct.distRNA = NaN;
            param_struct.szNUC = NaN;
            param_struct.cytoRNA = false;
            param_struct.xabsloc = NaN;
            param_struct.yabsloc = NaN;
            param_struct.normdist = NaN;
            param_struct.nucR = NaN;
        end
        
        %%
        function fit_result = FitGaussians3(img_raw, img_bkg, rna_pos, cell_mask, nuc_mask, xy_rad)
            
            if nargin < 6
                xy_rad = 4;
            end
            z_rad = 2;
            
            Y = size(img_raw,1);
            X = size(img_raw,2);
            Z = size(img_raw,3);
            cell_count = max(cell_mask(:));
            
            %Pixel border filters
            xydim = (xy_rad * 2) + 1;
            zdim = (z_rad * 2) + 1;
            b_filter_A = double(ones(xydim,xydim));
            b_filter_A(2:xydim-1,2:xydim-1) = 0;
            b_filter_B = (1.0 - b_filter_A) * 0.5;
            b_filter_C = b_filter_B;
            b_filter_C(3:xydim-2,3:xydim-2) = 2.0;
            b_filter_C = b_filter_C ./ 2.0;

            for c = cell_count:-1:1
                %Alloc output
                cell_spots = rna_pos{c};
                spot_count = size(cell_spots);
                cell_fits(spot_count) = RNAQuant.genGaussFitResultStruct();
                
                %Filter down to just this cell
                cell_filter = uint16(cell_mask == c);
                cell_box = regionprops(cell_filter,'BoundingBox','Area').BoundingBox;
                
                X0 = round(cell_box(1)) - 4;
                Y0 = round(cell_box(2)) - 4;
                X1 = round(cell_box(1)+ cell_box(3)) + 4;
                Y1 = round(cell_box(2)+ cell_box(4)) + 4;
                
                if X0 < 1; X0 = 1; end
                if Y0 < 1; Y0 = 1; end
                if X1 > X; X1 = X; end
                if Y1 > Y; Y1 = Y; end
                box_filter = cell_filter(Y0:Y1,X0:X1);
                
                for z = 1:Z       
                    cell_dbl(:,:,z) = double(box_filter);
                    cell_rna_raw(:,:,z) = immultiply(img_raw(Y0:Y1,X0:X1,z),box_filter);
                    cell_rna_bkg(:,:,z) = immultiply(img_bkg(Y0:Y1,X0:X1,z),box_filter);
                    cell_rna_nuc(:,:,z) = immultiply(double(nuc_mask(Y0:Y1,X0:X1,z)), cell_dbl);
                end
                
                cell_rna_raw = immultiply(double(cell_rna_raw),cell_dbl);
                cell_rna_bkg = immultiply(double(cell_rna_bkg),cell_dbl);
                
                if spot_count > 0
                    for s = 1:spot_count
                        my_spot = cell_spots(s);
                        
                        %Isolate the pixels immediately around the spot
                        spot_data = RNAUtils.isolateSpotData(cell_rna_raw, my_spot.x, my_spot.y, my_spot.z, xy_rad, z_rad);
                        
                        spot_bkg(zdim) = 0.0;
                        spot_filtered = NaN(xydim, xydim, zdim);
                        for z = 1:zdim
                            spot_bkg(z) = sum(sum(immultiply(spot_data(:,:,z),double(b_filter_A)),1))./sum(b_filter_A(:));
                            spot_filtered(:,:,z) = immultiply(spot_data(:,:,z),double(b_filter_B)); %Not used? Is this a bug?
                            spot_filtered(:,:,z) = round(immultiply(spot_data(:,:,z),double(b_filter_C))./1);
                        end
                        
                        %Try to fit Gaussian to spot
                        gparams_vec = RNAQuant.genGaussFitParamVector(2, 3, 1, double(max(spot_filtered(:))));
                        [mxX, mxY] = meshgrid(1:1:xydim);
                        options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9);
                        
                        axis_ratios(zdim) =  0.0;
                        gauss_sim = NaN(xydim, xydim, zdim);
                        fitted_params = NaN(5,zdim);
                        for z = 1:zdim
                            slice_filtered = spot_filtered(:,:,z);
                            slice_scaled = TMRimm./max(slice_filtered(:)); 
                            slice_scaled = im2bw(slice_scaled,0.5);
                            iprops = regionprops(double(slice_scaled),'MinorAxisLength','MajorAxisLength');
                            axis_ratios(z) = iprops.MajorAxisLength ./ iprops.MinorAxisLength; %try/catch?
                            
                            lb = [2 2 gparams_vec(3) gparams_vec(4)];
                            ub = [xydim-2 xydim-2  gparams_vec(3) 5*gparams_vec(4)] + 0.001;
                            [resparams,~,~,~,~,~,~] = lsqcurvefit(@RNAQuant.Gauss3DFit_Vec, gparams_vec, xydim, slice_scaled, lb, ub, options);
                            
                            lb = [resparams(1) resparams(2) 0 0];
                            ub = [resparams(1) resparams(2) 5 2*resparams(4)] + 0.001;
                            [resparams,~,~,~,~,~,~] = lsqcurvefit(@RNAQuant.Gauss3DFit_Vec, gparams_vec, xydim, slice_scaled, lb, ub, options);
                            
                            lb = [resparams(1) resparams(2) 0 0 0];
                            ub = [resparams(1) resparams(2) resparams(3)*3 resparams(3)*3 resparams(4)*3] + 0.001;
                            resparams([1 2 3 5]) = resparams;
                            resparams(4) = resparams(3);
                            [resparams,~,~,~,~,~,~] = lsqcurvefit(@RNAQuant.Gauss3DFitS1S2_Vec, resparams, xydim, slice_scaled, lb, ub, options);
                            
                            gauss_sim(:,:,z) = resparams(5)*exp(-((mxX-resparams(1)).^2)./(2*resparams(3)^2)-((mxY-resparams(2)).^2)./(2*resparams(4)^2));
                            fitted_params(:,z) = resparams;
                            fitted_params(1) = fitted_params(1) + my_spot.x - xy_rad - 1;
                            fitted_params(2) = fitted_params(2) + my_spot.y - xy_rad - 1;
                        end
                        
                        %Save some data to the return struct
                        spot_params = cell_fits(s);
                        spot_params.xfit = fitted_params(1,3);
                        spot_params.yfit = fitted_params(2,3);
                        spot_params.xinit = my_spot.x;
                        spot_params.yinit = my_spot.y;
                        spot_params.xgw = fitted_params(3,3);
                        spot_params.ygw = fitted_params(4,3);
                        spot_params.xFWHM = 2*sqrt(2*log(2)).*fitted_params(3,3);
                        spot_params.yFWHM = 2*sqrt(2*log(2)).*fitted_params(4,3);
                        spot_params.expMInt = double(max(spot_filtered(:)));
                        spot_params.fitMInt = fitted_params(5,3);
                        spot_params.TotExpInt = nansum(nansum(spot_filtered(:,:,3)));
                        spot_params.TotFitInt = nansum(nansum(gauss_sim(:,:,3)));
                        spot_params.r = axis_ratios(3);
                        if spot_params.xFWHM > spot_params.yFWHM
                            spot_params.rFit = spot_params.xFWHM./spot_params.yFWHM;
                        else
                            spot_params.rFit = spot_params.yFWHM./spot_params.xFWHM;
                        end
                        spot_params.back = spot_bkg(3);
                        
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
                        [~,minidx] = min(fitted_params(3,:));
                        if (minidx > 2 & minidx < 5)
                            df_params = double([1 3 2 0.67 0.3]);
                            lb = [0 0 0 0 0]; % lower bounds for finding the gaussian width
                            ub = [10 10 10 10 10] + 0.001; % upper bounds for finding the gaussian width
                            z_vec = 1:5;
                            options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9); % fit conditions
                            [df_params] = lsqcurvefit(@RNAQuant.Defoc4Order_Vec, df_params, z_vec, fitted_params(3,:), lb, ub, options);
                            z_vec = [1:0.01:5];
                            tempvec = df_params(1) .* sqrt(1 + ((z_vec - df_params(2)) ./ df_params(3)).^2 + ...
                                df_params(4) .* ((z_vec - df_params(2)) ./ df_params(3)).^3 + ...
                                df_params(5) .* ((z_vec - df_params(2)) ./ df_params(3)).^4 );
                            [~,minidx] = min(tempvec);
                            spot_params.zxfit = z_vec(minidx); %gm1
                        else
                            spot_params.zxfit = NaN; %gm1
                        end
                        
                        %2. Using y width
                        [~,minidx] = min(fitted_params(4,:));
                        if (minidx > 2 & minidx < 5)
                            df_params = double([1 3 2 0.67 0.3]);
                            lb = [0 0 0 0 0]; % lower bounds for finding the gaussian width
                            ub = [10 10 10 10 10] + 0.001; % upper bounds for finding the gaussian width
                            z_vec = 1:5;
                            options = optimset('Display', 'none','MaxFunEvals', 200,'MaxIter', 200,'TolX', 1E-9); % fit conditions
                            [df_params] = lsqcurvefit(@RNAQuant.Defoc4Order_Vec, df_params, z_vec, fitted_params(4,:), lb, ub, options);
                            z_vec = [1:0.01:5];
                            tempvec = df_params(1) .* sqrt(1 + ((z_vec - df_params(2)) ./ df_params(3)).^2 + ...
                                df_params(4) .* ((z_vec - df_params(2)) ./ df_params(3)).^3 + ...
                                df_params(5) .* ((z_vec - df_params(2)) ./ df_params(3)).^4 );
                            [~,minidx] = min(tempvec);
                            spot_params.zyfit = z_vec(minidx); %gm2
                        else
                            spot_params.zyfit = NaN; %gm2
                        end
                        
                        %3. Using intensity
                        [~,minidx] = min(fitted_params(5,:));
                        if (minidx > 2 & minidx < 5)
                            p3 = polyfit([1:5], fitted_params(5,:), 4);
                            tempvec = 1:0.01:5;
                            [~,minidx] = max(polyval(p3,tempvec));
                            spot_params.zqFit = max(diff(diff(polyval(p3,tempvec),2)));
                            spot_params.zint = tempvec(minidx);
                            spot_params.zrel = nanmean([spot_params.zxfit spot_params.zyfit spot_params.zint]);
                            spot_params.zstd = nanstd([spot_params.zxfit spot_params.zyfit spot_params.zint]);
                        else
                            spot_params.zint = NaN; %gm3
                            spot_params.zrel = NaN; %gm
                            spot_params.zstd = NaN; %gms
                            spot_params.zqFit = NaN; %r1
                        end
                        
                        spot_params.zinit = my_spot.z;
                        spot_params.zabs = spot_params.zrel + my_spot.z - 3;
                        
                        %Check if nuclear RNA, and calculate related
                        %properties.
                        spot_params.nucRNA = cell_rna_nuc(round(fitted_params(1,3)), round(fitted_params(2,3)), my_spot.z) > 0;
                        if spot_params.nucRNA
                        end
                        
                        cell_fits(s) = spot_params;
                    end
                else
                    fit_result{c} = [];
                    continue;
                end
                
                %Save output
                fit_result{c} = cell_fits;
            end
            
        end
        
    end
    
end