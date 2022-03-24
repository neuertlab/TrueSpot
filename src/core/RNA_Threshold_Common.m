%Common functions for RNA thresholding
%Blythe Hospelhorn
%Modified from code written by Ben Kesler & Gregor Neuert
%Version 2.2.2
%Updated September 9, 2021

%Modified from ABs_Threshold3Dim
%Copied from bgh_3DThresh_Common

%Update Log:
%   1.0.0 | 19.09.10
%
%   1.1.0 | 19.09.19
%       Fixed bug in unique spot detection method
%   1.1.1 | 19.09.20
%       Adapted to new file, tweaked some message outputs
%   1.2.0 | 19.09.27
%       Tweaked functions to ignore zeros in image (for masked bkg)
%       Also added funcs for getting stdev and average ignoring zeros
%   1.3.0 | 19.10.24
%       Added mask spot count scaling functions
%       More comments

%   2.0.0 | 21.02.18
%       Finally updated documentation...
%       Quite a few debugs and changes since 1.3.0
%   2.1.0 | 21.07.19
%       Added alternative to imhistc finally 
%       Background mask coord table filter was omitting z coords, lovely
%   2.2.0 | 21.08.17
%       Thresholder function added
%       Changed save format for dead pixel detection
%   2.2.1 | 21.08.31
%       mask_spots now takes either a 2D or 3D mask
%   2.2.2 | 21.09.09
%       Kept debug plot for window score in auto threshold method, moved to
%       its own method.
%       

%%
%
classdef RNA_Threshold_Common
    methods (Static)
        
        %%-------------------------[Filtering]-----------------------------
        
        %%
        %Generate a Gaussian filter matrix
        %ARGS
        %   fs (int) - ??
        %
        %RETURN
        %   gauss_filter 
        %
        function gauss_filter = getGaussianFilter(fs)
            %Set some variable values...
            mu1 = 7;
            mu2 = 7;
            s1 = 2;
    
            [X1,Y1] = meshgrid(1:1:fs);
            gauss_filter = 10 * exp(-((X1-mu1).^2)./(2*s1^2)-((Y1-mu2).^2)./(2*s1^2));
            gauss_filter = gauss_filter./sum(gauss_filter(:));
        end

        %%
        %Apply a gaussian filter to a 3D image matrix
        %ARGS
        %   in_img (num[Y][X][Z]) - Input image
        %   fs (int) - ??
        %
        %RETURN
        %   img_filtered (num[Y][X][Z]) - Filtered image
        %
        %ASSUMPTIONS
        %   -> The image matrix is 3D
        %
        %FLEXIBILITIES
        %   -> "num" can refer to any numerical type
        %
        function img_filtered = applyGaussianFilter(in_img,fs)
            %Get the filter matrix
            filter_mx = RNA_Threshold_Common.getGaussianFilter(fs);
    
            %figure(10101);
            %imshow(filter_mx, []);
    
            img_filtered = in_img;
    
            %Do slice by slice
            Z = size(in_img,3);
            for k=1:Z
                f_slice = imfilter(in_img(:,:,k),filter_mx);
                img_filtered(:,:,k) = f_slice;
            end

        end

        %%
        %Apply an edge detection filter to a 3D image matrix
        %ARGS
        %   in_img (num[Y][X][Z]) - Input image
        %
        %RETURN
        %   img_filtered (num[Y][X][Z]) - Filtered image
        %
        %ASSUMPTIONS
        %   -> The image matrix is 3D
        %
        %FLEXIBILITIES
        %   -> "num" can refer to any numerical type
        %
        function img_filtered = applyEdgeDetectFilter(in_img)
            %Get the filter matrix
            filter_mx = [-1 -1 -1;...
                         -1 +8 -1;...
                         -1 -1 -1];
    
            img_filtered = in_img;
    
            %Do slice by slice
            Z = size(in_img,3);
            for k=1:Z
                f_slice = imfilter(in_img(:,:,k),filter_mx);
                img_filtered(:,:,k) = f_slice;
            end
    
        end
        
        %%
        %-------------------------[Cleanup]-----------------------------
        
        %%
        %Set the borders of each 2D frame in a 3D image to black
        %ARGS
        %   in_img (uint16[Y][X][Z]) - Input image
        %   border_width (int) - Number of pixels to trim (blackout) from each side
        %                       of the image
        %   trimZ(bool) - Whether or not to also trim in the Z direction (zeros-out
        %                   planes 1 to border and Z-border to Z)
        %
        %RETURN
        %   img_filtered (uint16[Y][X][Z]) - Filtered image
        %
        %ASSUMPTIONS
        %   -> The image matrix is 3D
        %
        %FLEXIBILITIES
        %   -> "int" can refer to any integer type
        %
        function img_filtered = blackoutBorders(in_img,border_width,trimZ)
    
            img_filtered = in_img;
    
            if border_width > 0
                %Do slice by slice
                X = size(in_img,2);
                Y = size(in_img,1);
                Z = size(in_img,3);
        
                for k=1:Z
                    if trimZ
                        if (k <= border_width) || (k >= (Z - border_width))
                            img_filtered(:,:,k) = zeros(Y,X);
                            continue;
                        end
                    end
                    slice = in_img(:,:,k);
                    slice(1:border_width,:) = 0; %Top
                    slice(:,1:border_width) = 0; %Left
                    slice(Y-border_width+1:Y,:) = 0; %Bottom
                    slice(:,X-border_width+1:X) = 0; %Right
                    img_filtered(:,:,k) = slice;
                end
            end
    
        end

        %%
        %Detect dead pixels in an image and save their coordinates (1D) to
        %a file in the workspace for later reloading.
        %ARGS
        %   in_img (num[Y][X][Z]) - Input image
        %   savepath (string) - Path to save dead pixel coords to
        %       Defaults to 'recurring pixels 3 out of 6' if no/empty arg
        %
        %RETURN
        %   px_counts (int[3][1]) - Random dead seeming pixel count(2),
        %                           detected dead pixel count (3)
        %
        %SAVES
        %   'recurring pixels 3 out of 6' (to workspace)
        %       List of 1D (down cols) coordinates of dead pixels
        %
        %ASSUMPTIONS
        %   -> The image matrix is 3D
        %
        %FLEXIBILITIES
        %   -> "num" can refer to any numerical type
        %
        function px_counts = saveDeadPixels(in_img, savepath)
        %Copied wholesale from AB_FindThreshold_TMR_AF594_CY5_B.m

            if isempty(savepath)
                savepath = 'recurring pixels 3 out of 6';
            end
        
            Y = size(in_img,1);
            X = size(in_img,2);

            %testmat = zeros(3,1);
            %counter_test = 1;
            px_counts = zeros(3,1);
            slice_cut = 4/8;
            slice_check = 6;
            counter99 = 0;
            hi_pixels_all = zeros(1,1);
            %'determining recurring pixels'
            fprintf("Determining recurring pixels...\n")
    
            w2 = [-1 -1 -1;...
                  -1 +8 -1;...
                  -1 -1 -1;];
      
            if size(in_img,1) > 1
                for i = 1:slice_check
                    counter99 = counter99 + 1;
                    slice = imfilter(in_img(:,:,i),w2);
                    cutoff = mean(slice(:)) + 3 * std(slice(:));
                    temp_pix = find(slice > cutoff);
                    hi_pixels_all(counter99,1:size(temp_pix,1)) = temp_pix;
                end
            end
    
            slice_num = size(hi_pixels_all,1);
            rand_hi = zeros(size(hi_pixels_all));
    
            for j = 1:slice_num
                temp_hi = hi_pixels_all(j,:);
                temp_hi1 = temp_hi;
                temp_hi1(temp_hi1 == 0) = [];
                rand_hi_temp = randsample(Y * X, size(temp_hi1,2));
                rand_hi_temp = rand_hi_temp';
                rand_hi(j, 1:size(rand_hi_temp,2)) = rand_hi_temp;
            end
            
            imhistc_exists = exist('imhistc');
            if imhistc_exists == 3
                %this is a binary mex file - the one I have is only gonna
                %run on win x64 I think
                fprintf("imhistc was found. Using imhistc.\n");
                hist_pix_rand = imhistc(double(rand_hi(:)), double(Y*X), 0, double(Y*X));
            else
                fprintf("imhistc was NOT found. Using imhist.\n");
                hist_pix_rand = imhist(double(rand_hi(:)), double(Y*X));
            end
            %hist_pix_rand = imhist(double(rand_hi(:)), double(Y*Y), 0, double(Y*Y));
            rand_recur_pix = find(hist_pix_rand >= slice_num*slice_cut);
            %testmat(2, counter_test) = size(rand_recur_pix,1);
            px_counts(2, 1) = size(rand_recur_pix,1);
            %[num2str(size(rand_recur_pix,1)) ' randomly recurring pixels' ]
            fprintf("%d randomly recurring pixels\n", size(rand_recur_pix,1));
            hi_pixels_all(hi_pixels_all == 0) = [];
            if imhistc_exists == 3
                hist_pix = imhistc(double(hi_pixels_all(:)), double(Y * X), 0, double(Y * X));
            else
                hist_pix = imhist(double(hi_pixels_all(:)), double(Y * X));
            end
            recurring_pixels = find(hist_pix >= slice_num*slice_cut);
            %[num2str(size(recurring_pixels,1)) ' recurring pixels' ]
            fprintf("%d recurring pixels\n", size(recurring_pixels,1));
            %testmat(3, counter_test) = size(recurring_pixels,1);
            px_counts(3, 1) = size(recurring_pixels,1);
            %counter_test = counter_test + 1;
            
            %Convert back to 2D coordinates to save
            
            save(savepath, 'recurring_pixels')
   
        end

        %%
        %Clean dead pixels from image (provided list has been previously saved) by
        %resetting intensity values of dead pixels to average of surrounding
        %pixels. The function ignores pixels that fall on the image border.
        %ARGS
        %   in_img (num[Y][X][Z]) - Input image
        %
        %RETURN
        %   clean_img (num[Y][X][Z]) - Filtered image
        %
        %ASSUMPTIONS
        %   -> That a list of dead pixels (in 1D coordinates going down columns)
        %   has been previously saved to "recurring pixels 3 out of 6"
        %   -> The image matrix is 3D
        %
        %FLEXIBILITIES
        %   -> "num" can refer to any numerical type
        %
        function clean_img = cleanDeadPixels(in_img, savepath)

            if isempty(savepath)
                savepath = 'recurring pixels 3 out of 6';
            end
            
            load(savepath, 'recurring_pixels') %Load previously saved list of dead pixels
            rp_count = size(recurring_pixels,1);
    
            %Get size(s) and save. Cleans up code.
            dim1 = size(in_img,1); %Height
            dim2 = size(in_img,2); %Width
            dim3 = size(in_img,3); %Depth
    
            clean_img = in_img;
    
            %exmat is a matrix containing indices of border coordinates
            %Border pixels are excluded from cleaning
            exmat = zeros(max(dim1,dim2),4);
            sz = dim1 * dim2;
            exmat(1:dim1,1) = 1:dim1; %Left
            exmat(1:dim1,4) = (sz-dim1+1):sz; %Right
            exmat(1:2,2:3) = 1;
            exmat(3:dim1,2) = dim1*(2:(dim2-1)); %Bottom
            exmat(3:dim1,3) = dim1*(1:(dim2-2))+1; %Top
    
            %The notice messages are kept from previous code
            %'Averaging recurring pixels'
            fprintf("Averaging recurring pixels...\t")
            tic
            for k = 1:dim3
                slice = in_img(:,:,k); %[uint16[D][D]] 2D slice of image
                for c = 1:rp_count
                    p = recurring_pixels(c); %[int] 1D coordinate of bad pixel
                    exclude_check = (exmat == p); %[bool[D][D]] Whether each index in exmat == bad pix val
                    if max(exclude_check(:)) == 0 %No matches in exclusion matrix
                        west = p - dim1;
                        east = p + dim1;
                        surr_pix = [p-1, p+1, west, east, west-1, east-1, east+1, west+1];
                        %N S W E NW NE SE SW
                        slice(p) = mean(slice(surr_pix));
                    end
                end
                clean_img(:,:,k) = slice;
            end
            toc
    
            clear recurring_pixels
        end
        
        %%
        %-------------------------[Spot Detection]-----------------------------
        
        %%
        % (Method writing in progress)
        % Examine the derivative of a plot of spots detected vs. threshold 
        % to determine the best threshold to use for this image.
        % (deriv features tbd)
        %
        %ARGS
        %   spotcount_table (num[T][2]) - Table of counts of spots detected
        %       at T different threshold values. First column should be
        %       threshold, second column should be # of spots detected at
        %       that threshold.
        %   bkg_spotcount_table (num[T][2]) - OPTIONAL - May pass empty mtx
        %       Spot count table as above, but for background or control
        %   window_size (int) - Size of scan window (number of threshold
        %       points to consider at once)
        %   windowpos (double) - Position relative to scan window to call
        %       threshold when good range is found.
        %       0.0 - Back of the window (lowest threshold, most
        %           conservative)
        %       0.5 - Middle of the window
        %       1.0 - Front of the window (highest threshold)
        %   stdev_thresh (double) - Value window stdev must fall below to
        %       be called threshold. Default: 1.0
        %
        %RETURN
        %   threshold (int) - Suggested threshold derived from shape of
        %       spot count plot
        %
        function [threshold, win_out] = estimateThreshold(spotcount_table, bkg_spotcount_table, window_size, windowpos, stdev_thresh)
            
            deriv1 = diff(spotcount_table(:,2));
            deriv1 = smooth(deriv1);
            deriv1 = abs(deriv1);
            
            if ~isempty(bkg_spotcount_table)
                bkgderiv = diff(bkg_spotcount_table(:,2));
                bkgderiv = smooth(bkgderiv);
                bkgderiv = abs(bkgderiv);
            else
                bkgderiv = [];
            end
            
            %Find global max for bkg if present or data if not
            %Scan start point can only be AFTER this
            maxidx = -1;
            maxvalue = 0;
            checkcurve = deriv1;
            if ~isempty(bkgderiv)
                checkcurve = bkgderiv;
            end
            pointcount = size(checkcurve,1);
            for i = 1:pointcount
                if checkcurve(i) > maxvalue
                    maxidx = i;
                    maxvalue = checkcurve(i);
                end
            end
            
            %If bkg is present, find start point (when bkg deriv first = 0)
            %If it happens.
            startidx = maxidx;
            if ~isempty(bkgderiv)
                for i = maxidx:pointcount
                    if bkgderiv(i) < 1
                        startidx = i;
                        break;
                    end
                end
            end
            fprintf("Starting threshold scan from threshold %d\n", spotcount_table(startidx, 1));
            
            %Adjust starting point for window position
            winshift = 0;
            if windowpos > 0.0
                if windowpos >= 1.0
                    winshift = window_size;
                else
                    winshift = round(window_size * windowpos);
                end
                startidx = startidx - winshift;
            end
            
            
            %Window stdev
            %winstdev = NaN(pointcount, 1);
            winout = NaN(pointcount, 1);
            winmax = pointcount - window_size;
            for i = startidx:winmax
                w_back = i;
                w_front = i + window_size - 1;
                %winout(i) = std(deriv1(w_back:w_front,1));
                winout(i) = std(deriv1(w_back:w_front,1)) / mean(deriv1(w_back:w_front,1));
            end
            
            %Shift
            %win_stdevs = NaN(pointcount, 1);
            win_out = NaN(pointcount, 1);
            for i = 1:(pointcount - winshift)
                %win_stdevs(i+winshift) = winstdev(i);
                win_out(i+winshift) = winout(i);
            end
            
            %DEBUG plots
            %x = spotcount_table(1:pointcount,1);
            %y = winstdev;
            %xlinecolor = [0.608 0.621 0.628];
            %ylinecolor = [0.608 0.621 0.628];
            %handle1 = figure(88888);
            %ax = axes;
            %plot(x, y, 'LineWidth', 2, 'Color', 'k');
            %ax.FontSize = 14;
            %hold on;
            %line(get(ax,'XLim'), [stdev_thresh stdev_thresh],'Color',xlinecolor,'LineStyle','--');
            %xlabel('Threshold', 'FontSize', 16);
            %ylabel('Window Score', 'FontSize', 16);
            %title(['Window Size = ' num2str(window_size)], 'FontSize', 18);
            
            %figure(88889);
            %x = x(startidx:(pointcount-1));
            %y = deriv1(startidx:(pointcount-1));
            %ax2 = axes;
            %plot(x, y, 'LineWidth', 2, 'Color', 'k');
            %ax2.FontSize = 14;
            %hold on;
            %xlabel('Threshold', 'FontSize', 16);
            %ylabel('|diff(# Spots)|', 'FontSize', 16);
            
            %Threshold
            %figure(handle1);
            %win_out = smooth(win_out);
            threshold = 0;
            last_min = -1;
            for i = 1:pointcount
                if ~isnan(win_out(i))
                    %if win_out(i) <= stdev_thresh
                    %    threshold = spotcount_table(i,1);
                    %end
                    if (last_min < 0) | (win_out(i) < last_min)
                        threshold = spotcount_table(i,1);
                        last_min = win_out(i);
                    end
                end
                %if threshold > 0
                %    return;
                %end
            end
            
            %line([threshold threshold], get(ax,'YLim'),'Color',ylinecolor,'LineStyle','--');
        end
        
        %%
        %
        function fighandle = drawWindowscorePlot(th_table, win_stdevs, winscore_thresh, intensity_thresh)
            pointcount = size(win_stdevs, 1);
            x = th_table(1:pointcount);
            y = win_stdevs;
            %y = smooth(winstdev);
            xlinecolor = [0.608 0.621 0.628];
            ylinecolor = [0.608 0.621 0.628];
            fighandle = figure(88888);
            ax = axes;
            plot(x, y, 'LineWidth', 2, 'Color', 'k');
            ax.FontSize = 14;
            hold on;
            line(get(ax,'XLim'), [winscore_thresh winscore_thresh],'Color',xlinecolor,'LineStyle','--');
            xlabel('Threshold', 'FontSize', 16);
            ylabel('Window Score', 'FontSize', 16);
            %title(['Window Size = ' num2str(window_size)], 'FontSize', 18);
            title('Auto Threshold Selection', 'FontSize', 18);
            line([intensity_thresh intensity_thresh], get(ax,'YLim'),'Color',ylinecolor,'LineStyle','--');
        end
        
        %%
        %Apply threshold filters to a 2D image (or 3D image slice)
        %ARGS
        %   in_slice (num[Y][X]) - Input image
        %   th1 (num) - Raw threshold value. Pixels below this get filtered first.
        %   th2 (num) - Additional threshold value. Spot nucleating pixels get
        %   filtered if below this.
        %
        %RETURN
        %   img_filtered (num[Y][X]) - Threshold-filtered image.
        %   xx (int[N]) - List of x coordinates of spot-nucleating pixels passed
        %   yy (int[N]) - List of y coordinates of spot-nucleating pixels passed
        %   minVal (num) - Minimum pixel value found in filtered image
        %   maxVal (num) - Maximum pixel value found in filtered image
        %
        %ASSUMPTIONS
        %   -> The image matrix is 2D
        %
        %FLEXIBILITIES
        %   -> "num" can refer to any numerical type
        %   -> "int" can refer to any integer type
        %
        function [img_filtered,xx,yy,minVal,maxVal] = testThreshold_slice(in_slice,th1,th2)
            %Initialize outputs
            xx = NaN(10);
            yy = NaN(10);
            minVal = 0;
            maxVal = 0;

            %Initial filter
            IM = immultiply(in_slice, in_slice > th1); % Filters out pixels below threshold
    
            if sum(IM(:)) > 0 %Are there any spots at all?
                IM2 = imregionalmax(IM,8); %Filters down to only brightest pix in 2D for each spot (so can spot count)
                IM3 = immultiply(in_slice, IM2);
                minVal = min(IM3(:));
                maxVal = max(IM3(:))./2;
        
                [yy,xx] = find(IM3 > round(th2));
                img_filtered = IM3;
            else
                img_filtered = IM; %Nothing to see here
            end

        end

        %%
        %Apply spot detection to 3D image as a whole.
        %ARGS
        %   in_slice (uint16[Y][X][Z]) - Input image
        %   th1 (num) - Raw threshold value. Pixels below this get filtered first.
        %   th2 (num[Z]) - Additional threshold value/plane. Spot nucleating pixels get
        %               filtered if below this.
        %   zBorder (int) - Number of slices (from bottom and top of stack) to
        %                   ignore. These are not cut or blacked out, just ignored.
        %
        %RETURN
        %   img_filtered (uint16[Y][X][Z]) - Threshold-filtered image.
        %   xx (cell{int[N]}) - List of x coordinates of spot-nucleating pixels passed
        %   yy (cell{int[N]}) - List of y coordinates of spot-nucleating pixels passed
        %   minVal (double[Z]) - Minimum pixel value found in filtered image
        %   maxVal (double[Z]) - Maximum pixel value found in filtered image
        %
        %ASSUMPTIONS
        %   -> The image matrix is 3D of 16-bit int values
        %
        %FLEXIBILITIES
        %   -> "int" can refer to any integer type
        %
        function [img_filtered,xx,yy,minVal,maxVal] = testThreshold_3D(in_img,th1,th2,zBorder)
            %Initialize outputs
            %fprintf("Preparing spot testing...\t")
            %tic
            X = size(in_img, 2);
            Y = size(in_img, 1);
            Z = size(in_img, 3);
    
            xx = cell(Z);
            yy = cell(Z);
            minVal = NaN(Z);
            maxVal = NaN(Z);
    
            %Dump unwanted slices
            minZ = zBorder+1;
            maxZ = Z-zBorder;
            %img_filtered = in_img(:,:,minZ:maxZ);
            img_filtered = in_img;
            for z = 1:(minZ-1)
                    img_filtered(:,:,z) = zeros(Y,X);
            end
            for z = (maxZ+1):Z
                    img_filtered(:,:,z) = zeros(Y,X);
            end
            %toc
    
            %Initial filter
            %fprintf("Filter 1\t")
            %tic
            IM = immultiply(img_filtered, img_filtered > th1); % Filters out pixels below threshold
            %toc
    
            if sum(IM(:)) > 0 %Are there any spots at all?
                %fprintf("Filter 2\t")
                %tic
whos
                IM2 = imregionalmax(IM,26); %Filters down to only brightest pix in each spot in 3D (so can spot count)
                %toc
                %size(IM2)
                %size(in_img)
        
                %fprintf("Filter 3\t")
                %tic
                IM3 = immultiply(in_img, IM2);
                %toc
        
                %fprintf("Finding points...\t")
                %tic
                for z = minZ:maxZ
                    %IM3 = immultiply(in_img(:,:,z), IM2);
                    %idx = z - minZ + 1;
                    minVal(z) = min(IM3(:,:,z), [], [1, 2]);
                    maxVal(z) = max(IM3(:,:,z), [], [1, 2])./2;
                    [yz,xz] = find(IM3(:,:,z) > round(th2(z)));
                    xx{z} = xz;
                    yy{z} = yz;
                    img_filtered(:,:,z) = IM3(:,:,z);
                end
                %toc
            else
                img_filtered = IM; %Nothing to see here
            end
        end

        %%
        %Performs spot detection using a 2D projection of the max intensities at
        %each xy pixel position across all planes.
        %ARGS
        %   in_slice (num[Y][X][Z]) - Input image
        %   threshold (num) - Raw threshold value.

        %RETURN
        %   img_filtered (num[Y][X]) - Threshold-filtered image. (Projection)
        %   xx (int[N]) - List of x coordinates of spot-nucleating pixels passed
        %   yy (int[N]) - List of y coordinates of spot-nucleating pixels passed
        %   minVal (num) - Minimum pixel value found in filtered image
        %   maxVal (num) - Maximum pixel value found in filtered image
        %
        %ASSUMPTIONS
        %   -> The image matrix is 3D of 16-bit int values
        %
        %FLEXIBILITIES
        %   -> "num" can refer to any numerical type
        %   -> "int" can refer to any integer type
        %
        function [img_filtered,xx,yy,minVal,maxVal] = testThreshold_maxZ(in_img,threshold)
   
            %Prepare slice
            max_proj = max(in_img,[],3); % uint16[X][Y] - A 2D projection of the brightest pixels in each Z column
            p_std = std2(max_proj);
    
            %Run slice
            [img_filtered,xx,yy,minVal,maxVal] = RNA_Threshold_Common.testThreshold_slice(max_proj, threshold, p_std);
    
        end

        %%
        %Performs spot detection on an image using only the 2D plane with
        %the highest average intensity.
        %ARGS
        %   in_slice (num[Y][X][Z]) - Input image
        %   threshold (num) - Raw threshold value.
        %   zBorder (int) - Number of planes from top and bottom to ignore.   
        %
        %RETURN
        %   img_filtered (num[Y][X]) - Threshold-filtered image. (Single plane)
        %   xx (int[N]) - List of x coordinates of spot-nucleating pixels passed
        %   yy (int[N]) - List of y coordinates of spot-nucleating pixels passed
        %   minVal (num) - Minimum pixel value found in filtered image
        %   maxVal (num) - Maximum pixel value found in filtered image
        %
        %ASSUMPTIONS
        %   -> The image matrix is 3D
        %
        %FLEXIBILITIES
        %   -> "num" can refer to any numerical type
        %   -> "int" can refer to any integer type
        %
        %CAVEATS
        %   -> If there are multiple planes with the same maximum intensity
        %   average,only the first will be used.
        %
        function [img_filtered,xx,yy,minVal,maxVal] = testThreshold_bestAvgZ(in_img,threshold,zBorder)
            %Little more than a wrapper of testThreshold_selectZ that recalculates
            %   the averages every time it is called. 
            %TBH probably shouldn't use.
            Z = size(in_img, 3);
            plane_avgs = NaN(Z);
    
            minZ = zBorder+1;
            maxZ = Z-zBorder;
    
            for z = minZ:maxZ
                plane_avgs(z) = mean2(in_img(:,:,z));
            end
    
            [~,I] = nanmax(plane_avgs(:));
            [img_filtered,xx,yy,minVal,maxVal] = RNA_Threshold_Common.testThreshold_selectZ(in_img,threshold,plane_avgs,I);

        end

        %%
        %Performs spot detection on an image using only a single selected 2D plane.
        %
        %ARGS
        %   in_slice (num[Y][X][Z]) - Input image
        %   threshold (num) - Raw threshold value.
        %   img_avg (double[Z]) - Precalculated average intensities for each plane.
        %   use_planes (int OR int[N]) - Indices of selected plane(s). Only first will be used.
        %
        %RETURN
        %   img_filtered (num[Y][X]) - Threshold-filtered image. (Single plane)
        %   xx (int[N]) - List of x coordinates of spot-nucleating pixels passed
        %   yy (int[N]) - List of y coordinates of spot-nucleating pixels passed
        %   minVal (num) - Minimum pixel value found in filtered image
        %   maxVal (num) - Maximum pixel value found in filtered image
        %
        %ASSUMPTIONS
        %   -> The image matrix is 3D
        %
        %FLEXIBILITIES
        %   -> "num" can refer to any numerical type
        %   -> "int" can refer to any integer type
        %
        %CAVEATS
        %   -> If there are multiple plane indices provided, only the first will be
        %   used.
        %
        function [img_filtered,xx,yy,minVal,maxVal] = testThreshold_selectZ(in_img,threshold,img_avg,use_planes)
    
            use_z = use_planes(1);
            fprintf("Using plane: %d\n", use_z)
            [img_filtered,xx,yy,minVal,maxVal] = RNA_Threshold_Common.testThreshold_slice(in_img(:,:,use_z),threshold,img_avg(use_z));
    
        end
        
         %%
        %Run spot detection on the provided 3D image channel using the
        %provided threshold list and strategy
        %
        %ARGS
        %   img_filter (num[Y][X][Z]) - Pre-filtered image channel
        %   th_list (num[1][T]) - Table containing list of threshold values
        %                         to try
        %   th_strategy (string) - Detection algorithm to use
        %                           'max_proj' - Use max projection
        %                           'max_avg' - Use the plane with max avg
        %                           'all_3d' - 3D analysis (default)
        %   zBorder (int) - Number of planes from each end of stack to ignore
        %   save_filtered (bool) - Whether to save intermediate images
        %   saveStem (string) - Path stem for saving results
        %   collapse3D (bool) - If the all3D strat is set, this specifies
        %           whether to collapse the return coordinates/counts to 2D
        %   verbose (bool) - Whether to print update messages.
        %
        %RETURN
        %   spot_table (double[T][2]) - Table of spot counts at each
        %                               threshold level. First column is
        %                               threshold, second is spot count.
        %   coord_table (cell{int[N][M]}) - Cell array with coordinate
        %                                   table for each threshold level.
        %
        function [spot_table, coord_table] = run_spotDetectOnThresholdList(img_filter, th_list, th_strategy, zBorder,... 
                                                            save_filtered, saveStem, collapse3D, verbose)
    
                                                        
            if nargin < 8
                verbose = true;
            end
            
            %Save the basics
            Z = size(img_filter, 3);
            T = size(th_list, 2);
            minZ = zBorder+1;
            maxZ = Z-zBorder;

            %Pre-allocate some vars we gonna use later...
            plane_avgs = NaN(Z);
            for z = minZ:maxZ
                %plane_avgs(z) = mean2(img_filter(:,:,z));
                plane_avgs(z) = RNA_Threshold_Common.mean_noZeros(img_filter(:,:,z));
            end
    
            spot_table = NaN(T,2); %Columns are threshold, spot count. Rows are entries.
            coord_table = cell(T,1); %Cell vector. Each cell is an x,y table.
    
            %Move this list population up here
            for c = 1:T
                spot_table(c,1) = th_list(1,c);
            end
    
            %Generate coord & spot tables (do spot detection)
            if strcmp(th_strategy, 'max_proj')
                %Use only the max projection slice for spot detection
                for c = 1:T
                    th = th_list(1,c);
                    %Spot detect
                    [f_img,xx,yy,~,~] = RNA_Threshold_Common.testThreshold_maxZ(img_filter,th);
                    %Spot count
                    spot_count = RNA_Threshold_Common.countSpots_total(xx);
                    %spot_table(c,1) = th;
                    spot_table(c,2) = spot_count;
                    %Save coordinates
                    t_coords = [xx(:) yy(:)];
                    coord_table{c,1} = t_coords;
                    if save_filtered
                        savepath = sprintf("%s%s%s", saveStem, '_t', num2str(c));
                        save(savepath, 'f_img');
                    end
                end
            else
                if strcmp(th_strategy, 'max_avg')
                    %Use only the slice w/ highest avg intensity for spot detection
                    [~,I] = nanmax(plane_avgs(:));
                    for c = 1:T
                        th = th_list(1,c);
                        %Spot detect
                        [f_img,xx,yy,~,~] = RNA_Threshold_Common.testThreshold_selectZ(img_filter,th,plane_avgs,I);
                        %Spot count
                        spot_count = RNA_Threshold_Common.countSpots_total(xx);
                        %spot_table(c,1) = th;
                        spot_table(c,2) = spot_count;
                        %Save coordinates
                        t_coords = [xx(:) yy(:)];
                        coord_table{c,1} = t_coords;
                        if save_filtered
                            save([saveStem '_t' num2str(c)], 'f_img');
                        end
                    end
                else
                    %Run full 3D image. Yay!
                    if save_filtered
                        for c = 1:T
                            if verbose; tic; end
                            th = th_list(1,c);
                            if verbose; fprintf("Processing image using threshold = %f\t", th); end
                            %Spot detect
                            [f_img,xx,yy,~,~] = RNA_Threshold_Common.testThreshold_3D(img_filter,th,plane_avgs,zBorder);
                            %Whittle down spots
                            if verbose; toc; end

                            if collapse3D
                                if verbose; tic; end
                                if verbose; fprintf("Determining unique XY spots...\t"); end
                                [spot_count, t_coords] = RNA_Threshold_Common.countSpots_xyUnique(xx,yy);
                                spot_table(c,2) = spot_count;
                                coord_table{c,1} = t_coords;
                                if verbose; toc; end
                            else
                                [spot_count, t_coords] = RNA_Threshold_Common.gen3DCoordTable(xx,yy,minZ,maxZ);
                                spot_table(c,2) = spot_count;
                                coord_table{c,1} = t_coords;
                            end
                            save([saveStem '_t' num2str(c)], 'f_img');
                        end
                    else
                        for c = 1:T
                        %parfor c = 1:T
                            th = th_list(1,c);
                            if verbose
                                tic;
                                fprintf("Processing image using threshold = %f\n", th);
                            end
                            %Spot detect
whos
                            [~,xx,yy,~,~] = RNA_Threshold_Common.testThreshold_3D(img_filter,th,plane_avgs,zBorder);
                            %Whittle down spots
                            if collapse3D
                                if verbose; fprintf("Determining unique XY spots...\n"); end
                                [spot_count, t_coords] = RNA_Threshold_Common.countSpots_xyUnique(xx,yy);
                                spot_table(c,2) = spot_count;
                                coord_table{c,1} = t_coords;
                            else
                                [spot_count, t_coords] = RNA_Threshold_Common.gen3DCoordTable(xx,yy,minZ,maxZ);
                                spot_table(c,2) = spot_count;
                                coord_table{c,1} = t_coords;
                            end
                            if verbose; toc; end
                        end
                    end
                end
            end
                  
        end
        
        %%
        %-------------------------[Spot Functions]-----------------------------
        
        %%
        %Performs a spot count provided a list/matrix/cell vector of coordinates.
        %Given an input with multiple planes, spot count is determined by the 
        %   number of spots in the plane with the most spots
        %
        %ARGS
        %   xx (int[N] OR cell{int[N]}) - List of x (or y) coordinates of spots.
        %
        %RETURN
        %   spot_count (int) - Number of spots counted.
        %
        function spot_count = countSpots_maxSlice(xx)
            
            if iscell(xx)
                Z = size(xx,1);
    
                sizes = zeros(Z,1);
                for z = 1:Z
                    sizes(z) = size(xx{z},1);
                end
    
                spot_count = max(sizes(:,1));
            else
                spot_count = max(size(xx,1));
            end
        end

        %%
        %Performs a spot count provided a list/matrix/cell vector of coordinates.
        %Given an input with multiple planes, spot count is a sum of all spots on
        %all planes.
        %
        %ARGS
        %   xx (int[N] OR cell{int[N]}) - List of x (or y) coordinates of spots.
        %
        %RETURN
        %   spot_count (int) - Number of spots counted.
        %
        function spot_count = countSpots_total(xx)
            
            if iscell(xx)
                Z = size(xx,1);
    
                spot_count = 0;
                for z = 1:Z
                    spot_count = spot_count + size(xx{z},1);
                end
            else
                spot_count = size(xx,1);
            end
        
        end

        %%
        %Performs a spot count provided a list/matrix/cell vector of coordinates.
        %Given an input with multiple planes, spot count is a sum of all spots with
        %a unique x,y coordinate - in other words, spots found in the same position
        %on different planes are only counted once.
        %
        %ARGS
        %   xx (int[N] OR cell{int[N]}) - List of x coordinates of spots.
        %   yy (int[N] OR cell{int[N]}) - List of y coordinates of spots.
        %
        %RETURN
        %   spot_count (int) - Number of spots counted.
        %   unique_spots (double[N][2]) - List of unique spot coordinates. 
        %                                   x in column 1, y in column 2
        %
        function [spot_count, unique_spots] = countSpots_xyUnique(xx,yy)
    
            %Count total points
            Z = size(xx);
    
            total_points = 0;
            for z = 1:Z
                xz = xx{z};
                total_points = total_points + size(xz,1);
            end
    
            %Merge
            full_coord_table = zeros(total_points, 2);
            %fprintf("full_coord_table size: %d", size(full_coord_table));
            main_idx = 1;
            for z = 1:Z
                xz = xx{z};
                yz = yy{z};
                C = size(xz,1);
                if C > 0
                    for c = 1:C
                        %fprintf("main_idx = %d\n", main_idx)
                        full_coord_table(main_idx, 1) = xz(c);
                        full_coord_table(main_idx, 2) = yz(c);
                        main_idx = main_idx + 1;
                    end 
                end  
            end
    
            unique_spots = unique(full_coord_table, 'rows');
            spot_count = size(unique_spots,1);
            %fprintf("Unique points: %d\n", size(unique_coords,1))
        end
        
        %%
        % Convert a stack of find indices to a 3 column xyz coordinate
        % table.
        % This is called in the 3D spot detection process.
        %
        %ARGS
        %   xx (int[Z][N]) - List of x coordinates of found spots for each
        %       z slice
        %   yy (int[Z][N]) - List of y coordinates of found spots for each
        %       z slice
        %   minZ (int) - Index of bottom z slice analyzed
        %   maxZ (int) - Index of top z slice analyzed
        %
        %RETURN
        %   spot_count (int) - Total number of spots counted given provided
        %       parameters.
        %   coord_table (int[N][3]) - Output 3-column coordinate table:
        %       list of xyz coordinates of spots.
        %
        function [spot_count, coord_table] = gen3DCoordTable(xx,yy,minZ,maxZ)
            
            %First, count spots
            spot_count = 0;
        
            for zi = minZ:maxZ
                xlist = xx{zi};
                spot_count = spot_count + size(xlist,1);
            end
            
            %Allocate & copy coordinates
            coord_table = zeros(spot_count, 3);
            cidx = 1;
            for zi = minZ:maxZ
                xlist = xx{zi};
                ylist = yy{zi};
                
                xycount = size(xlist);
                for xyi = 1:xycount
                    %coord_table(cidx) = [xlist(xyi), ylist(xyi), zi];
                    coord_table(cidx,1) = xlist(xyi);
                    coord_table(cidx,2) = ylist(xyi);
                    coord_table(cidx,3) = zi;
                    cidx = cidx+1;
                end
                
            end
            
        end
        
        %%
        % Collapse a list of 3D coordinates to a list of 2D coordinates
        % with spots that only differ in the z direction merged.
        %
        %ARGS
        %   coord_table_3D ([N][3+]) - List of 3D coordinates - first 
        %       column for each row is x coord, second is y, third is z
        %
        %RETURN
        %   coord_table_2D ([M][2]) - List of 2D collapsed coordinates 
        %       where first column is x coord and second is y
        %
        function coord_table_2D = collapse3DCoordTable(coord_table_3D)
            
            xytable = coord_table_3D(:,1:2); %Strip the z column
            coord_table_2D = unique(xytable, 'rows'); %Now can reduce to unique rows!
            
        end
        
        %%-------------------------[Clustering]-----------------------------
        
        %%
        % (Don't use this method - it hasn't been updated to work in 3D and
        % I haven't really tested it)
        %
        % Generate a logical matrix describing which spots are within a
        % certain radius of others.
        %
        %ARGS
        %   coord_table_th ([N][2+]) - Table of coordinates of spots to
        %       analyze
        %   rad (num) - Radius threshold: two spots whose distance is below
        %       this threshold are considered "adjacent"
        %
        %RETURN
        %   adj_matrix (bool[N][N]) - Boolean matrix that is the number of spots tall
        %       and wide, each row and column correlating to a spot in the
        %       order it appears in the input table. "True" values in the
        %       matrix indicates that the two spots are within radius of
        %       each other.
        %
        function adj_matrix = generateAdjacencyMatrix(coord_table_th, rad)
            
            spotcount = size(coord_table_th, 1);
            adj_matrix = zeros(spotcount, spotcount);
            
            %Tests whether spots are within radius of each other
            for i = 1:spotcount
                x0 = coord_table_th(i,1); y0 = coord_table_th(i,2);
                minx = x0 - rad; maxx = x0 + rad;
                miny = y0 - rad; maxy = y0 + rad;
                for j = 1:spotcount
                    x1 = coord_table_th(j,1); y1 = coord_table_th(j,2);
                    if (x1 >= minx && x1 <= maxx)
                        if(y1 >= miny && y1 <= maxy)
                            adj_matrix(i,j) = 1;
                        end
                    end
                end
            end
            
            adj_matrix = logical(adj_matrix);
        end
        
        %%
        %(Don't use this method - it hasn't been updated to work in 3D and
        % I haven't really tested it)
        %   This was a draft method for clustering that may still be
        %   useful, but at the moment I am not using or supporting it.
        %
        %ARGS
        %
        %RETURN
        %
        function [groups, membership_matrix] = clusterPoints(adj_matrix)
            
            %Get j indices for each i...
            spotcount = size(adj_matrix, 1);
            hitlist = cell(spotcount);
            for i = 1:spotcount
                hitlist{i} = find(adj_matrix(i));
            end
           
            %For each i, look through j combos
            %Ignore any i's where the only hit is itself
            %Ignore any j's that are <= i
            %Generate a struct list that is i,j, and the intersect size/union size
            ccount = 0;
            for i = 1:spotcount
                myhits = hitlist{i};
                hitcount = size(myhits,1);
                if hitcount >= 2
                    for k = 1:hitcount
                        j = myhits(k);
                        if j > i
                            ccount = ccount+1;
                        end
                    end
                end
            end
    
            
            idx = 1;
            combos(ccount) = struct('i', -1, 'j', -1, 'jaccard', 0.000, 'isect', cell(1));
            jaccard_mtx = NaN(spotcount, spotcount);
            for i = 1:spotcount
                myhits = hitlist{i};
                hitcount = size(myhits,1);
                if hitcount >= 2
                    for k = 1:hitcount
                        j = myhits(k);
                        if j > i
                            jhits = hitlist{j};
                            ij_union = union(myhits, jhits);
                            ij_isect = intersect(myhits, jhits);
                            jaccard = size(ij_isect,1) / size(ij_union,1);
                            
                            %F@#$ing matlab just put it in a f#$%ing linked list
                            %This is why I hate interpreted languages
                            combos(idx).i = i;
                            combos(idx).j = j;
                            combos(idx).jaccard = jaccard;
                            combos(idx).isect{1} = ij_isect;
                            
                            jaccard_mtx(i,j) = jaccard;
                            
                            idx = idx + 1;
                        end
                    end
                end
            end
            
            %Sort combos by jaccard index, highest at front
            %https://www.mathworks.com/matlabcentral/answers/397385-how-to-sort-a-structure-array-based-on-a-specific-field
            comboTable = struct2table(combos);
            comboTable = sortrows(comboTable, 'jaccard');
            
            %Go through each one to assign groups...
            groups = zeros(spotcount);
            groupid = 1;
            for i = 1:ccount
                cisect_cells = comboTable(i,'isect');
                cisect = cisect_cells{1};
                
                %See what groups the isect spots are already in..
                %Are there any ungrouped isect spots?
                ungrouped = 0;
                bestgroup = 0;
                bestcount = 0;
                icount = size(cisect,1);
                for j = 1:icount
                    g = groups(j);
                    if g == 0
                        ungrouped = ungrouped + 1;
                    else
                        if bestgroup > 0
                            if g == bestgroup
                                bestcount = bestcount+1;
                            else
                                bestcount = 0;
                                bestgroup = -1;
                            end
                        elseif bestgroup == 0
                            bestgroup = g;
                            bestcount = bestcount+1;
                        end
                    end
                end
                
                if ungrouped > 0
                    %Either group them into best group or new group
                    if (bestgroup) > 0 && (bestcount >= icount/4)
                        %bestgroup
                        for j = 1:icount
                            if groups(j) == 0
                                groups(j) = bestgroup;
                            end
                        end
                    else
                        %new group
                        for j = 1:icount
                            if groups(j) == 0
                                groups(j) = groupid;
                            end
                        end
                        groupid = groupid + 1;
                    end
                end
            end
            
            %Rearrange to have list of members by group
            groupcount = groupid-1;
            group_members = cell(groupcount);
            for g = 1:groupcount
                mcount = 0;
                for i = 1:spotcount
                    if groups(i) == g
                        mcount = mcount + 1;
                    end
                end
                
                memlist = zeros(mcount);
                idx = 0;
                for i = 1:spotcount
                    if groups(i) == g
                        memlist(idx) = i;
                        idx = idx+1;
                    end
                end
                
                group_members{g} = memlist;
            end
            
            %Assign membership scores & adjust groups
            membership_matrix = NaN(spotcount, groupcount);
            threshold = 0.5;
            for i = 1:spotcount
                for g = 1:groupcount
                    %Go through every member of the group that isn't the same spot
                    %I'm just gonna average jaccards...
                    gmems = group_members{g};
                    sum = 0.0;
                    usemems = 0;
                    mcount = size(gmems, 1);
                    for m = 1:mcount
                        if gmems(m) ~= i
                            usemems = usemems + 1;
                            sum = sum + jaccard_mtx(i, gmems(m));
                        end
                    end
                    avg = sum/usemems;
                    
                    membership_matrix(i, g) = avg;
                end
                
                %Reassign membership...
                bestgroup = 0;
                bestscore = 0.0;
                for g = 1:groupcount
                    gscore = membership_matrix(i, g);
                    if gscore > bestscore
                        bestgroup = g;
                        bestscore = gscore;
                    end
                end
                
                if bestscore >= threshold
                    %Set group
                    groups(i) = bestgroup;
                end
            end
            
        end
        
        %%
        % Sort detected spots into clusters based on location relative to
        % other spots. Each cluster group has a "nucleating" spot that can
        % be used for downstream processing to eliminate false positive
        % signals or merge clusters into RNA clouds.
        % (this method is unfinished)
        %
        %ARGS
        %   img_filter ([Y][X][Z]) - Pre-processed (filtered) image stack
        %       of channel to analyze. 3D is optional, this method accepts
        %       2D images as well.
        %   coord_table_th (int[N][2+]) - List of spot coordinates (x,y,
        %       and optionally z) for all spots detected at threshold of
        %       interest
        %   r1 (num) - Radius of sub-image to be sample around every
        %       putative spot to determine average spot radius
        %
        %RETURN
        %   merged_coords structs:
        %       same number of spots as input table with new column
        %       Each spot is now a struct with the following fields:
        %           x (int), y (int), z (int), nuc_idx (int), weight(double),
        %           devscore(double), brightscore(double)
        %       If nuc_idx is 0, the spot is unmerged.
        %       If nuc_idx is -1, the spot is a nucleating spot
        %       Non-zero positive values of nuc_idx are the index of the
        %       nucleating spot this spot has been merged to.
        %
        %NOTES
        %   This method will not process tables containing more than
        %   1000000 spots as it's unecessarily memory intensive, it is
        %   likely that spots are so close together they will all cluster
        %   together anyway, and such a threshold that would produce so
        %   many hits is unlikely to be a good threshold to begin with.
        %
        function merged_coords = mergeSpots(img_filter, coord_table_th, r1)
            
            %Flag spots that are "weird"
            %Build an average profile of what a spot looks like in x-y and x-z
            %For each spot coord, take a subimage and scale intensity
            jclassname = 'hospelhornbg_rnafish.SpotMerger';
            spotcount = size(coord_table_th, 1);
            dimcount = size(coord_table_th, 2); %If 2D, add a dummy z coord
            %spotcount
            
            SCOUNT_THRESH = 1000000;
            if spotcount > SCOUNT_THRESH
                %This is to save memory.
                printf("Number of spots exceeds threshold. Cloud detection will not proceed.\n");
                merged_coords(1) = struct('x', -1, 'y', -1, 'z',-1, 'nuc_idx', -1, 'weight', 0.0, 'devscore', 0.0, 'brightscore', 0.0);
                return;
            end
            
            max_proj = max(img_filter,[],3);
            max_proj_d = double(max_proj);
            LminF = median(max_proj_d(:)) - round(0 * std(max_proj_d(:)));
            LmaxF = median(max_proj_d(:)) + round(10 * std(max_proj_d(:)));
            %figure(123);
            %imshow(max_proj_d, [LminF, LmaxF]);
            
            %Extract spot images and take average/stdev
            %Initialize Java object
            dia = 2*r1 + 1;
            if dimcount == 2
                s_imgs = NaN(dia, dia, spotcount);
                for s = 1:spotcount
                    s_imgs(:,:,s) = RNA_Threshold_Common.get_image_around(img_filter, coord_table_th(s,1), coord_table_th(s,2), 0, r1, 1);
                end
                avg_spot = nanmean(s_imgs,3);
                std_spot = nanstd(s_imgs, 0, 3);
                jmerger = javaMethod('importFromU16_2D_noimg', jclassname, coord_table_th, r1);
                javaMethod('setAverageSpot_2D', jmerger, avg_spot);
            elseif dimcount == 3
                s_imgs = NaN(dia, dia, dia, spotcount);
                for s = 1:spotcount
                    s_imgs(:,:,:,s) = RNA_Threshold_Common.get_image_around(img_filter, coord_table_th(s,1), coord_table_th(s,2), coord_table_th(s,3), r1, 1);
                end
                avg_spot = nanmean(s_imgs,4);
                std_spot = nanstd(s_imgs, 0, 4);
                jmerger = javaMethod('importFromU16_3D_noimg', jclassname, coord_table_th, r1);
                javaMethod('setAverageSpot_3D', jmerger, avg_spot);
            end
            
            %Initialize java object
            %avg_spot_xy = uint16(javaMethod('getAverageSpot_double', jmerger));
            %figure(121);
            %imshow(avg_spot_xy, [LminF, LmaxF]);
            
            %Determine adj matrix radius from average image...
            maxdbl = max(avg_spot, [], 'all');
            avg_spot16 = uint16(avg_spot);
            max16 = uint16(round(maxdbl));
            binsz = 5;
            binct = (65535/binsz) + 1;
            [histo_counts] = imhist(avg_spot16, binct);
            sm = smooth(histo_counts);
            sm = smooth(sm);
            %figure(123);
            %plot(sm);
            deriv = diff(sm);
            %find first max
            max1 = -1;
            for i = 2:max16
                if deriv(i) < 0 && deriv(i-1) > 0
                    max1 = i;
                    break;
                end
            end
            
            %find first dip below (a threshold) after max
            binthresh = 25;
            autot = 0;
            for i = max1:max16
                if sm(i) <= binthresh
                    autot = i * binsz;
                    break;
                end
            end
            %return;
            autot
            
            %autot = Background_Extractor.autoestStdevThreshold_lowweight(avg_spot_xy, max16);
            javaMethod('setBackgroundThreshold_AvgSpot', jmerger, autot);
            avgrad = javaMethod('calculateAverageSpotRadius', jmerger);
            
            %Generate adj matrix
            javaMethod('generateDistMatrix', jmerger);
           
            %Move these calculations to matlab then pass results back to java
            %javaMethod('scoreDeviation', jmerger); 
            %javaMethod('scoreBrightness', jmerger);
            brightScores = NaN(spotcount,1);
            devScores = NaN(spotcount,1);
            if dimcount == 2
                for s = 1:spotcount
                    simg = double(s_imgs(:,:,s));
                    brightScores(s,1) = nanmean(simg, 'all');
                    devScores(s,1) = nanmean(abs(avg_spot - simg)./std_spot);
                end
            elseif dimcount == 3
                for s = 1:spotcount
                    simg = double(s_imgs(:,:,:,s));
                    brightScores(s,1) = nanmean(simg, 'all');
                    devScores(s,1) = nanmean(abs(avg_spot - simg)./std_spot, 'all');
                end
            end
            javaMethod('setBrightnessScores', jmerger, brightScores(:,1)); 
            javaMethod('setDevScores', jmerger, devScores(:,1));
            
            %Cluster
            if dimcount == 3
                %Square the initial merge radius
                avgrad = avgrad * avgrad;
                fprintf("3D input - cluster radius increased to %f\n", avgrad);
                javaMethod('setMergeRadius', jmerger, avgrad);
            end
            javaMethod('clusterSpots', jmerger);
            javaMethod('mergeGroups', jmerger);
            
            spotgroups = javaMethod('getGroupAssignments', jmerger);
            javaMethod('getMembershipMatrix', jmerger);
            %spotgroups
            %memmtx

            javaMethod('determineNucleatingSpots', jmerger);
            nucleating = javaMethod('getNucleatingList', jmerger);
            %nucleating
            
            %DEBUG - lop a chunk out of the image to examine more closely
            %Also build dist table between all points in this subimage
            
            %RNA_Threshold_Common.saveClusteringSubImage(double(img_filter), LminF, LmaxF, coord_table_th, spotgroups, 908, 1142, 1240, 1474);
            %return;
            
            figure(124);
            imshow(max_proj_d, [LminF, LmaxF]);
            for i = 1:spotcount
                javai = i-1;
                hold on;
                x = coord_table_th(i, 1);
                y = coord_table_th(i, 2);
                g = spotgroups(i);
                if g == 0
                    plot(x, y, 'or', 'markersize', 10); 
                else
                    if (nucleating(g+1) == javai)
                        text(x, y, num2str(g), 'Color', 'y');
                    else
                        text(x, y, num2str(g), 'Color', 'g');
                    end
                end
            end
            %'2D view complete'
            
            if dimcount == 3
                figure(33);
                %hold on;
                
                %Set axes...
                X = size(img_filter,2);
                Y = size(img_filter,1);
                Z = size(img_filter,3);
            
                %ax = gca;
                ax = axes();
                ax.XLim = [1, X];
                ax.YLim = [1, Y];
                ax.ZLim = [1, Z];
                ax.XAxisLocation = 'top';
                ax.YDir = 'reverse';
                ax.XLabel.String = 'X Axis';
                ax.YLabel.String = 'Y Axis';
                ax.ZLabel.String = 'Z Axis';
                
                hold on;
                set(gca, 'Color', 'k');
                for i = 1:spotcount
                    javai = i-1;
                    x = coord_table_th(i, 1);
                    y = coord_table_th(i, 2);
                    z = coord_table_th(i, 3);
                    g = spotgroups(i);
                    if g == 0
                        %plot3(x, y, z, 'or', 'markersize', 10); 
                        text(x, y, z, num2str(i), 'Color', 'c');
                    else
                        if (nucleating(g+1) == javai)
                            %text(x, y, z, num2str(g), 'Color', 'y');
                            text(x, y, z, num2str(i), 'Color', 'y');
                        else
                            %text(x, y, z, num2str(g), 'Color', 'g');
                            text(x, y, z, num2str(i), 'Color', 'g');
                        end
                    end
                end
                
                RNA_Threshold_Common.render3DImage(img_filter, autot);
                
            end
            
            return;
            
            %Save and return
            %devscores = javamethod('getDevScores', jmerger);
            %brightscores = javamethod('getBrightnessScores', jmerger);
            merged_coords(spotcount) = struct('x', 0, 'y', 0, 'z', 0, 'nuc_idx', 0, 'weight', 0.0, 'devscore', 0.0, 'brightscore', 0.0);
            for s = 1:spotcount
                merged_coords(s).x = coord_table_th(s,1);
                merged_coords(s).y = coord_table_th(s,2);
                
                if dimcount == 2
                    merged_coords(s).z = 1;
                elseif dimcount == 3
                    merged_coords(s).z = coord_table_th(s,3);
                end
                
                merged_coords(s).weight = 1.0;
                merged_coords(s).devscore = devscores(s-1); %Java index
                merged_coords(s).brightscore = brightscores(s-1); %Java index
                
                merged_coords(s).nuc_idx = 0;
                g = spotgroups(i);
                if g ~= 0
                    nspot = nucleating(g+1); %Java returns array w/ 0 slot
                    if nspot == s
                        merged_coords(s).nuc_idx = -1;
                    else
                        merged_coords(s).nuc_idx = nspot;
                    end
                end
            end
            
        end
        
        %%
        %  [Debug method, use if you need. You will have to change the path though.]
        %   Renders sub-images (three sets of z slices) marking spots,
        %   first cicled, than by index in table, then by group. These
        %   images are saved to the output path (hard-coded in script,
        %   change if needed).
        %   Also calculates distance between spots and outputs distance
        %   tables as csv files.
        %
        %ARGS
        %   img (num[Y][X][Z]) - Original image channel to use for output
        %       image displays
        %   LminF (num) - Min scale value to pass to imshow
        %   LmaxF (num) - Max scale value to pass to imshow
        %   coord_table (num[N][2+]) - List of coordinates of detected
        %       spots (columns being x,y and optionally z). Any columns
        %       beyond the third will be ignored.
        %   groups (int[N]) - Array the size of the number of spots
        %       where the element at the index correlating to a spot in the
        %       coord table is the group the spot is part of. 0 means it is
        %       ungrouped.
        %   x1 (num) - Left bound x coordinate of area in image to sample
        %   x2 (num) - Right bound x coordinate of area in image to sample
        %   y1 (num) - Top bound y coordinate of area in image to sample
        %   y2 (num) - Bottom bound y coordinate of area in image to sample
        %
        %
        function saveClusteringSubImage(img, LminF, LmaxF, coord_table, groups, x1, x2, y1, y2)
            %Since this is a debug function, I'll leave some stuff kind of
            %inefficient
            
            outstem = 'D:\usr\bghos\labdat\imgproc\cluster_debug';
            %spotcount = size(coord_table, 1);
            
            %Get subset of spots and remap the coords
            inbox = coord_table(:,1) > x1 & coord_table(:,1) < x2 & coord_table(:,2) > y1 & coord_table(:,2) < y2; 
            [boxed_idxs, ~] = find(inbox);
            boxed_count = size(boxed_idxs,1);
            boxed_coords = zeros(boxed_count,3);
            boxed_coords(:,1) = coord_table(boxed_idxs,1) - x1;
            boxed_coords(:,2) = coord_table(boxed_idxs,2) - y1;
            
            coord_info_tbl = zeros(boxed_count, 4);
            coord_info_tbl(:,1) = boxed_idxs(:,1);
            coord_info_tbl(:,2:4) = coord_table(boxed_idxs,1:3);
            csvwrite([outstem '\spot_coords.csv'], coord_info_tbl);
            
            %Dist tables (X,Y) (X,Y,Z) (X,Y, scaled Z)
            dist_tbl = NaN(boxed_count, boxed_count);
            for i = 2:boxed_count
                ii = boxed_idxs(i);
                for j = 1:i
                    jj = boxed_idxs(j);
                    xdist = (coord_table(ii,1) - coord_table(jj,1))^2;
                    ydist = (coord_table(ii,2) - coord_table(jj,2))^2;
                    dist_tbl(i,j) = sqrt(xdist + ydist);
                end
            end
            csvwrite([outstem '\disttbl_xy.csv'], dist_tbl);
            
            for i = 2:boxed_count
                ii = boxed_idxs(i);
                for j = 1:i
                    jj = boxed_idxs(j);
                    xdist = (coord_table(ii,1) - coord_table(jj,1))^2;
                    ydist = (coord_table(ii,2) - coord_table(jj,2))^2;
                    zdist = (coord_table(ii,3) - coord_table(jj,3))^2;
                    dist_tbl(i,j) = sqrt(xdist + ydist + zdist);
                end
            end
            csvwrite([outstem '\disttbl_xyz.csv'], dist_tbl);
            
            for i = 2:boxed_count
                ii = boxed_idxs(i);
                for j = 1:i
                    jj = boxed_idxs(j);
                    xdist = (coord_table(ii,1) - coord_table(jj,1))^2;
                    ydist = (coord_table(ii,2) - coord_table(jj,2))^2;
                    zdist = 4 * (coord_table(ii,3) - coord_table(jj,3))^2;
                    dist_tbl(i,j) = sqrt(xdist + ydist + zdist);
                end
            end
            csvwrite([outstem '\disttbl_xy4z.csv'], dist_tbl);
            
            %Images: normal (w/circles), group marked, number marked
            subimg = img(y1:y2, x1:x2, :);
            Z = size(img,3);
            for z = 1:Z
                
                %Get subset of spots on this plane
                inplane = inbox & (coord_table(:,3) == z);
                if isempty(inplane)
                    f_handle = figure(7);
                    imshow(subimg(:,:,z), [LminF, LmaxF]);
                    hold on;
                    saveas(f_handle, [outstem '\normal\subimg_z' num2str(z) '.png']);
                    saveas(f_handle, [outstem '\group\subimg_z' num2str(z) '.png']);
                    saveas(f_handle, [outstem '\numbers\subimg_z' num2str(z) '.png']);
                    close(f_handle);
                else 
                    [plane_idxs, ~] = find(inplane);
                    spotcount = size(plane_idxs, 1);
                    
                    f_handle = figure(7);
                    imshow(subimg(:,:,z), [LminF, LmaxF]);
                    hold on;
                    for i = 1:spotcount
                        ii = plane_idxs(i,1);
                        hold on;
                        x = coord_table(ii, 1) - x1;
                        y = coord_table(ii, 2) - y1;
                        plot(x, y, 'or', 'markersize', 10); 
                    end 
                    saveas(f_handle, [outstem '\normal\subimg_z' num2str(z) '.png']);
                    close(f_handle);
                
                    f_handle = figure(7);
                    imshow(subimg(:,:,z), [LminF, LmaxF]);
                    hold on;
                    for i = 1:spotcount
                        ii = plane_idxs(i,1);
                        hold on;
                        x = coord_table(ii, 1) - x1;
                        y = coord_table(ii, 2) - y1;
                        g = groups(ii);
                        if g == 0
                            plot(x, y, 'or', 'markersize', 10); 
                        else
                            text(x, y, num2str(g), 'Color', 'g');
                        end
                    end
                    saveas(f_handle, [outstem '\group\subimg_z' num2str(z) '.png']);
                    close(f_handle);
                    
                    f_handle = figure(7);
                    imshow(subimg(:,:,z), [LminF, LmaxF]);
                    hold on;
                    for i = 1:spotcount
                        ii = plane_idxs(i,1);
                        hold on;
                        x = coord_table(ii, 1) - x1;
                        y = coord_table(ii, 2) - y1;
                        text(x, y, num2str(ii), 'Color', 'c');
                    end
                    saveas(f_handle, [outstem '\numbers\subimg_z' num2str(z) '.png']);
                    close(f_handle);
                    
                end

            end
            
        end
        
        %%-------------------------[Mask]-----------------------------
        
        %%
        %Scale spot counts generated from analysis of a masked image to
        %account for relative coverage of mask.
        %
        %ARGS
        %   mask_raw (bool[Y][X]) - Raw image mask
        %   spot_table_raw (double[T][2]) - Raw spot count table
        %
        %RETURN
        %   scaled_spot_table (double[T][2]) - Table with scaled counts
        %
        function scaled_spot_table = scaleMaskSpotTable(mask_raw, spot_table_raw)
            
            total_pix = size(mask_raw,1).*size(mask_raw,2);
            nz = find(mask_raw);
            unmasked_pix = size(nz,1);
            
            scale_factor = total_pix./unmasked_pix;
            t_count = size(spot_table_raw, 1);
            
            fprintf("Total Pixels: %d\n", total_pix)
            fprintf("Unmasked Pixels: %d\n", unmasked_pix)
            fprintf("Scale Factor: %f\n", scale_factor)
            
            scaled_spot_table = NaN(t_count, 2);
            scaled_spot_table(1:t_count,1) = spot_table_raw(1:t_count,1);
            
            for i = 1:t_count
                scaled_spot_table(i,2) = spot_table_raw(i,2) .* scale_factor;
            end
            
        end
        
        %%
        %Scale spot counts generated from analysis of a masked image to
        %account for relative coverage of mask. Save new spot table to
        %analysis directory in place of old table (old table is moved to a
        %different file).
        %
        %ARGS
        %   mask_path (string) - Path to mat file containing raw mask
        %   mask_analysis_path_stem (string) - Path stem of previous analysis
        %
        function scale_mask_spot_table(mask_path, mask_analysis_path_stem)
            
            spot_table_path = [mask_analysis_path_stem '_spotTable'];
            load(spot_table_path, 'spot_table');
            raw_table = spot_table;
            
            load(mask_path, 'background_mask');
            scaled_spot_table = RNA_Threshold_Common.scaleMaskSpotTable(background_mask, raw_table);
            
            spot_table = raw_table;
            save([mask_analysis_path_stem '_spotTableRaw'], 'spot_table');
            
            spot_table = scaled_spot_table;
            save([mask_analysis_path_stem '_spotTable'], 'spot_table');
        end
        
        %%
        % Generate a new spot count table and coordinate table that masks
        % out any spots that are outside the provided mask. This is meant
        % to be used with a background mask to filter for spots that fall
        % in the background area of the image.
        %
        %ARGS
        %   mask_raw (logical[Y][X]([Z])) - Raw 2D or 3D background mask
        %   coord_table (cell{int[N][2+]}) - Coordinate table of spots detected
        %       from original image. Each cell correlates to a tried
        %       threshold.
        %   th_list (int[1][T]) - [Optional] Table of threshold values.
        %       Defaults to 1:1:(coord_table size)
        %
        %RETURN
        %   masked_spot_table (int[T][2]) - Table of thresholds and the
        %       number of spots detected at that threshold in the masked
        %       region only.
        %   masked_coord_table (cell{int[N][INDIM]}) - Coordinate table of
        %       only detected spots in the masked region.
        %
        function [masked_spot_table, masked_coord_table] = mask_spots(mask_raw, coord_table, th_list)
            tcount = size(coord_table,1);
            
            if nargin < 3
                th_list = [1:1:tcount];
            end
            
            masked_spot_table = NaN(tcount,2); %Columns are threshold, spot count. Rows are entries.
            masked_coord_table = cell(tcount,1); %Cell vector. Each cell is an x,y,z table.
            
            for c = 1:tcount
                masked_spot_table(c,1) = th_list(1,c);
            end
            
            mdim = ndims(mask_raw);
            ndim = size(coord_table{1},2);
            tdims = max(mdim, ndim);
            
   
            for t = 1:tcount
                ttable = coord_table{t};
                s_count = size(ttable,1);
                tmptbl = NaN(s_count,tdims);
                for i = 1:s_count
                    x = ttable(i,1);
                    y = ttable(i,2);
                    z = 1;
                    if ndim == 3
                        z = ttable(i,3);
                    end
                    
                    if mdim == 2
                        if mask_raw(y,x) 
                            tmptbl(i,1) = x;
                            tmptbl(i,2) = y;
                            if ndim == 3
                                tmptbl(i,3) = z;
                            end
                        end
                    elseif mdim == 3
                        if mask_raw(y,x,z)
                        	tmptbl(i,1) = x;
                        	tmptbl(i,2) = y;
                         	tmptbl(i,3) = z;
                        end
                    end
                end
                nanrows = any(isnan(tmptbl),2);
                ttable2 = tmptbl(~nanrows, :);
                n_count = size(ttable2,1);
                masked_coord_table{t} = ttable2;
                masked_spot_table(t,2) = n_count;
            end
        end
        
        %%-------------------------[Visualization]-----------------------------
        
        %%
        %   (Debug function)
        %   Attempt to render an approximation of the 3D image as a MATLAB
        %   3D plot. Pixels above the background threshold will be mapped
        %   as points and colored according to how they are binned by
        %   intensity.
        %   The result has poor contrast, but it seems to do the trick if
        %   you know what you're looking for.
        %
        %ARGS
        %   img (num[Y][X][Z]) - 3D image to try to render in plot.
        %   bkg_thresh (num) - Intensity value below which a pixel should
        %       be considered background and not plotted.
        %
        %RETURN
        %   fhandle (figure handle) - Handle of resulting MATLAB figure.
        %
        function fhandle = render3DImage(img, bkg_thresh)
            'Now rendering 3D image (debug)'
            fhandle = figure(333);
            
            maxval = max(img, [], 'all');

            %Set axes...
            X = size(img,2);
            Y = size(img,1);
            Z = size(img,3);
            
            %ax = gca;
            ax = axes();
            ax.XLim = [1, X];
            ax.YLim = [1, Y];
            ax.ZLim = [1, Z];
            ax.XAxisLocation = 'top';
            ax.YDir = 'reverse';
            ax.XLabel.String = 'X Axis';
            ax.YLabel.String = 'Y Axis';
            ax.ZLabel.String = 'Z Axis';
                
            hold on;
            set(gca, 'Color', 'k');
            %Bin values above threshold in intervals of 1/10 maxval - bkg
            bincount = 10;
            binsize = (maxval - bkg_thresh)/bincount;
            binbot = bkg_thresh;
            bintop = binbot + binsize;
            for bin = 1:bincount
                shade = bin/bincount;
                %https://www.mathworks.com/matlabcentral/answers/789-using-find-in-a-3d-matrix-in-matlab
                if bin == bincount
                    bintop = maxval+1;
                end
                [y,x,z] = ind2sub(size(img),find(img >= binbot & img < bintop));
                scatter3(x, y, z, 8, [shade, shade, shade], '.'); 
                
                binbot = bintop;
                bintop = binbot + binsize;
            end
        end
        
        %%-------------------------[Misc. Utilities]-----------------------------
        
        %%
        %Get the standard deviation of all non-zero values in a matrix
        %ARGS
        %   my_matrix (num[][][...]) - A numerical matrix of any dimensions
        %
        %RETURN
        %   stdev (double) - Standard deviation of all non-zero values
        %
        function stdev = std_noZeros(my_matrix)
            [~,~,nz_vals] = find(my_matrix);
            stdev = std(nz_vals);
        end
        
        %%
        %Get the mean of all non-zero values in a matrix
        %ARGS
        %   my_matrix (num[][][...]) - A numerical matrix of any dimensions
        %
        %RETURN
        %   average (double) - Mean of all non-zero values
        %
        function average = mean_noZeros(my_matrix)
            [~,~,nz_vals] = find(my_matrix);
            average = mean(nz_vals);
        end
        
        %%
        %   Get a rectangular (prism) sub-image of size 2r+1 x 2r+1 around the
        %   point specified by x,y, and z. If the image ends before r in
        %   any direction from x,y,z can be reached, those pixels are set
        %   to NaN.
        %
        %ARGS
        %   img (num[Y][X][Z]) - Image to sample from
        %   x (int) - x coordinate of spot to sample around
        %   y (int) - y coordinate of spot to sample around
        %   z (int) - z coordinate of spot to sample around
        %   r (int) - "Radius" or number of pixels in each direction from
        %       spot to sample
        %   z_scale (double) - Ratio of distance between z pixels to xy
        %       pixels.
        %
        %RETURN
        %   sub_img (num[2r+1][2r+1][2r+1]) - Requested sub-image.
        %
        function sub_img = get_image_around(img, x, y, z, r, z_scale)
            %Generates a sub-image centered around the specified coordinates
            z_rad = r.*z_scale;
            X = 2*r + 1;
            Y = X;
            Z = 2*z_rad + 1;
            dims = ndims(img);
            
            if(dims == 3)
                sub_img = NaN(Y,X,Z);
            else
                sub_img = NaN(Y,X);
            end
            
            %Determine ranges to copy
            %---- x
            xs0 = x-r; xs1 = x+r;
            xt0 = 1; xt1 = X;
            dlo = 1 - xs0;
            dhi = xs1 - size(img, 2);
            if dlo > 0
                xs0 = xs0 + dlo;
                xt0 = xt0 + dlo;
            end
            if dhi > 0
                xs1 = xs1 - dhi;
                xt1 = xt1 - dhi;
            end
            
            %---- y
            ys0 = y-r; ys1 = y+r;
            yt0 = 1; yt1 = Y;
            dlo = 1 - ys0;
            dhi = ys1 - size(img, 1);
            if dlo > 0
                ys0 = ys0 + dlo;
                yt0 = yt0 + dlo;
            end
            if dhi > 0
                ys1 = ys1 - dhi;
                yt1 = yt1 - dhi;
            end
            
            %---- z
            if dims == 3
                zs0 = z-z_rad; zs1 = z+z_rad;
                zt0 = 1; zt1 = Z;
                dlo = 1 - zs0;
                dhi = zs1 - size(img, 3);
                if dlo > 0
                    zs0 = zs0 + dlo;
                    zt0 = zt0 + dlo;
                end
                if dhi > 0
                    zs1 = zs1 - dhi;
                    zt1 = zt1 - dhi;
                end
            end
            
            %Do the copy
            if(dims == 3)
                sub_img(yt0:yt1,xt0:xt1,zt0:zt1) = img(ys0:ys1,xs0:xs1,zs0:zs1);
            else
                sub_img(yt0:yt1,xt0:xt1) = img(ys0:ys1,xs0:xs1);
            end
            
            
        end
              
        %%
        %   (In progress)
        %   Check whether there are any background pixels along the line
        %   between two points as specified by the background threshold.
        %
        %ARGS
        %   img (num[Y][X][Z]) - Image coordinates refer to.
        %   xyz1 (int[3]) - Coordinates (x,y,z) of first point
        %   xyz2 (int[3]) - Coordinates (x,y,z) of second point
        %   bkg_thresh (num) - Intensity threshold below which pixels in
        %       image are considered "background".
        %
        %RETURN
        %   hasbkg (bool) - True if there are background pixels between the
        %       two spots, false if not.
        %
        function hasbkg = detectBkgBetweenPoints(img, xyz1, xyz2, bkg_thresh)
            %TODO
        end
        
        %%
        %   (In progress)
        %   Generate a matrix that indicates whether two given spots with
        %   the coordinates provided in the table are adjoined, that is,
        %   have no background pixels between them as determined by the
        %   provided background threshold.
        %
        %ARGS
        %   img (num[Y][X][Z]) - Image coordinates refer to.
        %   coords (int[N][2+]) - List of spot coordinates. 2D xy tables
        %       are acceptable.
        %   bkg_thresh (num) - Intensity threshold below which pixels in
        %       image are considered "background".
        %
        %RETURN
        %   adjmtx (bool[N][N]) - Matrix indicating which spots are
        %       adjoined - "true" means the spots with indices correlating
        %       to the row and column are adjoined (there are no background
        %       pixels between them), "false" means there is background
        %       between the two spots.
        %
        function adjmtx = checkSpotAdj(img, coords, bkg_thresh)
               %Checks to see for each combination of two spots
               % if those two spots have any background color pixels
               % between them when a line is drawn from one to the other
               
               %A TRUE in the output matrix says no bkg between these spots
               spotcount = size(coords, 1);
               adjmtx = false(spotcount, spotcount);
               
               
        end
        
    end
end