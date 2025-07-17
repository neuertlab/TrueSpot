%%
%Image background extraction methods & GUI
%Blythe Hospelhorn
%Version 2.0.1
%Modified March 12, 2024

%Update Log:
%   1.0.0 | 21.03.11
%       Initital proper documentation
%   2.0.0 | 21.12.13
%       Changed initialize() to take the channel directly instead -
%       converted old initialize() to read_and_initialize()
%   2.0.1 | 24.03.12
%       Now calls CellSeg.openCellMask to open cell mask file. That way,
%       can take more formats.
%

%%
classdef Background_Extractor
    
    %GUI Controls:
            %   space (or any key) - Ready listener
            %   
            %   Left Click - Set filter minimum
            %   Middle Click - Exit
            %   Right Click - Set filter maximum
            
            %   x - Exit
            %   s - Save mask
            %   h - Close holes (imopen)
            %   e - Erode holes (imerode)
            %   r - Reset image
            
            %   [] - Decrease/Increase open radius
            %   {} - Decrease/Increase erode radius
            %   +- - Increase/Decrease fine tune amount
            %   <> - Decrease/Increase last moved threshold by fine tune
                        
    properties
        img_ch; %The 3D image stack of light channel (raw)
        cell_mask; %Mask for cutting out pre-detected cells
        max_proj; %Max projection (for viewing)
        
        stdev_proj; %Stdev projection
        stdev_max;
        bkg_mask; %Final background mask
        
        stdev_thresh_min; %Minimum stdev (pixels below are masked out)
        stdev_thresh_max; %Maximum stdev (pixels above are masked out)
        
        exit_ok; %Looper for GUI
        inc_value; %Current increment value for GUI
        last_thresh; %Last threshold modified (for GUI)
        img_handles; %Handles for figures so they can be closed or modded easily
        
        save_path; %Save path to use for bkg mask if called thru GUI
        
        open_rad; %Current radius for imopen function
        erode_rad; %Current radius for imerode function
    end
    
    methods
        
        %%
        %Load Background_Extractor object with tiff image,
        %providing the channel number and index of the light channel.
        %
        %ARGS
        %   tiffpath (string) - Path to source TIFF
        %   total_ch (int) - Total number of channels in source image
        %   light_ch_idx (int) - Index of light channel in source image
        %   seg_file_path (string) - Path to mat file containing cell
        %                            segmentation data for image
        %   verbose (bool) - Whether to print TIFF read warnings
        %
        function obj = read_and_initialize(obj, tiffpath, total_ch, light_ch_idx, seg_file_path, verbose)
            
            %Read in stack
            [tif, ~] = LoadTif(tiffpath, total_ch, [light_ch_idx], ~verbose);
            obj = obj.initialize(tif{light_ch_idx,1}, seg_file_path);
            
        end
        
        %%
        %Load Background_Extractor object with provided channel data.
        %
        %ARGS
        %   light_ch (double[][][]) - 3D light channel data
        %   seg_file_path (string) - Path to mat file containing cell
        %                            segmentation data for image
        %
        function obj = initialize(obj, light_ch, seg_file_path)
            cells = CellSeg.openCellMask(seg_file_path); %More formats.
            obj.cell_mask = (cells == 0);
            obj = obj.initializeMaskLoaded(light_ch);
        end

        function obj = initializeMaskLoaded(obj, light_ch)
            obj.img_ch = light_ch;

            %Generate stdev projection
            obj.stdev_proj = std(obj.img_ch, 0, 3);
            
            %Generate max projection (for viewing)
            obj.max_proj = max(obj.img_ch,[],3);
            
            %Set defaults
            obj.stdev_thresh_min = 0;
            obj.stdev_thresh_max = 65535;
            obj.save_path = 'background_mask';
            obj.open_rad = 8;
            obj.erode_rad = 8;
            
            %Debug
            %figure(2348);
            %imshow(obj.max_proj,[]);
            
            %figure(256894);
            %imshow(obj.stdev_proj,[]);
            
            %Autodetect threshold
            stdev16 = uint16(obj.stdev_proj);
            max_stdev = max(max(stdev16(:,:)));
            autot = Background_Extractor.autoestStdevThreshold_bimodal(stdev16, max_stdev);
            
            if autot < 0
                %No bimodal valley found. Try enhancing contrast
                obj = obj.enhanceContrast();
                stdev16 = uint16(obj.stdev_proj);
                max_stdev = max(max(stdev16(:,:)));
                autot = Background_Extractor.autoestStdevThreshold_lowweight(stdev16, max_stdev);
                obj.stdev_thresh_max = autot;
                %obj.open_rad = 12;
                %obj.erode_rad = 12;
            else
                obj.stdev_thresh_max = autot;
            end
            
            obj = obj.updateBkgMask();
        end
        
        %%
        %Launch the background threshold selection GUI.
        %
        function obj = launchThresholdSelector(obj)
            
            %obj.bkg_mask = obj.cell_mask;
            
            %Set defaults
            %obj.stdev_thresh_min = 0;
            %obj.stdev_thresh_max = 65535;
            
            obj.inc_value = 10; %Fine tuning increment defo
            obj.last_thresh = 0; %0 if min, 1 if max
            obj.exit_ok = 0;
            
            obj = obj.renderThSelectFigures(); 
            
            while obj.exit_ok == 0
                w = waitforbuttonpress;
                if w
                    obj = obj.onReadyKey();
                end
            end
            
            obj = obj.updateBkgMask();
        end
        
        %%
        %Render the threshold selection figures.
        %
        function obj = renderThSelectFigures(obj)
            
            obj.img_handles = cell(3);
            color3 = [0.608 0.621 0.628]; %#A8ABAD
            dim = [.45, .4, .5, .5];
            
            %First figure
            %(Left) Max proj, (Right) stdev proj - both unmasked
            obj.img_handles{1} = figure(333);
            clf;
            
            subplot(1,2,1);
            imshow(obj.max_proj,[]);
            hold on;
            title("Max Projection");
            impixelinfo;
            
            subplot(1,2,2);
            imshow(obj.stdev_proj,[]);
            hold on;
            title("Stack STDEV Projection");
            impixelinfo;
            
            %Second figure, max proj masked
            %obj = obj.updateBkgMask();
            masked_proj = immultiply(obj.max_proj, obj.bkg_mask);
            obj.img_handles{2} = figure(334);
            clf;
            imshow(masked_proj,[]);
            hold on;
            title("Max Projection (Masked)");
            impixelinfo;
            
            %Third figure, histogram of stdev proj
            stdev16 = uint16(obj.stdev_proj);
            max_stdev = max(max(stdev16(:,:)));
            obj.stdev_max = max_stdev;
            if obj.stdev_thresh_max > max_stdev
                obj.stdev_thresh_max = max_stdev;
            end
            th_min = obj.stdev_thresh_min;
            th_max = obj.stdev_thresh_max;
            obj.img_handles{3} = figure(335);
            clf;
            ax = axes;
            [counts,~] = imhist(stdev16, max_stdev);
            imhist(stdev16, max_stdev);
            hold on;
            max_bin_sz = max(counts);
            axis([0 max_stdev 0 max_bin_sz]);
            title("STDEV Histogram");
            line([th_min th_min], get(ax,'YLim'),'Color',color3,'LineStyle','--');
            line([th_max th_max], get(ax,'YLim'),'Color',color3,'LineStyle','--');
            %line([th_min th_min], [0, max_bin_sz],'Color',color3,'LineStyle','--');
            %line([th_max th_max], [0, max_bin_sz],'Color',color3,'LineStyle','--');
            t1_str = ['Minimum: ' num2str(th_min)];
            t2_str = ['Maximum: ' num2str(th_max)];
            anno_str = {t1_str, t2_str};
            annotation('textbox',dim,'String',anno_str,'FitBoxToText','on');
            %set(obj.img_handles{3}, 'KeyPressFcn', @(obj)obj.onReadyKey());
           
        end
        
        %%
        %Update the mask to the new parameters.
        %
        function obj = updateBkgMask(obj)
            
            mymask = obj.cell_mask;
            
            %Cut pixels below min...
            low_mask = (obj.stdev_proj >= obj.stdev_thresh_min);
            hi_mask = (obj.stdev_proj <= obj.stdev_thresh_max);
            
            mymask = immultiply(low_mask, mymask);
            mymask = immultiply(hi_mask, mymask);
            
            obj.bkg_mask = mymask;
            
            %figure(11111111);
            %clf;
            %imshow(mymask, []);
            
        end
        
        %%
        %Set the default mask save path.
        %
        %ARGS
        %   savepath (string) - Path to set for saving mask from GUI
        %
        function obj = setBkgMaskSavePath(obj, savepath)
            obj.save_path = savepath;
        end
        
        %%
        %Save the background mask to a mat file
        %
        %ARGS
        %   savepath (string) - Path to save mask data to
        %
        function obj = saveMaskTo(obj, savepath)
            background_mask = obj.bkg_mask;
            th_min = obj.stdev_thresh_min;
            th_max = obj.stdev_thresh_max;
            save(savepath, 'background_mask', 'th_min', 'th_max');
        end
        
        %%
        %Save the entire background extractor instance to a mat file
        %
        %ARGS
        %   savepath (string) - Path to save object data to
        %
        function obj = saveMeTo(obj, savepath)
            save(savepath, 'obj');
        end
        
        %%
        %Listener key press callback function - waits for another key or
        %mouse press then deploys next callback.
        %
        function obj = onReadyKey(obj)
            [x,~,but] = ginput(1);
            
            if but == 1 %Left click, set min
                obj.last_thresh = 0;
                newval = x;
                if newval < 0
                    newval = 0;
                end
                if newval >= obj.stdev_thresh_max
                    newval = obj.stdev_thresh_max - 1;
                end
                    
                obj.stdev_thresh_min = newval;
                obj = obj.updateBkgMask();
                obj = obj.renderThSelectFigures();
            elseif but == 2 %Middle click, exit loop
                obj.exit_ok = 1;
                obj = obj.closeFigures();
            elseif but == 3 %Right click, set max
                obj.last_thresh = 1;
                newval = x;
                if newval <= obj.stdev_thresh_min
                    newval = obj.stdev_thresh_min + 1;
                end
                if newval > obj.stdev_max
                    newval = obj.stdev_max;
                end
                    
                obj.stdev_thresh_max = newval;
                obj = obj.updateBkgMask();
                obj = obj.renderThSelectFigures();
            elseif but == 120 %'x' - exit loop
                obj.exit_ok = 1;
                obj = obj.closeFigures();
            elseif but == 115 %'s' - save
                fprintf("Saving mask data to %s...\n", obj.save_path);
                obj.saveMaskTo(obj.save_path);
                fprintf("Mask data saved!\n");
            elseif but == 104 %'h' - try close holes
                obj = obj.closeMaskHoles();
                obj = obj.renderThSelectFigures();
            elseif but == 101 %'e' - try erode
                obj = obj.erodeMask();
                obj = obj.renderThSelectFigures();
            elseif but == 114 %'r' - reset image
                obj = obj.updateBkgMask();
                obj = obj.renderThSelectFigures();
            elseif but == 93 %']' - increase open radius
                obj.open_rad = obj.open_rad + 1;
                obj = obj.closeMaskHoles();
                obj = obj.renderThSelectFigures();
            elseif but == 91 %'[' - decrease open radius
                if obj.open_rad > 0
                    obj.open_rad = obj.open_rad - 1;
                end
                obj = obj.closeMaskHoles();
                obj = obj.renderThSelectFigures();
            elseif but == 125 %'}' - increase erode radius
                obj.erode_rad = obj.erode_rad + 1;
                obj = obj.erodeMask();
                obj = obj.renderThSelectFigures();
            elseif but == 123 %'{' - decrease erode radius
                if obj.erode_rad > 0
                    obj.erode_rad = obj.erode_rad - 1;
                end
                obj = obj.erodeMask();
                obj = obj.renderThSelectFigures();
            elseif but == 43 %'+' - Increase increment size
                %Increment size cannot exceed 1000
                if obj.inc_value < 1000
                    obj.inc_value = obj.inc_value + 1;
                    fprintf("Increment size increased to %d\n", obj.inc_value);
                end
            elseif but == 45 %'-' - Decrease increment size
                %Increment size cannot be below 1
                if obj.inc_value > 1
                    obj.inc_value = obj.inc_value - 1;
                    fprintf("Increment size decreased to %d\n", obj.inc_value);
                end
            elseif but == 62 %'>' - Increment last thresh moved
                if obj.last_thresh
                    newval = obj.stdev_thresh_max + obj.inc_value;
                    if newval > obj.stdev_max
                        newval = obj.stdev_max;
                    end
                    obj.stdev_thresh_max = newval;
                else
                    newval = obj.stdev_thresh_min + obj.inc_value;
                    if newval >= obj.stdev_thresh_max
                        newval = obj.stdev_thresh_max - 1;
                    end
                    obj.stdev_thresh_min = newval;
                end
                obj = obj.updateBkgMask();
                obj = obj.renderThSelectFigures();
            elseif but == 60 %'<' - Decrement last thresh moved
                if obj.last_thresh
                    newval = obj.stdev_thresh_max - obj.inc_value;
                    if newval <= obj.stdev_thresh_min
                        newval = obj.stdev_thresh_min + 1;
                    end
                    obj.stdev_thresh_max = newval;
                else
                    newval = obj.stdev_thresh_min - obj.inc_value;
                    if newval < 0
                        newval = 0;
                    end
                    obj.stdev_thresh_min = newval;
                end
                obj = obj.updateBkgMask();
                obj = obj.renderThSelectFigures();
            end
            
        end
        
        %%
        %Use stored figure handles to close all figures in preparation for
        %exiting.
        %
        function obj = closeFigures(obj)
            img_count = size(obj.img_handles, 1);
            for i = 1:img_count
                close(obj.img_handles{i});
            end
        end
        
        %%
        %Close small holes in current mask using the radius set for this
        %instance and the imopen function.
        %
        function obj = closeMaskHoles(obj)
            obj.bkg_mask = imopen(obj.bkg_mask, strel('disk', obj.open_rad));
            obj.bkg_mask = ~imfill(~obj.bkg_mask, 'holes');
        end
        
        %%
        %Erode current mask using the radius set for this
        %instance and the imerode function.
        %
        function obj = erodeMask(obj)
            obj.bkg_mask = imerode(obj.bkg_mask, strel('disk', obj.erode_rad));
        end
        
        %%
        %Perform an automated contrast enhancement on the currently loaded
        %image stack of this Background_Extractor. This method is mainly
        %for use in automated threshold detection if input image contrast
        %is low.
        %
        function obj = enhanceContrast(obj)
            obj.img_ch = Background_Extractor.enhanceStackContrast(obj.img_ch);
            
            %Generate stdev projection
            obj.stdev_proj = std(obj.img_ch, 0, 3);
            
            %Generate max projection (for viewing)
            obj.max_proj = max(obj.img_ch,[],3);
        end
        
        %%
        %Call closeMaskHoles() which runs a quick imopen on the current
        %mask, print information about the parameters used, and save the
        %mask to the path specified in the object instance variables.
        %
        function obj = applyAndSave(obj)
            %Also print some info on the threshold selected and radius used
           %for reference
           
           fprintf("Saving mask to %s\n", obj.save_path);
           fprintf("Pixel stdevpass range: %d - %d\n", obj.stdev_thresh_min, obj.stdev_thresh_max);
           fprintf("imopen Radius: %d\n", obj.open_rad);
           
           obj = obj.closeMaskHoles();
           obj = obj.saveMaskTo(obj.save_path);
            
        end
        
        %%
        %Render this instance's image with the background mask applied and
        %save the image to disk at the specified path.
        %
        %ARGS
        %   outpath (string) - Path on local file system to save image to
        %
        function obj = exportMaskedImage(obj, outpath)
            masked_proj = immultiply(obj.max_proj, obj.bkg_mask);
            fig = figure(1);
            imshow(masked_proj, []);
            saveas(fig, outpath);
            close(fig);
            %imwrite(masked_proj, outpath);
        end
        
        %%
        %(Debug function)
        %Render this a histogram derived from the max z projection of this
        %instance's associated image. 
        %
        %ARGS
        %   outpath (string) - Path on local file system to save image to
        %
        function obj = exportHistogramImage(obj, outpath)
            
            color3 = [0.608 0.621 0.628];
            dim = [.45, .4, .5, .5];
            
            stdev16 = uint16(obj.stdev_proj);
            max_stdev = max(max(stdev16(:,:)));
            obj.stdev_max = max_stdev;
            if obj.stdev_thresh_max > max_stdev
                obj.stdev_thresh_max = max_stdev;
            end
            th_min = obj.stdev_thresh_min;
            th_max = obj.stdev_thresh_max;
            fig = figure(1);
            clf;
            ax = axes;
            [counts,~] = imhist(stdev16, max_stdev);
            imhist(stdev16, max_stdev);
            hold on;
            max_bin_sz = max(counts);
            axis([0 max_stdev 0 max_bin_sz]);
            title("STDEV Histogram");
            line([th_min th_min], get(ax,'YLim'),'Color',color3,'LineStyle','--');
            line([th_max th_max], get(ax,'YLim'),'Color',color3,'LineStyle','--');
            t1_str = ['Minimum: ' num2str(th_min)];
            t2_str = ['Maximum: ' num2str(th_max)];
            anno_str = {t1_str, t2_str};
            annotation('textbox',dim,'String',anno_str,'FitBoxToText','on');
            
            
            saveas(fig, outpath);
            close(fig);
        end
        
    end
    
    methods(Static)
        
        %%
        %Load a Background_Extractor instance from a saved mat file.
        %
        %ARGS
        %   class_savepath (string) - Path instance was previously saved to
        %
        function obj = loadExtractorFrom(class_savepath)
            load(class_savepath, 'obj');
        end
        
        %%
        %Load mask data previously saved to a mat file.
        %
        %ARGS
        %   savepath (string) - Path mask was previously saved to
        %
        %RETURN
        %   bkg_mask (boolean[Y][X]) - Loaded mask as a 2D logical matrix
        %
        function bkg_mask = loadMaskFrom(savepath)
            load(savepath, 'background_mask');
            bkg_mask = background_mask;
        end
        
        %%
        %Attempt to auto-detect a background standard deviation threshold
        %(stdev of xy pixel across z) with a bimodal stdev distribution. 
        %
        %ARGS
        %   stdev16 (uint16[Y][X]) - StDev projection of image channel 
        %       (scaled to 16-bit unsigned ints) across z
        %   max_stdev (num) - Maximum standard devation value to use for
        %       scaling
        %
        %RETURN
        %   thresh (num) - Suggested auto-detected stdev threshold
        %       below which a pixel can be called "background"
        %       This method returns -1 if a good threshold was not
        %       detected.
        %
        function thresh = autoestStdevThreshold_bimodal(stdev16, max_stdev)  
            %Returns -1 if the distribution isn't bimodal
            
            %binct = max_stdev/4;
            binct = max_stdev;
            binsz = 65535/(binct-1);
            [histo_counts] = imhist(stdev16, binct);
            
            %LPF
            %lowpass_thresh = 0.25;
            %smooth = lowpass(histo_counts, lowpass_thrhjnnnn/esh); Looks
            %like one of my cats had a good time...
            
            %Smooth
            sm = smooth(histo_counts');
            sm = smooth(sm);
            
            %floor to 0
            sm(sm < 0) = 0;
            % for i=1:size(sm,1)
            %     if sm(i) < 0
            %         sm(i) = 0;
            %     end
            % end
            %transpose(1:size(histo_counts))
            %smgauss = fit(transpose(1:size(histo_counts)), smooth, 'gauss2');
            
            %Plot smoothed function -- Debug only
            %figure(88888);
            %plot(histo_counts);
            %hold on;
            %plot(smooth);
            %plot(raw_counts);
            %max_bin_sz = max(histo_counts);
            %axis([0 binct 0 max_bin_sz]);
            %title("STDEV Histogram (4X Bins, Smoothed)");
            %legend({'Quad Bins', 'Quad Bins + Smooth', 'Raw'})
            %title("STDEV Histogram");
            %legend({'Raw', 'Smoothed', 'Derivative'})
            
            %Deriv
            deriv = diff(sm);
            
            %Plot deriv approx. -- Debug only
            %figure(99999);
            %plot(sm);
            %hold on;
            %plot(deriv);
            
            %Find two maxes and a min
            max1_val = -1;
            max1 = -1;
            max2 = -1;
            min1 = -1;
            
            ptcount = size(deriv,2);
            for i=2:ptcount
                if (deriv(i-1) > 0) && (deriv(i) < 0)
                    %Maximum
                    if max1 == -1
                        max1 = i;
                        max1_val = sm(i);
                    elseif max2 == -1
                        if(sm(i) > (max1_val/20)) %second peak must be at least 5% height of first
                            max2 = i;
                        end
                    end
                elseif (deriv(i-1) < 0) && (deriv(i) > 0)
                    %Minimum
                    if (max1 ~= -1) && (min1 == -1)
                        min1 = i;
                    end
                end
            end
            
            %Map (Debug)
           % color3 = [0.608 0.621 0.628];
            %line([max1 max1], [0 max_bin_sz],'Color',color3,'LineStyle','--');
            %line([min1 min1], [0 max_bin_sz],'Color',color3,'LineStyle','--');
            %line([max2 max2], [0 max_bin_sz],'Color',color3,'LineStyle','--');
            %max1
            %min1
            %max2
            
            %Return -1 if 
            if max2 == -1
                thresh = -1;
                return;
            end
            if min1 == -1
                thresh = -1;
                return;
            end
            
            thresh = min1 * binsz;
            %max1
            %min1
            %max2
            
        end
        
        %%
        %Attempt to auto-detect a background standard deviation threshold
        %(stdev of xy pixel across z) with a low-weighted stdev distribution. 
        %
        %ARGS
        %   stdev16 (uint16[Y][X]) - StDev projection of image channel 
        %       (scaled to 16-bit unsigned ints) across z
        %   max_stdev (num) - Maximum standard devation value to use for
        %       scaling
        %
        %RETURN
        %   thresh (num) - Suggested auto-detected stdev threshold
        %       below which a pixel can be called "background"
        %       This method returns -1 if a good threshold was not
        %       detected.
        %
        function thresh = autoestStdevThreshold_lowweight(stdev16, max_stdev)  

            binct = max_stdev/10;
            binsz = 65535/(binct-1);
            [histo_counts] = imhist(stdev16, binct);
            
            %Smooth
            sm = smooth(histo_counts');
            sm = smooth(sm);
            
            %floor to 0
            sm(sm < 0) = 0;
            
            %Deriv
            deriv1 = smooth(diff(sm));
            deriv2 = smooth(diff(deriv1));
            
            %Plot -- Debug only
            %figure(88888);
            %plot(histo_counts);
            %hold on;
            %plot(sm);
            %plot(deriv1);
            %plot(deriv2);
            %max_bin_sz = max(histo_counts);
            %axis([0 binct 0 max_bin_sz]);
            %title("STDEV Histogram");
            %legend({'Raw', 'Smoothed', '1st Derivative', '2nd Derivative'})
            
            %Find two maxes and a min
            max1 = -1;
            min1 = -1;
            pttarg = 0.0;
            mypt = -1;
            
            ptcount = size(deriv1,2);
            for i=2:ptcount
                if (deriv1(i-1) > 0) && (deriv1(i) <= 0)
                    %Maximum
                    if max1 == -1
                        max1 = i;
                    end
                end
                if (deriv2(i-1) < 0) && (deriv2(i) >= 0)
                    %Minimum in 1st deriv after max in function
                    if (max1 ~= -1) && (min1 == -1)
                        min1 = i;
                        pttarg = deriv1(i)/20;
                        %pttarg = -20;
                    end
                end
                if (min1 ~= -1) && (mypt == -1)
                    if deriv1(i) >= pttarg
                        mypt = i;
                        %mypt
                    end
                    %if deriv2(i) < pttarg
                    %    mypt = i;
                    %end
                end
            end

            if mypt == -1
                thresh = 0;
                return;
            end
            
            thresh = mypt * binsz;

        end
        
        %%
        %Generate a copy of the input image with the contrast enhanced
        %using imadjust.
        %
        %ARGS
        %   inimg (num[Y][X][Z]) - Input image as 3D matrix
        %
        %RETURN
        %   outimg (num[Y][X][Z]) - Output image as 3D matrix
        %
        function outimg = enhanceStackContrast(inimg)
            X = size(inimg, 2);
            Y = size(inimg, 1);
            Z = size(inimg, 3);
            outimg = NaN(Y, X, Z);
            
            %maxraw = max(inimg, [], 'all');
            %maxval = maxraw/65536.0;
            
            for z = 1:Z
                %outimg(:,:,z) = adapthisteq(inimg(:,:,z), 'ClipLimit', 0.9, 'NBins', 1024);
                %outimg(:,:,z) = imadjust(inimg(:,:,z), stretchlim(inimg(:,:,z)), [0 1], 0.9);
                nplane = inimg(:,:,z)/65535.0;
                outimg(:,:,z) = imadjust(nplane, stretchlim(nplane), [0 1], 0.5);
                %outimg(:,:,z) = imadjust(nplane, [0 maxval], [0 1], 0.5);
            end
            outimg = outimg * 65536.0;
            outimg = RNA_Threshold_Common.applyGaussianFilter(outimg,13,2);
        end
        
    end
    
end