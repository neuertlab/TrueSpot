%Module for visualizing RNA spot counts at various thresholds for an image
%upon which spot detection has been previously run.
%Blythe Hospelhorn
%Modified from code written by Ben Kesler & Gregor Neuert
%Version 1.0.0
%Updated Mar 11, 2021

%Update Log:
%   1.0.0 | 21.03.11
%       Init doc
%   1.1.0 | 21.09.01
%       Modified so can synchronously print plots to output without user
%       input


%%
classdef RNA_Threshold_Plotter
     
    %GUI Controls:
            %   s - Save
            %   Left Click - Select threshold to view
            %   Right Click - Exit
    
    %%
    methods (Static)
        
        %%
        %Generate set of plots and image viewer from pre-processed stack
        %and spot detection data as an interactive MATLAB GUI.
        %This method calls the other methods in this class, so from the
        %outside, this is the only one you will probably need to call.
        %
        %ARGS
        %   save_stem_signal (string) - Path prefix for spot detection
        %       output data of signal channel to use as plotter input.
        %   save_stem_ctrls (cell[N-1]) - Path prefixes (stored as an array
        %       of string cells) for spot detection output data of all
        %       control channels/images to use.
        %   probeNames (cell[N]) - List of probe names for all samples,
        %       including signal and controls (with signal first), stored
        %       as ab array of string cells, to use for labeling.
        %   st_thresh (int) - Initial threshold value to mark plot and
        %       render projections with
        %   img_output_dir (string) - Path to directory to output any data
        %       or images saved from GUI to.
        %   count_from_coords (bool) - Whether to recount spots at each
        %       threshold from the coordinate tables (true) or just use the
        %       spot count tables (false)
        %   zmin (int) - Index of lowest z slice to include in count. Only
        %       used if coord recount is done.
        %   zmax (int) - Index of highest z slice to include in count. Only
        %       used if coord recount is done.
        %   async (bool) - Whether to launch plots as asynchronous (wait
        %       for user input) or not (render and save plots then return)
        %
        %EXAMPLE INPUT (Pseudosyntax)
        %   save_stem_signal = "D:\Users\Me\Science\RNA1\Signal3D"
        %   save_stem_ctrls = {"D:\Users\Me\Science\RNA1\SCR3D",
        %                      "D:\Users\Me\Science\RNA1\BKG3D",
        %                      "D:\Users\Me\Science\RNA1\NoProbe3D"}
        %   probeNames = {"My RNA",
        %                 "Scramble",
        %                 "Background",
        %                 "No Probe"}
        %   st_thresh = 50
        %   img_output_dir = "D:\Users\Me\Science\RNA1\SpotPlots"
        %   count_from_coords = false
        %   zmin = 1
        %   zmax = 69
        %   async = true
        %
        function plotPreprocessedData(save_stem_signal, save_stem_ctrls, probeNames,...
                st_thresh, img_output_dir, count_from_coords, zmin, zmax, async)
            
            async = Force2Bool(async);
            
            ctrl_count = size(save_stem_ctrls, 1);
            idx_rna = 1;
            
            spot_table_suffix = '_spotTable.mat';
            tbl_path_RNA = [save_stem_signal spot_table_suffix];

            load(tbl_path_RNA, 'spot_table');
            th_table = spot_table(:,1);
            reg_tables(idx_rna).spot_table = double(spot_table(:,2));
            %reg_tables(idx_rna).deriv1 = diff(reg_tables(idx_rna).spot_table);
            %reg_tables(idx_rna).deriv1 = smooth(diff(reg_tables(idx_rna).spot_table));
            for c = 1:ctrl_count
                %fprintf("-DEBUG- Loading: %s...\n", [save_stem_ctrls{c} spot_table_suffix]);
                tidx = c+1;
                load([save_stem_ctrls{c} spot_table_suffix], 'spot_table');
                reg_tables(tidx).spot_table = double(spot_table(:,2));
                %reg_tables(tidx).deriv1 = smooth(diff(reg_tables(tidx).spot_table));
            end
            

            %Load coordinates
            coord_table_suffix = '_coordTable.mat';
            tbl_path_RNA = [save_stem_signal coord_table_suffix];

            load(tbl_path_RNA, 'coord_table');
            reg_tables(idx_rna).ctable = coord_table;
            for c = 1:ctrl_count
                %fprintf("-DEBUG- Loading: %s...\n", [save_stem_ctrls{c} coord_table_suffix]);
                load([save_stem_ctrls{c} coord_table_suffix], 'coord_table');
                reg_tables(c+1).ctable = coord_table;
            end
            
            if count_from_coords
                %Adjust counts to match coordinates
                tcount = ctrl_count + 1;
                T = size(th_table,1);
                for c = 1:tcount
                    coordtbl = reg_tables(c).ctable;
                    for t = 1:T
                        ctt = coordtbl{t};
                        ctt_sz = size(ctt,1);
                        count = 0;
                        if size(ctt,2) > 2
                            for s = 1:ctt_sz
                                if ctt(s,3) >= zmin
                                    if ctt(s,3) <= zmax
                                        count = count+1;
                                    end
                                end
                            end
                        else
                            fprintf("WARNING: Coordinate table for sample %d (threshold %d) is not 3D!\n",c, t);
                            count = ctt_sz;
                        end
                        reg_tables(c).spot_table(t) = count;
                    end
                end
            end
            
            %Derivs
            totalcount = ctrl_count + 1;
            for c = 1:totalcount
                %I put the smooths on different lines so they can be easily
                %commented out if needed
                d1 = diff(reg_tables(c).spot_table);
                d1 = smooth(d1);
                d2 = diff(d1);
                d2 = smooth(d2);
                reg_tables(c).deriv1 = d1;
                reg_tables(c).deriv2 = d2;
            end

            %Load and prepare images
            imgstructs_suffix = '_imgviewstructs';
            path_RNA = [save_stem_signal imgstructs_suffix];

            load(path_RNA, 'my_images');
            all_imgs(1, 1) = my_images(1);
            all_imgs(1, 2) = my_images(2);
            for c = 1:ctrl_count
                %fprintf("-DEBUG- Loading: %s...\n", [save_stem_ctrls{c} imgstructs_suffix]);
                load([save_stem_ctrls{c} imgstructs_suffix], 'my_images');
                all_imgs(c+1, 1) = my_images(1);
                all_imgs(c+1, 2) = my_images(2);
            end

            %Differences & Ratios
            %for c = 1:ctrl_count
            %    reg_tables(c+1).diffs = abs(reg_tables(idx_rna).spot_table(:) - reg_tables(c+1).spot_table(:));
            %    reg_tables(c+1).ratio = abs(reg_tables(idx_rna).spot_table(:)./reg_tables(c+1).spot_table(:));
            %end

            T = size(th_table,1);
            Tmax = th_table(T);
            if st_thresh < 1
                th_idx = int16(T./2);
            else
                th_i_list = find(th_table(:) == st_thresh,1);
                th_idx = th_i_list(1);
            end
            %th_idx = T./2;
            %th = th_table(th_idx);

            %Logs
            %log_tables(idx_rna).spots = log(reg_tables(idx_rna).spot_table(:));
            log_tables(idx_rna).spots = log10(reg_tables(idx_rna).spot_table(:));
            for c = 1:ctrl_count
                %log_tables(c+1).spots = log(reg_tables(c+1).spot_table(:));
                log_tables(c+1).spots = log10(reg_tables(c+1).spot_table(:));
                %log_tables(c+1).diffs = log(reg_tables(c+1).diffs(:));
                %log_tables(c+1).ratio = log(reg_tables(c+1).ratio(:));
            end
            
            img_handles = RNA_Threshold_Plotter.drawPlots(th_table, reg_tables, log_tables, th_idx, probeNames, all_imgs);

            if async
                buu = 5;
                while buu == 5
                    [x,~,but] = ginput(1);
                    if but == 1 %Left click, set thresh
                        %fprintf("Left click! x = %f\n", x)
                        if x > 0.0 && x < Tmax
                            rounded = round(x);
                            th_i_list = find(th_table(:) == rounded,1);
                            th_idx = th_i_list(1);
                            %th = spot_table_tsix(th_idx,1);
                            img_handles = RNA_Threshold_Plotter.drawPlots(th_table, reg_tables, log_tables, th_idx, probeNames, all_imgs);
                        end
                    elseif but == 2 %Middle click, exit loop
                        buu = 2;
                    elseif but == 3 %Right click, exit loop
                        buu = 2;
                    elseif but == 115 %"s" - save
                        fprintf("Saving plots and images to %s...\n", img_output_dir)
                        mkdir(img_output_dir);
                        %Normal plot
                        saveas(img_handles{1}, [img_output_dir filesep 'SpotPlot_th_' num2str(th_table(th_idx)) '.png']);
                        saveas(img_handles{1}, [img_output_dir filesep 'SpotPlot_th_' num2str(th_table(th_idx))], 'epsc');
                        %Log plot
                        saveas(img_handles{2}, [img_output_dir filesep 'LogSpotPlot_th_' num2str(th_table(th_idx)) '.png']);
                        saveas(img_handles{2}, [img_output_dir filesep 'LogSpotPlot_th_' num2str(th_table(th_idx))], 'epsc');
                        %Images
                        img_count = ctrl_count+1;
                        for i = 1:img_count
                            istruct(1) = all_imgs(i,1);
                            istruct(2) = all_imgs(i,2);
                            out_path = [img_output_dir filesep probeNames{i} '_filteredIMG_th_' num2str(th_table(th_idx))];
                            RNA_Threshold_Plotter.saveFigureImage(istruct, 1, reg_tables(i).ctable{th_idx}, out_path);
                            out_path = [img_output_dir filesep probeNames{i} '_rawIMG_th_' num2str(th_table(th_idx))];
                            RNA_Threshold_Plotter.saveFigureImage(istruct, 2, reg_tables(i).ctable{th_idx}, out_path);
                        end

                        fprintf("Plots saved!\n")
                    end
                end
            else
                %Just render and save.
                mkdir(img_output_dir);
                %Normal plot
                saveas(img_handles{1}, [img_output_dir filesep 'SpotPlot_th_' num2str(th_table(th_idx)) '.png']);
                saveas(img_handles{1}, [img_output_dir filesep 'SpotPlot_th_' num2str(th_table(th_idx))], 'epsc');
                %Log plot
                saveas(img_handles{2}, [img_output_dir filesep 'LogSpotPlot_th_' num2str(th_table(th_idx)) '.png']);
                saveas(img_handles{2}, [img_output_dir filesep 'LogSpotPlot_th_' num2str(th_table(th_idx))], 'epsc');
                %Images
                img_count = ctrl_count+1;
                for i = 1:img_count
                    istruct(1) = all_imgs(i,1);
                    istruct(2) = all_imgs(i,2);
                    out_path = [img_output_dir filesep probeNames{i} '_filteredIMG_th_' num2str(th_table(th_idx))];
                    RNA_Threshold_Plotter.saveFigureImage(istruct, 1, reg_tables(i).ctable{th_idx}, out_path);
                    out_path = [img_output_dir filesep probeNames{i} '_rawIMG_th_' num2str(th_table(th_idx))];
                    RNA_Threshold_Plotter.saveFigureImage(istruct, 2, reg_tables(i).ctable{th_idx}, out_path);
                end
            end
            
            %Delete figures
            img_count = size(img_handles, 1);
            for i = 1:img_count
                close(img_handles{i});
            end
            
        end
        
        %%
        %Draw plots of RNA spot count vs. threshold & ln(spot count) vs.
        %threshold for all requested samples and render to MATLAB figures.
        %Additionally render 2D max projections of images with spots
        %circled for reference.
        %
        %This is an internal method - best if only called by other class
        %methods.
        %
        %ARGS
        %   th_table (int[T]) - Array of threshold values 
        %   reg_tables (struct[N]) - Data for non-log scaled plots
        %   log_tables (struct[N]) - Data for log-scaled plots
        %   th_idx (int) - Index of threshold in threshold table to use as
        %       the focus to draw images and plots
        %   probeNames (cell<string>[N]) - List of probe names for all samples 
        %   all_images (cell<img_struct>[N,2]) - Data containing image
        %       projections and scaling values for each probe
        %
        %RETURN
        %   img_handles (cell<figure_handle>[N+2]) - Handles to figures
        %       rendered by the plotter
        %
        function img_handles = drawPlots(th_table, reg_tables, log_tables, th_idx, probeNames, all_images)
            
            img_count = size(reg_tables, 2);
            img_handles = cell(img_count + 2 + 2,1);
            th = th_table(th_idx);
            
            idx_rna = 1;
            
            %fprintf("-DEBUG- Image count: %d...\n", img_count);
            
            color1 = [0.849 0.633 0.778]; %#F2BBE0
            color2 = [0.759 0.665 0.802]; %#DBC3E6
            color3 = [0.416 0.427 0.678]; %#6A6DAD
            color4 = [0.416 0.588 0.678]; %#6A96AD
            grey = [0.608 0.621 0.628]; %#A8ABAD
            dim = [.45, .4, .5, .5];
            line_colors = {color1;
                           color2; 
                           color3; 
                           color4};
            
            %Draw images
            for i = 1:img_count
                istruct(1) = all_images(i,1);
                istruct(2) = all_images(i,2);
                img_handles{4+i} = RNA_Threshold_Plotter.draw_images(istruct, 2+i, reg_tables(i).ctable{th_idx}, [probeNames{i} ' Image | Threshold = ' num2str(th)]);
            end
            
            %TODO: Maybe resize to make plots a bit more horizontal
            %plot_count = (img_count.*3) - 2;
            plot_count = img_count;
            legend_names = cell(1, plot_count);
            plots = cell(plot_count);
            img_handles{1} = figure(989);
            clf;
            ax = axes;
            plots{1} = plot(th_table(:),reg_tables(idx_rna).spot_table(:),'LineWidth',2);
            hold on;
            legend_names{1,1} = probeNames{1};
            idx = 2;
            for i = 2:img_count
                plots{i} = plot(th_table(:),reg_tables(i).spot_table(:),'-.','LineWidth',2);
                legend_names{1,i} = probeNames{i};
                idx = idx + 1;
            end
            %for i = 2:img_count
            %    plots{idx} = plot(th_table(:),reg_tables(i).diffs(:),'Color',line_colors{i-1},'LineStyle',':','LineWidth',2);
            %    legend_names{1,idx} = [probeNames{1} ' - ' probeNames{i}];
            %    idx = idx + 1;
            %end
            %for i = 2:img_count
            %    plots{idx} = plot(th_table(:),reg_tables(i).ratio(:),'Color',line_colors{i-1},'LineStyle','--','LineWidth',2);
            %    legend_names{1,idx} = [probeNames{1} ' / ' probeNames{i}];
            %    idx = idx + 1;
            %end
            line([th th], get(ax,'YLim'),'Color',grey,'LineStyle','--');
            legend(legend_names);
            
            %Annotation showing values
            anno_str = cell(1, img_count+2);
            anno_str{1, 1} = ['Threshold: ' num2str(th)];
            anno_str{1, 2} = '';
            for i = 1:img_count
                anno_str{1, i+2} = [probeNames{i} ' Spots: ' num2str(reg_tables(i).spot_table(th_idx))];
            end
            annotation('textbox',dim,'String',anno_str,'FitBoxToText','on');
            xlabel('Threshold');
            ylabel('# Spots Detected');
            
            %Derivatives...
            thcount = size(th_table,1);
            legend_names = cell(1, img_count);
            plots = zeros(img_count);
            img_handles{3} = figure(888);
            clf;
            ax = axes;
            plots(1) = plot(th_table(1:thcount-1),reg_tables(1).deriv1(:),'LineWidth',2);
            hold on;
            legend_names{1,1} = ['d(' probeNames{1} ')/dx'];
            idx = 2;
            for i = 2:img_count
                plots(i) = plot(th_table(1:thcount-1),reg_tables(i).deriv1(:),'-.','LineWidth',2);
                legend_names{1,i} = ['d(' probeNames{i} ')/dx'];
                idx = idx + 1;
            end
            line([th th], get(ax,'YLim'),'Color',grey,'LineStyle','--');
            legend(legend_names);
            xlabel('Threshold');
            ylabel('d(# Spots Detected)/dx');
            
            anno_str = cell(1, img_count+2);
            anno_str{1, 1} = ['Threshold: ' num2str(th)];
            anno_str{1, 2} = '';
            for i = 1:img_count
                anno_str{1, i+2} = [probeNames{i} ' d(Spots)/dx: ' num2str(reg_tables(i).deriv1(th_idx))];
            end
            annotation('textbox',dim,'String',anno_str,'FitBoxToText','on');
            
            %And second deriv...
            legend_names = cell(1, img_count);
            plots = zeros(img_count);
            img_handles{4} = figure(899);
            clf;
            ax = axes;
            plots(1) = plot(th_table(1:thcount-2),reg_tables(1).deriv2(:),'LineWidth',2);
            hold on;
            legend_names{1,1} = ['d2(' probeNames{1} ')/dx'];
            idx = 2;
            for i = 2:img_count
                plots(i) = plot(th_table(1:thcount-2),reg_tables(i).deriv2(:),'-.','LineWidth',2);
                legend_names{1,i} = ['d2(' probeNames{i} ')/dx'];
                idx = idx + 1;
            end
            line([th th], get(ax,'YLim'),'Color',grey,'LineStyle','--');
            legend(legend_names);
            xlabel('Threshold');
            ylabel('d2(# Spots Detected)/dx');
            
            anno_str = cell(1, img_count+2);
            anno_str{1, 1} = ['Threshold: ' num2str(th)];
            anno_str{1, 2} = '';
            for i = 1:img_count
                anno_str{1, i+2} = [probeNames{i} ' d2(Spots)/dx: ' num2str(reg_tables(i).deriv1(th_idx))];
            end
            annotation('textbox',dim,'String',anno_str,'FitBoxToText','on');

            %log plot
            legend_names = cell(1, plot_count);
            plots = zeros(plot_count);
            img_handles{2} = figure(898);
            clf;
            ax = axes;
            plots(1) = plot(th_table(:),log_tables(idx_rna).spots(:),'LineWidth',2);
            hold on;
            legend_names{1,1} = ['log10(' probeNames{1} ')'];
            idx = 2;
            for i = 2:img_count
                plots(i) = plot(th_table(:),log_tables(i).spots(:),'-.','LineWidth',2);
                legend_names{1,i} = ['log10(' probeNames{i} ')'];
                idx = idx + 1;
            end
            %for i = 2:img_count
            %    plots(idx) = plot(th_table(:),log_tables(i).diffs(:),'Color',line_colors{i-1},'LineStyle',':','LineWidth',2);
            %    legend_names{1,idx} = ['ln(' probeNames{1} ' - ' probeNames{i} ')'];
            %    idx = idx + 1;
            %end
            %for i = 2:img_count
            %    plots(idx) = plot(th_table(:),log_tables(i).ratio(:),'Color',line_colors{i-1},'LineStyle','--','LineWidth',2);
            %    legend_names{1,idx} = ['ln(' probeNames{1} ' / ' probeNames{i} ')'];
            %    idx = idx + 1;
            %end
            line([th th], get(ax,'YLim'),'Color',grey,'LineStyle','--');
            legend(legend_names);
            xlabel('Threshold');
            ylabel('log10(# Spots Detected)');
            
        end
        
        %%
        %Draw MATLAB figure of the images given and circle spot locations.
        %
        %This is an internal method - best if only called by other class
        %methods.
        %
        %ARGS
        %   istruct (image_struct[2]) - Image data (for unfiltered and
        %       filtered max z projection)
        %   figureNumber (int) - Number to use for figure handle
        %   spot_coords (int[S][2+]) - Spot coordinate table (x,y)
        %   title_str (string) - String to use to title figure
        %
        %RETURN
        %   f_handle (fig_handle) - Handle of generated MATLAB figure
        %
        function f_handle = draw_images(istruct, figureNumber, spot_coords, title_str)
            %Draw images
            f_handle = figure(figureNumber);
            clf;
    
            %Render left image (usually filtered image) w/ spots circled
            subplot(1,2,1);
            img1 = istruct(1);
            imshow(img1.image,[img1.Lmin img1.Lmax]);
            hold on;
            if(~isempty(spot_coords))
                plot(spot_coords(:,1), spot_coords(:,2),'or','markersize',10);
            end
            %title([' plane: ' num2str(planeNumber) ', th: ' num2str(threshold)]);
            impixelinfo;
    
            %Render right image (usually max projection) w/ spots circled
            subplot(1,2,2);
            img2 = istruct(2);
            imshow(img2.image,[img2.Lmin img2.Lmax]);
            hold on;
            if(~isempty(spot_coords))
                plot(spot_coords(:,1), spot_coords(:,2),'or','markersize',10);
            end
            %title([' plane: ' num2str(planeNumber) ', th: ' num2str(threshold)]);
            impixelinfo;
    
            sgtitle(title_str);
    
        end
        
        %%
        %Save the specified images with spots circled to the specified path
        %on local file system.
        %
        %This is an internal method - best if only called by other class
        %methods.
        %
        %ARGS
        %   istruct (image_struct[2]) - Image data (for unfiltered and
        %       filtered max z projection)
        %   istruct_idx(int) - Index in istruct array of image to use
        %   spot_coords (int[S][2+]) - Spot coordinate table (x,y)
        %   out_path(string) - Image output path stem
        %
        function saveFigureImage(istruct, istruct_idx, spot_coords, out_path)
            
            f_handle = figure(4444);
            clf
            
            my_img = istruct(istruct_idx);
            imshow(my_img.image,[my_img.Lmin my_img.Lmax]);
            hold on;
            plot(spot_coords(:,1), spot_coords(:,2),'or','markersize',10);
            
            saveas(f_handle, [out_path '.png']);
            saveas(f_handle, out_path, 'epsc');
            
            close(f_handle);
            
        end

    end
end