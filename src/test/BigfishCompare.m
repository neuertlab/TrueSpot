%
%%

classdef BigfishCompare

    methods (Static)
        
        % ========================== Process ==========================
        
        function csv_path = getCSVPath(dirpath, thresh_val)
            csv_path = [dirpath filesep 'spots_' sprintf('%04d', thresh_val) '.csv'];
        end

        function [coord_table, spot_table] = importBigFishCsvs(bfdir, output_stem, z_offset)
            %Figure out t range to preallocate.
            t = 1;
            filetblname = BigfishCompare.getCSVPath(bfdir, t);
            while ~isfile(filetblname)
                t = t+1;
                filetblname = BigfishCompare.getCSVPath(bfdir, t);
            end
            tmin = t;
            
            while isfile(filetblname)
                t = t+1;
                filetblname = BigfishCompare.getCSVPath(bfdir, t);
            end
            tmax = t - 1;
            T = tmax - tmin + 1;
            
            %Preallocate structures
            spot_table = zeros(T,2);
            coord_table = cell(T,1);
            spot_table(:,1) = [tmin:1:tmax];
            
            %TODO It's complaining about something being out of bounds.
            %WHAT that is, who the fuck knows
            %Import csvs
            i = 1;
            for t = tmin:tmax
                filetblname = BigfishCompare.getCSVPath(bfdir, t);
                if isfile(filetblname)
                    
                    %If it's empty, skip.
                    finfo = dir(filetblname);
                    if finfo.bytes < 1
                        break;
                    end
                    
                    %Load in.
                    raw_coord_table = uint16(csvread(filetblname));
                    spot_table(i,1) = t;
                    spot_table(i,2) = size(raw_coord_table,1);
                    scount = spot_table(i,2);
        
                    tcoords = uint16(zeros(scount,3));
                    tcoords(:,1) = raw_coord_table(:,3) + 1;
                    tcoords(:,2) = raw_coord_table(:,2) + 1;
                    tcoords(:,3) = raw_coord_table(:,1) + z_offset;
        
                    coord_table{i} = tcoords;
                else
                    break;
                end
                i = i+1;
            end
            
            %Save to mat files
            save([output_stem '_coordTable'], 'coord_table', '-v7.3');
            save([output_stem '_spotTable'], 'spot_table');
        end

        function annoobj = loadSpotSelector(stem, ref_stem, bf_spot_table, bf_coord_table, overwrite)
            
            %UPDATE SNAP STOPAT to make sure it is not above threshold
            snapstop = 20;
            [bfdir, ~, ~] = fileparts(stem);
            sumpath = [bfdir filesep 'summary.txt'];
            if isfile(sumpath)
                [~, ~, bfthresh] = BigfishCompare.readSummaryTxt(sumpath);
                thidx = RNAUtils.findThresholdIndex(bfthresh, transpose(bf_spot_table(:,1)));
                if thidx < snapstop; snapstop = thidx; end
            end
           
            %Checks if needs to make a new one or load existing one.
            %Needs to be able to check if ref has been updated since last
            %snap...
            if ~overwrite & RNA_Threshold_SpotSelector.selectorExists(stem)
                %Open
                annoobj = RNA_Threshold_SpotSelector.openSelector(stem, true);
                
                %Update refset if needed.
                ts_me = annoobj.ref_last_modified;
                ts_ref = RNA_Threshold_SpotSelector.getRefsetTimestamp(ref_stem);
                if ~isempty(ts_ref)
                    if ts_ref > ts_me
                        %Needs to be updated...
                        fprintf("Truthset needs to be updated...\n");
                        annoobj.toggle_del_unsnapped = false;
                        
                        src_selector = RNA_Threshold_SpotSelector.openSelector(ref_stem, true);
                        annoobj.ref_coords = src_selector.ref_coords;
                        clear src_selector;
                        
                        %All z toggle - OFF. (Should only snap to
                        %NEARBY z options!!)
                        annoobj.toggle_singleSlice = true;
                        annoobj.toggle_allz = ~annoobj.toggle_singleSlice;
                        annoobj = annoobj.refSnapToAutoSpots(snapstop);
                        annoobj = annoobj.updateFTable();
                        annoobj.ref_last_modified = datetime;
                    end
                else
                    fprintf("Truthset needs to be updated...\n");
                    annoobj.toggle_del_unsnapped = false;
                    
                    src_selector = RNA_Threshold_SpotSelector.openSelector(ref_stem, true);
                    annoobj.ref_coords = src_selector.ref_coords;
                    clear src_selector;

                    annoobj.toggle_singleSlice = true;
                    annoobj.toggle_allz = ~annoobj.toggle_singleSlice;
                    annoobj = annoobj.refSnapToAutoSpots(snapstop);
                    annoobj = annoobj.updateFTable();
                    annoobj.ref_last_modified = datetime;
                end
                
                if annoobj.f_scores_dirty
                    annoobj = annoobj.updateFTable();
                end
            else
                src_selector = RNA_Threshold_SpotSelector.openSelector(ref_stem, true);
                [~, annoobj] = src_selector.makeCopy();
                clear src_selector;
                
                %Create new one
                annoobj.save_stem = stem;
                annoobj = annoobj.loadNewSpotset(bf_spot_table, bf_coord_table);
                
                %Update refset
                annoobj.toggle_del_unsnapped = false;
                annoobj.toggle_singleSlice = true;
                annoobj.toggle_allz = ~annoobj.toggle_singleSlice;
                annoobj = annoobj.refSnapToAutoSpots(snapstop);
                annoobj = annoobj.updateFTable();
                annoobj.ref_last_modified = datetime;
                
                %Return
            end
        end

        function thresh_res = updateThresholdResults(spots_table, ref_spots_run)
            if ~isempty(ref_spots_run)
                hbbf_params = RNAThreshold.paramsFromSpotsrun(ref_spots_run);
                hbbf_params.sample_spot_table = spots_table;
                thresh_res = RNA_Threshold_Common.estimateThreshold(hbbf_params);
            else
                thresh_res = RNAThreshold.runDefaultParameters(spot_counts);
            end
        end
        
        function [zmin, zmax, bfthresh] = readSummaryTxt(filepath)
            %fgetl
            fhandle = fopen(filepath);
            line = fgetl(fhandle);
            sspl = split(line, ':');
            sspl = split(sspl{2}, '-');
            zmin = str2num(strtrim(sspl{1})) + 1;
            zmax = str2num(strtrim(sspl{2}));

            line = fgetl(fhandle);
            sspl = split(line, ':');
            bfthresh = round(str2num(strtrim(sspl{2})));
            fclose(fhandle);
        end
        
        function coord_table = fixCoordinates(zmin, zmax, X, Y, coord_table, bf_stem)
            T = size(coord_table,1);
            
            xoff = 0;
            yoff = 0;
            zoff = 0;
            
            %1 find offsets
            for t = 1:T
                ttbl = coord_table{t,1};
                S = size(ttbl,1);
                for s = 1:S
                    if ttbl(s,1) < 1
                        offset = 1 - ttbl(s,1);
                        if offset > abs(xoff)
                            xoff = offset;
                        end
                    end
                    
                    if ttbl(s,1) > X
                        offset = X - ttbl(s,1);
                        if abs(offset) > abs(xoff)
                            xoff = offset;
                        end
                    end
                    
                    if ttbl(s,2) < 1
                        offset = 1 - ttbl(s,2);
                        if offset > abs(yoff)
                            yoff = offset;
                        end
                    end
                    
                    if ttbl(s,2) > Y
                        offset = Y - ttbl(s,2);
                        if abs(offset) > abs(yoff)
                            yoff = offset;
                        end
                    end
                    
                    if ttbl(s,3) < zmin
                        offset = zmin - ttbl(s,3);
                        if abs(offset) > abs(zoff)
                            zoff = offset;
                        end
                    end
                    
                    if ttbl(s,3) > zmax
                        offset = zmax - ttbl(s,3);
                        if abs(offset) > abs(zoff)
                            zoff = offset;
                        end
                    end
                end
            end
            
            %2 apply offsets
            if (xoff == 0) & (yoff == 0) & (zoff == 0)
                return; 
            end
            
            if xoff ~= 0
                poff = xoff * (-1);
                fprintf("WARNING: X coordinates appear to be shifted %d. Adjusting...\n", poff);
            end
            if yoff ~= 0
                poff = yoff * (-1);
                fprintf("WARNING: Y coordinates appear to be shifted %d. Adjusting...\n", poff);
            end
            if zoff ~= 0
                poff = zoff * (-1); %This yields 0? How??
                fprintf("WARNING: Z coordinates appear to be shifted %d. Adjusting...\n", poff);
            end
            
            for t = 1:T
                ttbl = coord_table{t,1};
                S = size(ttbl,1);
                for s = 1:S
                    ttbl(s,1) = ttbl(s,1) + xoff;
                    ttbl(s,2) = ttbl(s,2) + yoff;
                    ttbl(s,3) = ttbl(s,3) + zoff;
                end
                coord_table{t,1} = ttbl;
            end
            
            save([bf_stem '_coordTable'], 'coord_table');
        end
        
        % ========================== Plots ==========================
        
        function res_struct = generateEmptyPlotable()
            res_struct = struct('x', []);
            res_struct.spot_plot = [];
            res_struct.fscores = []; %If 1 col, as-is. If 3 col, avg, bot, top.
            res_struct.thresh_val = 0;
            res_struct.thresh_info = [];
            res_struct.thresh_val_alt = 0;
            res_struct.thresh_info_alt = [];
        end
        
        function res_struct = generateSensPrecPlotInfo()
            res_struct = struct('sensitivity', []);
            res_struct.precision = [];
            res_struct.fscores = [];
            res_struct.color_by_fscore = true;
            res_struct.color_groups = 20;
            res_struct.color_lo = [0.05 0.05 0.05];
            res_struct.color_hi = [0.95 0.95 0.95]; %Obviously full white doesn't show up
            res_struct.color_mid = [0.5 0.5 0.5];
            res_struct.thresh_idx = 0;
            res_struct.marker = 'x';
            res_struct.plot_name = 'MyPlot';
        end
        
        function drawFScoreSwash(my_plotable, line_style, line_color)
            if isempty(my_plotable.fscores); return; end
            fcols = size(my_plotable.fscores, 2);
            
            plot(my_plotable.x(:,1),my_plotable.fscores(:,1),line_style,'LineWidth',2,'Color',line_color);
            
            if fcols == 3
                plot(my_plotable.x(:,1),my_plotable.fscores(:,2),line_style,'LineWidth',1,'Color',line_color);
                plot(my_plotable.x(:,1),my_plotable.fscores(:,3),line_style,'LineWidth',1,'Color',line_color);
            end
        end
        
        function fig_handle = drawSpotPlot(plotable_hb, plotable_bf, plotable_rsbf)
            figno = round(rand() * 10000);
            legend_names = cell(1, 2);
            color1 = [0.071 0.094 0.529]; %#121887 - Blue
            color2 = [0.529 0.051 0.051]; %#870d0d - Red
            color3 = [0.098 0.776 0.710]; %#19c6b5 - Cyan(ish)
            color1_light = [0.627 0.647 0.969]; %#a0a5f7 - Blue
            color2_light = [0.969 0.588 0.588]; %#f79696 - Red
            color3_light = [0.671 0.827 0.812]; %#19c6b5 - Cyan(ish)
            color4_light = [0.910 0.588 0.969]; %#e896f7 - Purple
            color5_light = [0.706 0.969 0.776]; %#b4f7c6 - Green(?)

            fig_handle = figure(figno);
            clf;
            plot(plotable_bf.x(:,1),log10(double(plotable_bf.spot_plot(:,1))),'LineWidth',2,'Color',color1);
            hold on;
            if ~isempty(plotable_rsbf)
                plot(plotable_rsbf.x(:,1),log10(double(plotable_rsbf.spot_plot(:,1))),':','LineWidth',2,'Color',color3);
            end
            plot(plotable_hb.x(:,1),log10(double(plotable_hb.spot_plot(:,1))),'-.','LineWidth',2,'Color',color2);
            
            if plotable_bf.thresh_val > 0
                xline(plotable_bf.thresh_val, '--', 'Threshold - BigFISH', 'Color', color1_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end
            if plotable_hb.thresh_val > 0
                xline(plotable_hb.thresh_val, '--', 'Threshold - Homebrew', 'Color', color2_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end
            if plotable_bf.thresh_val_alt > 0
                xline(plotable_bf.thresh_val_alt, '--', 'Threshold - Homebrew on BigFISH Callset', 'Color', color4_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end
            legend_names{1,1} = 'BigFISH';
            
            if ~isempty(plotable_rsbf)
                if plotable_rsbf.thresh_val > 0
                    xline(plotable_rsbf.thresh_val, '--', 'Threshold - BigFISH (Rescaled)', 'Color', color3_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
                end
                if plotable_rsbf.thresh_val_alt > 0
                    xline(plotable_rsbf.thresh_val_alt, '--', 'Threshold - Homebrew on BigFISH (Rescaled)', 'Color', color5_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
                end
                legend_names{1,2} = 'BigFISH (Rescaled)';
                legend_names{1,3} = 'Homebrew';
            else
                legend_names{1,2} = 'Homebrew';
            end
            legend(legend_names);
            xlabel('Threshold');
            ylabel('log10(# Spots Detected)');
        end
        
        function fig_handle = drawFScorePlot(plotable_hb, plotable_bf, plotable_rsbf)
            figno = round(rand() * 10000);
            legend_names = cell(1, 2);
            color1 = [0.071 0.094 0.529]; %#121887 - Blue
            color2 = [0.529 0.051 0.051]; %#870d0d - Red
            color3 = [0.098 0.776 0.710]; %#19c6b5 - Cyan(ish)
            color1_light = [0.627 0.647 0.969]; %#a0a5f7 - Blue
            color2_light = [0.969 0.588 0.588]; %#f79696 - Red
            color3_light = [0.671 0.827 0.812]; %#19c6b5 - Cyan(ish)
            color4_light = [0.910 0.588 0.969]; %#e896f7 - Purple
            color5_light = [0.706 0.969 0.776]; %#b4f7c6 - Green(?)
            
            fig_handle = figure(figno);
            clf;
            BigfishCompare.drawFScoreSwash(plotable_bf, '-', color1);
            hold on;
            if ~isempty(plotable_rsbf)
                BigfishCompare.drawFScoreSwash(plotable_rsbf, ':', color3);
            end
            BigfishCompare.drawFScoreSwash(plotable_hb, '-.', color2);
            
            ylim([0.0 1.0]);
            if plotable_bf.thresh_val > 0
                xline(plotable_bf.thresh_val, '--', 'Threshold - BigFISH', 'Color', color1_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end
            if plotable_hb.thresh_val > 0
                xline(plotable_hb.thresh_val, '--', 'Threshold - Homebrew', 'Color', color2_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end
            if plotable_bf.thresh_val_alt > 0
                xline(plotable_bf.thresh_val_alt, '--', 'Threshold - Homebrew on BigFISH Callset', 'Color', color4_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            end
            legend_names{1,1} = 'BigFISH';
            
            if ~isempty(plotable_rsbf)
                if plotable_rsbf.thresh_val > 0
                    xline(plotable_rsbf.thresh_val, '--', 'Threshold - BigFISH (Rescaled)', 'Color', color3_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
                end
                if plotable_rsbf.thresh_val_alt > 0
                    xline(plotable_rsbf.thresh_val_alt, '--', 'Threshold - Homebrew on BigFISH (Rescaled)', 'Color', color5_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
                end
                legend_names{1,2} = 'BigFISH (Rescaled)';
                legend_names{1,3} = 'Homebrew';
            else
                legend_names{1,2} = 'Homebrew';
            end
            legend(legend_names);
            xlabel('Threshold');
            ylabel('FScore');
        end
        
        function sensPrecHighlightTh(plot_info, switch_xy)
            if plot_info.thresh_idx < 1; return; end
            thmax = size(plot_info.sensitivity,1);
            if plot_info.thresh_idx > thmax; return; end
            
            if switch_xy
                y = plot_info.sensitivity(plot_info.thresh_idx);
                x = plot_info.precision(plot_info.thresh_idx);
            else
                x = plot_info.sensitivity(plot_info.thresh_idx);
                y = plot_info.precision(plot_info.thresh_idx);
            end
            
            %fscore_text = ['fscore = ' num2str(plot_info.fscores(plot_info.thresh_idx))];
            plot(x, y, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', plot_info.color_lo, 'MarkerSize', 15);
            %annotation('textarrow', [(x-0.1) x], [(y-0.1) y],'String',{plot_info.plot_name, fscore_text});
        end
        
        function drawSingleSensPrecToCurrentPlot(plot_info, switch_xy)
            if (plot_info.color_by_fscore)
                cgroups = plot_info.color_groups;
                divs = cgroups + 1;
                mid = divs/2;
                last_bound = 0.0;
                for g = 1:cgroups
                    %Derive color
                    plot_color = NaN(1,3);
                    if g < mid
                        midweight = g/mid;
                        loweight = 1.0 - midweight;
                        plot_color(1,:) = (plot_info.color_mid .* midweight) + (plot_info.color_lo .* loweight);
                    elseif g == mid
                        plot_color(1,:) = plot_info.color_mid;
                    else
                        hiweight = (g - mid)/mid;
                        midweight = 1.0 - hiweight;
                        plot_color(1,:) = (plot_info.color_mid .* midweight) + (plot_info.color_hi .* hiweight);
                    end
                    
                    %Subset...
                    new_bound = g/cgroups;
                    if g == cgroups
                        plot_idxs = find(and((plot_info.fscores >= last_bound),(plot_info.fscores <= new_bound)));
                    else
                        plot_idxs = find(and((plot_info.fscores >= last_bound),(plot_info.fscores < new_bound)));
                    end
                    last_bound = new_bound;
                    if isempty(plot_idxs); continue; end
                    
                    %Add each group to plot.
                    if ~switch_xy
                        x = plot_info.sensitivity(plot_idxs);
                        y = plot_info.precision(plot_idxs);
                    else
                        y = plot_info.sensitivity(plot_idxs);
                        x = plot_info.precision(plot_idxs);
                    end
                    plot(x, y, 'LineStyle', 'none', 'Marker', plot_info.marker, 'MarkerEdgeColor', plot_color);
                    hold on;
                end
            else
                if ~switch_xy
                    x = plot_info.sensitivity;
                    y = plot_info.precision;
                else
                    y = plot_info.sensitivity;
                    x = plot_info.precision;
                end
                plot(x, y, 'LineStyle', 'none', 'Marker', plot_info.marker, 'MarkerEdgeColor', plot_info.color_mid);
            end
        end
        
        function fig_handle = drawSensPrecPlot(plot_info_hb, plot_info_bf, plot_info_rsbf, invert_x, invert_y, switch_xy)
            if (nargin < 4); invert_x = false; end
            if (nargin < 5); invert_y = false; end
            if (nargin < 6); switch_xy = false; end
            
            figno = round(rand() * 10000);
            %legend_names = cell(1, 2);
            
            fig_handle = figure(figno);
            clf;
            
            %ax_old = gca;
            ax = axes();
            ax.XLim = [0.0 1.0];
            ax.YLim = [0.0 1.0];
            if invert_x
            	ax.XDir = 'reverse';
            end
            if invert_y
            	ax.YDir = 'reverse';
            end
            
            hold on;
            if ~isempty(plot_info_hb)
                BigfishCompare.drawSingleSensPrecToCurrentPlot(plot_info_hb, switch_xy);
            end
            if ~isempty(plot_info_bf)
                BigfishCompare.drawSingleSensPrecToCurrentPlot(plot_info_bf, switch_xy);
            end
            if ~isempty(plot_info_rsbf)
                BigfishCompare.drawSingleSensPrecToCurrentPlot(plot_info_rsbf, switch_xy);
            end
            
            if ~isempty(plot_info_hb)
                BigfishCompare.sensPrecHighlightTh(plot_info_hb, switch_xy);
            end
            if ~isempty(plot_info_bf)
                BigfishCompare.sensPrecHighlightTh(plot_info_bf, switch_xy);
            end
            if ~isempty(plot_info_rsbf)
                BigfishCompare.sensPrecHighlightTh(plot_info_rsbf, switch_xy);
            end
            
            if switch_xy
                xlabel('Precision');
                ylabel('Sensitivity');
            else
                xlabel('Sensitivity');
                ylabel('Precision');
            end
            axes(ax);
        end
        
        % ========================== Images ==========================
        
        function [xy_match_rows, xyz_match_rows] = coordinateInTable(coord_table, x, y, z)
            xyz_match_rows = [];
            
            x_match = coord_table(:,1) == x;
            y_match = coord_table(:,2) == y;
            xy_match = x_match & y_match;
            
            xy_match_rows = find(xy_match);
            if ~isempty(xy_match_rows)
                xyz_match = coord_table(xy_match_rows,3) == z;
                xyz_match_rows = find(xyz_match);
            end
        end
        
        function fig_handle = drawSpotsOnImage(img_projection, coordset_1, coordset_2, lmin, lmax, figno)
            
            if nargin < 6
                figno = round(rand() * 10000);
            end
            
            colorA = [1.000, 0.000, 0.000];
            colorB = [0.000, 0.000, 1.000];
            colorC = [1.000, 0.000, 1.000]; %Overlap
            colorD = [0.000, 1.000, 1.000]; %XYOverlap
            
            %Split tables between 1 only, 2 only, and overlap
            size1 = size(coordset_1,1);
            size2 = size(coordset_2,1);
            alloc = size1;
            if size2 > alloc; alloc = size2; end
            
            tbl_b = coordset_2;
            tbl_a = zeros(alloc,3); ai = 1;
            tbl_c = zeros(alloc,3); ci = 1;
            tbl_d = zeros(alloc,3); di = 1;
            
            for i = 1:size1
                if isempty(tbl_b)
                    tbl_a(ai,:) = coordset_1(i,:);
                    ai = ai+1;
                    continue;
                end
                
                x =  coordset_1(i,1);
                y =  coordset_1(i,2);
                z =  coordset_1(i,3);
                [xy_match_rows, xyz_match_rows] = BigfishCompare.coordinateInTable(tbl_b, x, y, z);
                
                if ~isempty(xy_match_rows)
                    if ~isempty(xyz_match_rows)
                        tbl_c(ci,:) = coordset_1(i,:);
                        ci = ci+1;
                    else
                        tbl_d(di,:) = coordset_1(i,:);
                        di = di+1;
                    end
                    %Remove these from coordset2
                    tbl_b(xy_match_rows,:) = [];
                else
                    tbl_a(ai,:) = coordset_1(i,:);
                    ai = ai+1;
                end
            end
            
            if ai > 1
                tbl_a = tbl_a(1:ai-1,:);
            else
                tbl_a = [];
            end
            
            if ci > 1
                tbl_c = tbl_c(1:ci-1,:);
            else
                tbl_c = [];
            end
            
            if di > 1
                tbl_d = tbl_d(1:di-1,:);
            else
                tbl_d = [];
            end
            
            fig_handle = figure(figno);
            clf;
            if nargin >= 5
                imshow(img_projection,[lmin lmax]);
            else
                imshow(img_projection,[]);
            end
            hold on;
            
            if ~isempty(tbl_a)
                plot(tbl_a(:,1), tbl_a(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',colorA,'markersize',10);
            end
            if ~isempty(tbl_b)
                plot(tbl_b(:,1), tbl_b(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',colorB,'markersize',10);
            end
            if ~isempty(tbl_c)
                plot(tbl_c(:,1), tbl_c(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',colorC,'markersize',10);
            end
            if ~isempty(tbl_d)
                plot(tbl_d(:,1), tbl_d(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',colorD,'markersize',10);
            end
            
            impixelinfo;
        end
        
        % ========================== Overall ==========================
        
        function [success_bool, fighandles] = doBigfishCompare(hb_stem, bf_dir, bf_imgname, leave_fig, overwrite)
            success_bool = true;
            fighandles.spotsfig = [];
            fighandles.fscorefig = [];
            fighandles.bfhbimg = [];
            fighandles.bfhbbfimg = [];
            fighandles.sensprecfig = [];
             
            if nargin < 4
                leave_fig = false;
            end
            if nargin < 5
                overwrite = false;
            end
            
            %Prepare plotables
            plotable_hb = BigfishCompare.generateEmptyPlotable();
            plotable_bf = BigfishCompare.generateEmptyPlotable();
            plotable_rs = []; %Only gen if rs run found.
            
            plot_info_hb = BigfishCompare.generateSensPrecPlotInfo();
            plot_info_bf = BigfishCompare.generateSensPrecPlotInfo();
            plot_info_rs = []; %Only gen if rs run found.

            %Look for and read summary file.
            summary_file_path = [bf_dir filesep 'summary.txt'];
            if ~isfile(summary_file_path)
                success_bool = false;
                return;
            end
            [zmin, zmax, bfthresh] = BigfishCompare.readSummaryTxt(summary_file_path);
            fprintf("Z Range: %d - %d\n", zmin, zmax);
            fprintf("BF Threshold: %d\n", bfthresh);
            plotable_bf.thresh_val = bfthresh;
            
            bf_stem = [bf_dir filesep bf_imgname];
            
            %Check if need to import
            bf_ct_path = [bf_stem '_coordTable.mat'];
            if overwrite | ~isfile(bf_ct_path)
                fprintf("Importing Big-FISH output...\n");
                [coord_table, spot_table] = BigfishCompare.importBigFishCsvs(bf_dir, bf_stem, zmin);
                if isempty(coord_table); success_bool = false; return; end
            else
                %Load from MAT files
                load([bf_stem '_coordTable'], 'coord_table');
                load([bf_stem '_spotTable'], 'spot_table');
            end
            plotable_bf.x = spot_table(:,1);
            plotable_bf.spot_plot = spot_table(:,2);
            
            %Load HB Reference Run
            fprintf("Loading Reference...\n");
            ref_spots_run = RNASpotsRun.loadFrom(hb_stem);
            [~, ref_spot_table] = ref_spots_run.loadSpotsTable();
            plotable_hb.x = ref_spot_table(:,1);
            plotable_hb.spot_plot = ref_spot_table(:,2);
            plotable_hb.thresh_info = ref_spots_run.threshold_results;
            plotable_hb.thresh_val = ref_spots_run.threshold_results.threshold;
            
            X = ref_spots_run.idims_sample.x;
            Y = ref_spots_run.idims_sample.y;
            coord_table = BigfishCompare.fixCoordinates(zmin, zmax, X, Y, coord_table, bf_stem);

            %Run HB thresholder on BF callset
            fprintf("Thresholding Big-FISH callset...\n");
            bf_thresh_res = BigfishCompare.updateThresholdResults(spot_table, ref_spots_run);
            plotable_bf.thresh_info_alt = bf_thresh_res;
            plotable_bf.thresh_val_alt = bf_thresh_res.threshold;

            %Snap BF refset
            %TODO Multiple refsets?
            use_fscores = RNA_Threshold_SpotSelector.selectorExists(hb_stem);
            fscores_ref = [];
            fscores_bf = [];
            if use_fscores
                fprintf("Snapping Big-FISH callset to reference...\n");
                annoobj = BigfishCompare.loadSpotSelector(bf_stem, hb_stem, spot_table, coord_table, overwrite);
                annoobj.z_min = zmin;
                annoobj.z_max = zmax;
                annoobj = annoobj.updateFTable();
                annoobj.saveMe();
                plotable_bf.fscores = annoobj.f_scores(:,1);
                
                plot_info_bf.plot_name = 'Big-FISH';
                plot_info_bf.sensitivity = annoobj.f_scores(:,2) ./ (annoobj.f_scores(:,2) + annoobj.f_scores(:,4));
                plot_info_bf.precision = annoobj.f_scores(:,2) ./ (annoobj.f_scores(:,2) + annoobj.f_scores(:,3));
                plot_info_bf.fscores = annoobj.f_scores(:,1);
                plot_info_bf.marker = 'd';
                plot_info_bf.color_lo = [0.00 0.00 0.10];
                plot_info_bf.color_hi = [0.10 0.10 0.90];
                plot_info_bf.color_mid = [0.70 0.70 0.90];
                
                th_idx = find(plotable_bf.x == plotable_bf.thresh_val);
                if ~isempty(th_idx)
                    plot_info_bf.thresh_idx = th_idx;
                end
                
                %Load F tables
                src_anno = RNA_Threshold_SpotSelector.openSelector(hb_stem, true);
                %src_anno.z_min = zmin;
                %src_anno.z_max = zmax;
                src_anno.z_min = ref_spots_run.z_min_apply;
                src_anno.z_max = ref_spots_run.z_max_apply;
                src_anno = src_anno.updateFTable();
                plotable_hb.fscores = src_anno.f_scores(:,1);
                
                plot_info_hb.plot_name = 'Homebrew';
                plot_info_hb.sensitivity = src_anno.f_scores(:,2) ./ (src_anno.f_scores(:,2) + src_anno.f_scores(:,4));
                plot_info_hb.precision = src_anno.f_scores(:,2) ./ (src_anno.f_scores(:,2) + src_anno.f_scores(:,3));
                plot_info_hb.fscores = src_anno.f_scores(:,1);
                plot_info_hb.marker = 'd';
                plot_info_hb.color_lo = [0.10 0.00 0.00];
                plot_info_hb.color_hi = [0.90 0.10 0.10];
                plot_info_hb.color_mid = [0.90 0.70 0.70];
                
                th_idx = find(plotable_hb.x == plotable_hb.thresh_val);
                if ~isempty(th_idx)
                    plot_info_hb.thresh_idx = th_idx;
                end
            end
            
            %Check for RS run, and repeat process if found.
            rs_dir = replace(bf_dir, ['data' filesep 'bigfish'], ['data' filesep 'bigfish' filesep '_rescaled']);
            summary_file_path = [rs_dir filesep 'summary.txt'];
            if isfile(summary_file_path)
                rs_okay = true;
                fprintf("Rescaled run found!\n");
                [zmin_rs, zmax_rs, rsthresh] = BigfishCompare.readSummaryTxt(summary_file_path);
                fprintf("Z Range: %d - %d\n", zmin_rs, zmax_rs);
                fprintf("BF Threshold: %d\n", rsthresh);
                
                plotable_rs = BigfishCompare.generateEmptyPlotable();
                plotable_rs.thresh_val = rsthresh;
                rs_stem = [rs_dir filesep bf_imgname];
                
                rs_ct_path = [rs_stem '_coordTable.mat'];
                if overwrite | ~isfile(rs_ct_path)
                    fprintf("Importing RS Big-FISH output...\n");
                    [coord_table_rs, spot_table_rs] = BigfishCompare.importBigFishCsvs(rs_dir, rs_stem, zmin_rs);
                    if isempty(coord_table) 
                        rs_okay = false;
                    end
                else
                    %Load from MAT files
                    coord_table_temp = coord_table;
                    spot_table_temp = spot_table;
                    load([rs_stem '_coordTable'], 'coord_table');
                    load([rs_stem '_spotTable'], 'spot_table');
                    coord_table_rs = coord_table;
                    spot_table_rs = spot_table;
                    coord_table = coord_table_temp;
                    spot_table = spot_table_temp;
                    if isempty(coord_table_rs) 
                        rs_okay = false;
                    end
                end
                
                if rs_okay
                    plotable_rs.x = spot_table_rs(:,1);
                    plotable_rs.spot_plot = spot_table_rs(:,2);
                    
                    coord_table_rs = BigfishCompare.fixCoordinates(zmin_rs, zmax_rs, X, Y, coord_table_rs, rs_stem);
                
                    fprintf("Thresholding RS Big-FISH callset...\n");
                    rs_thresh_res = BigfishCompare.updateThresholdResults(spot_table_rs, ref_spots_run);
                    plotable_rs.thresh_info_alt = rs_thresh_res;
                    plotable_rs.thresh_val_alt = rs_thresh_res.threshold;
                    
                    if use_fscores
                        %TODO Multiple refsets?
                        fprintf("Snapping RS Big-FISH callset to reference...\n");
                        annoobj = BigfishCompare.loadSpotSelector(rs_stem, hb_stem, spot_table_rs, coord_table_rs, overwrite);
                        annoobj.z_min = zmin_rs;
                        annoobj.z_max = zmax_rs;
                        annoobj = annoobj.updateFTable();
                        annoobj.saveMe();
                        plotable_rs.fscores = annoobj.f_scores(:,1);
                        
                        plot_info_rs = BigfishCompare.generateSensPrecPlotInfo();
                        plot_info_rs.plot_name = 'Big-FISH (Rescaled)';
                        plot_info_rs.sensitivity = annoobj.f_scores(:,2) ./ (annoobj.f_scores(:,2) + annoobj.f_scores(:,4));
                        plot_info_rs.precision = annoobj.f_scores(:,2) ./ (annoobj.f_scores(:,2) + annoobj.f_scores(:,3));
                        plot_info_rs.fscores = annoobj.f_scores(:,1);
                        plot_info_rs.marker = 'd';
                        plot_info_rs.color_lo = [0.00 0.10 0.00];
                        plot_info_rs.color_hi = [0.10 0.90 0.10];
                        plot_info_rs.color_mid = [0.70 0.90 0.70];
                        
                        th_idx = find(plotable_rs.x == plotable_rs.thresh_val);
                        if ~isempty(th_idx)
                            plot_info_rs.thresh_idx = th_idx;
                        end
                    end
                else
                    plotable_rs = [];
                end
                
            end
            
            %Print graphs
            plotdir = [bf_dir filesep 'plots'];
            mkdir(plotdir);

            fighandles.spotsfig = BigfishCompare.drawSpotPlot(plotable_hb, plotable_bf, plotable_rs);
            saveas(fighandles.spotsfig, [plotdir filesep 'logspots.png']);
            if ~leave_fig; close(fighandles.spotsfig); fighandles.spotsfig = []; end

            fighandles.fscorefig = [];
            if use_fscores
                fighandles.fscorefig = BigfishCompare.drawFScorePlot(plotable_hb, plotable_bf, plotable_rs);
                saveas(fighandles.fscorefig, [plotdir filesep 'fscoreplot.png']);
                if ~leave_fig; close(fighandles.fscorefig); fighandles.fscorefig = []; end
                
                fighandles.sensprecfig = BigfishCompare.drawSensPrecPlot(plot_info_hb, plot_info_bf, plot_info_rs, true, false, false);
                saveas(fighandles.sensprecfig, [plotdir filesep 'sensitivity_precision.png']);
                if ~leave_fig; close(fighandles.sensprecfig); fighandles.sensprecfig = []; end
            end

            %Render spot images?
            [~, img_channel] = ref_spots_run.loadFilteredImage();
            if ~isempty(img_channel)
                %TODO HB and BF runs have different z min/max now?
                img_channel = img_channel(:,:,zmin:zmax);
                max_proj = double(max(img_channel,[],3));
                clear img_channel;
                Lmin = median(max_proj(:)) - round(0 * std(max_proj(:)));
                Lmax = median(max_proj(:)) + round(10 * std(max_proj(:)));
            
                %BF vs. HB each at chosen thresh
                [~, ref_ttbl] = ref_spots_run.loadThresholdTable();
                ref_tidx = find(ref_ttbl >= ref_spots_run.threshold_results.threshold, 1);
                bfhb_tidx = find(spot_table(:,1) >= bf_thresh_res.threshold, 1);
                bf_tidx = find(spot_table(:,1) >= bfthresh, 1);

                [~, ref_coord_table] = ref_spots_run.loadCoordinateTable();
                bf_coordset = coord_table{bf_tidx,1};
                hb_coordset = ref_coord_table{ref_tidx};
                fighandles.bfhbimg = BigfishCompare.drawSpotsOnImage(max_proj, hb_coordset, bf_coordset, Lmin, Lmax, 615);
                saveas(fighandles.bfhbimg, [plotdir filesep 'bf_v_hb_spots.png']);
                if ~leave_fig; close(fighandles.bfhbimg); fighandles.bfhbimg = []; end
                clear hb_coordset;
                clear ref_coord_table;
            
                %BF vs. BF curve/HB thresholder
                hb_coordset = coord_table{bfhb_tidx};
                fighandles.bfhbbfimg = BigfishCompare.drawSpotsOnImage(max_proj, hb_coordset, bf_coordset, Lmin, Lmax, 301);
                saveas(fighandles.bfhbbfimg, [plotdir filesep 'bf_v_hbbf_spots.png']);
                if ~leave_fig; close(fighandles.bfhbbfimg); fighandles.bfhbbfimg = []; end
            end
        end
        
    end

end