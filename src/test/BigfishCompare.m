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
                    tcoords(:,3) = raw_coord_table(:,1) + 1 + z_offset;
        
                    coord_table{i} = tcoords;
                else
                    break;
                end
                i = i+1;
            end
            
            %Save to mat files
            save([output_stem '_coordTable'], 'coord_table');
            save([output_stem '_spotTable'], 'spot_table');
        end

        function annoobj = loadSpotSelector(stem, ref_stem, bf_spot_table, bf_coord_table, overwrite)
           
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
                        annoobj = annoobj.refSnapToAutoSpots();
                        annoobj = annoobj.updateFTable();
                        annoobj.ref_last_modified = datetime;
                    end
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
                annoobj = annoobj.refSnapToAutoSpots();
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
            zmin = str2num(strtrim(sspl{1}));
            zmax = str2num(strtrim(sspl{2}));

            line = fgetl(fhandle);
            sspl = split(line, ':');
            bfthresh = round(str2num(strtrim(sspl{2})));
            fclose(fhandle);
        end
        
        % ========================== Plots ==========================
        
        function fig_handle = drawSpotPlot(ref_spots_table, bf_spot_table, bf_thresh, hbbf_thresh_res, hb_thresh_res)
            figno = round(rand() * 10000);
            legend_names = cell(1, 2);
            color1 = [0.071 0.094 0.529]; %#121887 - Blue
            color2 = [0.529 0.051 0.051]; %#870d0d - Red
            color1_light = [0.627 0.647 0.969]; %#a0a5f7
            color2_light = [0.969 0.588 0.588]; %#f79696
            color3_light = [0.910 0.588 0.969]; %#e896f7

            fig_handle = figure(figno);
            clf;
            plot(bf_spot_table(:,1),log10(bf_spot_table(:,2)),'LineWidth',2,'Color',color1);
            hold on;
            plot(ref_spots_table(:,1),log10(ref_spots_table(:,2)),'-.','LineWidth',2,'Color',color2);
            xline(bf_thresh, '--', 'Threshold - BIG-FISH', 'Color', color1_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            xline(hb_thresh_res.threshold, '--', 'Threshold - Homebrew', 'Color', color2_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            xline(hbbf_thresh_res.threshold, '--', 'Threshold - Homebrew on BIG-FISH Callset', 'Color', color3_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            legend_names{1,1} = 'BIG-FISH';
            legend_names{1,2} = 'Homebrew';
            legend(legend_names);
            xlabel('Threshold');
            ylabel('log10(# Spots Detected)');
        end
        
        function fig_handle = drawFScorePlot(ref_fscores, bf_fscores, bf_thresh, hbbf_thresh_res, hb_thresh_res)
            figno = round(rand() * 10000);
            legend_names = cell(1, 2);
            color1 = [0.071 0.094 0.529]; %#121887 - Blue
            color2 = [0.529 0.051 0.051]; %#870d0d - Red
            color1_light = [0.627 0.647 0.969]; %#a0a5f7
            color2_light = [0.969 0.588 0.588]; %#f79696
            color3_light = [0.910 0.588 0.969]; %#e896f7
            
            fig_handle = figure(figno);
            clf;
            plot(bf_fscores(:,1),bf_fscores(:,2),'LineWidth',2,'Color',color1);
            hold on;
            plot(ref_fscores(:,1),ref_fscores(:,2),'-.','LineWidth',2,'Color',color2);
            ylim([0.0 1.0]);
            xline(bf_thresh, '--', 'Threshold - BIG-FISH', 'Color', color1_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            xline(hb_thresh_res.threshold, '--', 'Threshold - Homebrew', 'Color', color2_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            xline(hbbf_thresh_res.threshold, '--', 'Threshold - Homebrew on BIG-FISH Callset', 'Color', color3_light,'LineWidth',2,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
            legend_names{1,1} = 'BIG-FISH';
            legend_names{1,2} = 'Homebrew';
            legend(legend_names);
            xlabel('Threshold');
            ylabel('FScore');
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
            
            if nargin < 4
                leave_fig = false;
            end
            if nargin < 5
                overwrite = false;
            end

            %Look for and read summary file.
            summary_file_path = [bf_dir filesep 'summary.txt'];
            if ~isfile(summary_file_path)
                success_bool = false;
                return;
            end
            [zmin, zmax, bfthresh] = BigfishCompare.readSummaryTxt(summary_file_path);
            fprintf("Z Range: %d - %d\n", zmin, zmax);
            fprintf("BF Threshold: %d\n", bfthresh);
            
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

            %Load HB Reference Run
            fprintf("Loading Reference...\n");
            ref_spots_run = RNASpotsRun.loadFrom(hb_stem);
            [~, ref_spot_table] = ref_spots_run.loadSpotsTable();

            %Run HB thresholder on BF callset
            fprintf("Thresholding Big-FISH callset...\n");
            bf_thresh_res = BigfishCompare.updateThresholdResults(spot_table, ref_spots_run);

            %Snap BF refset
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
                T = size(annoobj.f_scores, 1);
                fscores_bf = NaN(T,2);
                fscores_bf(:,1) = spot_table(:,1);
                fscores_bf(:,2) = annoobj.f_scores(:,1);

                %Load F tables
                src_anno = RNA_Threshold_SpotSelector.openSelector(hb_stem, true);
                src_anno.z_min = zmin;
                src_anno.z_max = zmax;
                src_anno = src_anno.updateFTable();
                T = size(src_anno.f_scores, 1);
                fscores_ref = NaN(T,2);
                fscores_ref(:,1) = ref_spot_table(:,1);
                fscores_ref(:,2) = src_anno.f_scores(:,1);
            end

            %Print graphs
            plotdir = [bf_dir filesep 'plots'];
            mkdir(plotdir);

            fighandles.spotsfig = BigfishCompare.drawSpotPlot(ref_spot_table, spot_table, bfthresh, bf_thresh_res, ref_spots_run.threshold_results);
            saveas(fighandles.spotsfig, [plotdir filesep 'logspots.png']);
            if ~leave_fig; close(fighandles.spotsfig); fighandles.spotsfig = []; end

            fighandles.fscorefig = [];
            if use_fscores
                fighandles.fscorefig = BigfishCompare.drawFScorePlot(fscores_ref, fscores_bf, bfthresh, bf_thresh_res, ref_spots_run.threshold_results);
                saveas(fighandles.fscorefig, [plotdir filesep 'fscoreplot.png']);
                if ~leave_fig; close(fighandles.fscorefig); fighandles.fscorefig = []; end
            end

            %Render spot images?
            [~, img_channel] = ref_spots_run.loadFilteredImage();
            if ~isempty(img_channel)
                img_channel = img_channel(:,:,zmin:zmax);
                max_proj = double(max(img_channel,[],3));
                clear img_channel;
                Lmin = median(max_proj(:)) - round(0 * std(max_proj(:)));
                Lmax = median(max_proj(:)) + round(10 * std(max_proj(:)));
            
                %BF vs. HB each at chosen thresh
                [~, ref_coord_table] = ref_spots_run.loadCoordinateTable();
                bf_coordset = coord_table{bfthresh,1};
                hb_coordset = ref_coord_table{ref_spots_run.threshold_results.threshold};
                fighandles.bfhbimg = BigfishCompare.drawSpotsOnImage(max_proj, hb_coordset, bf_coordset, Lmin, Lmax, 615);
                saveas(fighandles.bfhbimg, [plotdir filesep 'bf_v_hb_spots.png']);
                if ~leave_fig; close(fighandles.bfhbimg); fighandles.bfhbimg = []; end
                clear hb_coordset;
                clear ref_coord_table;
            
                %BF vs. BF curve/HB thresholder
                hb_coordset = coord_table{bf_thresh_res.threshold};
                fighandles.bfhbbfimg = BigfishCompare.drawSpotsOnImage(max_proj, hb_coordset, bf_coordset, Lmin, Lmax, 301);
                saveas(fighandles.bfhbbfimg, [plotdir filesep 'bf_v_hbbf_spots.png']);
                if ~leave_fig; close(fighandles.bfhbbfimg); fighandles.bfhbbfimg = []; end
            end
        end
        
    end

end