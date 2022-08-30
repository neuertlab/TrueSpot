%
%%

classdef Seglr2
    
    methods (Static)
        
        function linestruct = genEmptyLineStruct()
            linestruct = struct("slope", 0.0, "yintr", 0.0, "rsquared", 0.0, "segrsq", 0.0, "rsqwavg", 0.0);
        end
        
        function linestruct = genLineStruct(slope, yintr, rsq)
            linestruct = struct("slope", slope, "yintr", yintr, "rsquared", rsq, "segrsq", 0.0, "rsqwavg", 0.0);
        end
        
        function linfit = genEmptySegfitStruct()
            linfit = struct("left", [], "right", [], "break_index", -1, "segrsq", 0.0, "rsqwavg", 0.0, ...
                "break_x", 0, "break_y", 0.0, "input_x0", [], "input_y0", [], "input_bLeft", [], "input_bRight", []);
        end
        
        function fitp = genEmptyFitParamStruct()
            fitp = struct("min", 0.0, "max", 0.0, "start", 0.0);
        end
        
        function linfit = fitToSpotCurve(data, xmin, xmax, verbosity)
            
            linfit = Seglr2.genEmptySegfitStruct();
            min_points = 3;
            data_points = size(data,1);
            
            if nargin < 4
                verbosity = 1;
            end
            scanmin = min_points;
            scanmax = data_points - min_points;
            if nargin > 1
                %Min provided.
                usrmin = find(data(:,1) == xmin, 1);
                if usrmin > scanmin
                    scanmin = usrmin;
                end
            end
            if nargin > 2
                %Max provided.
                usrmax = find(data(:,1) == xmax, 1);
                if usrmax < scanmax
                    scanmax = usrmax;
                end
            end
            
            if verbosity > 1
                fprintf("Calculating scan ranges...\n");
            end
            
            %Find Y knot range
            ystd = std(data(:,2),0,'all','omitnan');
            y0_min = min(data(:,2),[],'omitnan') - ystd;
            y0_max = mean(data(:,2),'all','omitnan');
            y0_st = y0_min + (rand() * (y0_max - y0_min));
            
            %Find left intercept range
            [~, testb] = Seglr2.lineeqs(data(2:data_points,1),data(2:data_points,2),data(1,1), data(1,2));
            b1_min = min(testb,[],'all','omitnan') - ystd;
            b1_max = max(testb,[],'all','omitnan');
            b1_st = b1_min + (rand() * (b1_max - b1_min));
            
            %Find right intercept range
            [~,minidx] = min(data(:,2),[],'omitnan');
            [~, testb] = Seglr2.lineeqs(data(:,1),data(:,2),data(minidx,1), data(minidx,2));
            b2_min = min(testb,[],'all','omitnan');
            b2_max = max(testb,[],'all','omitnan');
            b2_st = b2_min + (rand() * (b2_max - b2_min));
            
            x0_st = scanmin + (rand() * (scanmax - scanmin));
            
            %(X0,Y0,b1,b2,x)
            if verbosity > 1
                fprintf("Performing fit...\n");
            end
            ftobj = fittype('((((x .* (Y0 - b1)) / X0) + b1) .* (x <= X0)) + ((((x .* (Y0 - b2)) / X0) + b2) .* (x > X0))');
            cfobj = fit(data(:,1),data(:,2),ftobj,'StartPoint',[x0_st,y0_st,b1_st,b2_st],'Lower',[scanmin,y0_min,b1_min,b2_min],'Upper',[scanmax,y0_max,b1_max,b2_max]);
            
            if verbosity > 1
                fprintf("Saving fit parameters...\n");
            end
            
            linfit.break_x = cfobj.X0;
            linfit.break_y = cfobj.Y0;
            for i = 1:data_points
                if data(i,1) >= linfit.break_x
                    linfit.break_index = i;
                    break;
                end
            end
            
            linfit.left = Seglr2.genEmptyLineStruct();
            linfit.left.yintr = cfobj.b1;
            linfit.left.slope = (linfit.break_y - linfit.left.yintr)/(linfit.break_x);
            
            linfit.right = Seglr2.genEmptyLineStruct();
            linfit.right.yintr = cfobj.b2;
            linfit.right.slope = (linfit.break_y - linfit.right.yintr)/(linfit.break_x);
            
            %rsq values
            linfit.left.rsquared = Seglr2.rsquared(data, linfit.left.slope, linfit.left.yintr);
            linfit.right.rsquared = Seglr2.rsquared(data, linfit.right.slope, linfit.right.yintr);
            linfit.left.segrsq = Seglr2.rsquared_seg(data, linfit.break_index, linfit.left, linfit.right);
            linfit.right.segrsq = linfit.left.segrsq;
            linfit.segrsq = linfit.left.segrsq;
            data_l = data(1:linfit.break_index,:);
            data_r = data(linfit.break_index:data_points,:);
            linfit.left.rsqwavg = Seglr2.rescoreWeighted(data_l, data_r, linfit.left, linfit.right);
            linfit.right.rsqwavg = linfit.left.rsqwavg;
            linfit.rsqwavg = linfit.left.rsqwavg;
            
            %For record keeping...
            linfit.input_x0 = Seglr2.genEmptyFitParamStruct();
            linfit.input_x0.start = x0_st;
            linfit.input_x0.min = scanmin;
            linfit.input_x0.max = scanmax;
            
            linfit.input_y0 = Seglr2.genEmptyFitParamStruct();
            linfit.input_y0.start = y0_st;
            linfit.input_y0.min = y0_min;
            linfit.input_y0.max = y0_max;
            
            linfit.input_bLeft = Seglr2.genEmptyFitParamStruct();
            linfit.input_bLeft.start = b1_st;
            linfit.input_bLeft.min = b1_min;
            linfit.input_bLeft.max = b1_max;
            
            linfit.input_bRight = Seglr2.genEmptyFitParamStruct();
            linfit.input_bRight.start = b2_st;
            linfit.input_bRight.min = b2_min;
            linfit.input_bRight.max = b2_max;
            
        end
        
        function linfit = fitTo(data, xmin, xmax, iterations, verbosity)
            min_points = 3;
            data_points = size(data,1);
            if data_points < ((min_points * 2) - 1)
                return;
            end
            
            if nargin < 4
                iterations = 3;
            end
            if nargin < 5
                verbosity = 1;
            end
            
            %Determine knot scan range
            scanmin = min_points;
            scanmax = data_points - min_points;
            if nargin > 1
                %Min provided.
                usrmin = find(data(:,1) == xmin, 1);
                if usrmin > scanmin
                    scanmin = usrmin;
                end
            end
            if nargin > 2
                %Max provided.
                usrmax = find(data(:,1) == xmax, 1);
                if usrmax < scanmax
                    scanmax = usrmax;
                end
            end
            
            %Preallocate result array.
            linearr(data_points) = Seglr2.genEmptySegfitStruct();
%             dbg_srsq = NaN(data_points,1);
%             dbg_rsq_l = NaN(data_points,1);
%             dbg_rsq_r = NaN(data_points,1);
%             dbg_wavg = NaN(data_points,1);
            
            for i = scanmin:scanmax
                if verbosity > 0
                    fprintf("Trying breakpoint x = %d...\n", data(i,1));
                end
                
                %Do initial fit to get y range.
                data_l = data(1:i,:);
                data_r = data(i:data_points,:);
                
                rawfit_l = fitlm(data_l(:,1), data_l(:,2));
                rawfit_r = fitlm(data_r(:,1), data_r(:,2));
                
                yrng = NaN(1,3);
                yrng(1,1) = data(i,2);
                yrng(1,2) = predict(rawfit_l, data(i,1));
                yrng(1,3) = predict(rawfit_r, data(i,1));
                
                yscan_min = min(yrng, [], 'all', 'omitnan');
                yscan_max = max(yrng, [], 'all', 'omitnan');
                
                if verbosity > 0
                    fprintf("Scanning for knot anchor between y = %d and %d...\n", yscan_min, yscan_max);
                end
                
                linearr(i).break_index = i;
                [linearr(i).left, linearr(i).right] = Seglr2.tryBreakpoint(data, i, yscan_min, yscan_max, iterations);
                if isempty(linearr(i).left) | isempty(linearr(i).right)
                    continue;
                end
                linearr(i).segrsq = linearr(i).left.segrsq;
                linearr(i).rsqwavg = linearr(i).left.rsqwavg;
                
%                 dbg_srsq(i,1) = Seglr2.rsquared_seg(data, i, linearr(i).left, linearr(i).right);
%                 dbg_rsq_l(i,1) = linearr(i).left.rsquared;
%                 dbg_rsq_r(i,1) = linearr(i).right.rsquared;
%                 dbg_wavg(i,1) = linearr(i).rsqwavg;
            end
            
            %Find segmented rsq maximum.
            best_score = 0.0;
            best_break = -1;
            
            for i = scanmin:scanmax
                if linearr(i).rsqwavg > best_score
                    best_score = linearr(i).rsqwavg;
                    best_break = i;
                end
            end
            
            linfit = [];
            if best_break > 0
                linfit = linearr(best_break);
            end
            
            if ~isempty(linfit)
                linfit.segrsq = Seglr2.rsquared_seg(data, linfit.break_index, linfit.left, linfit.right);
            end
            
        end
        
        function [line_left, line_right] = tryBreakpoint(data, break_idx, ymin, ymax, nitr)
            %Split pieces.
            total_points = size(data,1);
            data_l = data(1:break_idx,:);
            data_r = data(break_idx:total_points,:);
            
            x_anchor = data(break_idx,1);
            if nargin < 3
                ymin = data(break_idx,2) - 1000;
                ymax = data(break_idx,2) + 1000;
            end
            if nargin < 5
                nitr = 3;
            end
            
            %Prepare bookkeeping
            %map_l = containers.Map({0.0},{Seglr2.genEmptyLineStruct()});
            %map_r = containers.Map({0.0},{Seglr2.genEmptyLineStruct()});
            
            %Set up some variables
            testcvr = ymax - ymin;
            testamt = 32;
            testinv = testcvr/(testamt-1);
            
            bestscore = 0.0; %This floors it for speed. Toss breakpoints that don't look like they can gen a good fit.
            bestyval = NaN;
            bestline_l = [];
            bestline_r = [];
            
            %Force lines to intersect at breakx. Try various y values.
            for i = 1:nitr
                testyval = ymin;
                for j = 1:testamt
                    %See if this value has already been tested...
%                     if (testyval ~= 0.0) & (isKey(map_l, testyval))
%                         fit_l = map_l(testyval);
%                     else
%                         fit_l = Seglr2.fitAnchored(data_l, x_anchor, testyval);
%                         map_l(testyval) = fit_l;
%                     end
%                     
%                     if (testyval ~= 0.0) & (isKey(map_r, testyval))
%                         fit_r = map_r(testyval);
%                     else
%                         fit_r = Seglr2.fitAnchored(data_r, x_anchor, testyval);
%                         map_r(testyval) = fit_r;
%                     end

                    %Retests are so rare, this is actually faster.
                    fit_l = Seglr2.fitAnchored(data_l, x_anchor, testyval);
                    fit_r = Seglr2.fitAnchored(data_r, x_anchor, testyval);
                    
                    %Find weighted rsq for whole thing
%                     segrsq = Seglr2.rsquared_seg(data, break_idx, fit_l, fit_r);
%                     fit_l.segrsq = segrsq;
%                     fit_r.segrsq = segrsq;
                    rsqwavg = Seglr2.rescoreWeighted(data_l, data_r, fit_l, fit_r);
                    fit_l.rsqwavg = rsqwavg;
                    fit_r.rsqwavg = rsqwavg;
                    if rsqwavg > bestscore
                        bestscore = rsqwavg;
                        bestyval = testyval;
                        bestline_l = fit_l;
                        bestline_r = fit_r;
                    end
                    
                    testyval = testyval + testinv;
                end
                
                %Prepare to "zoom in" for next iteration
                if isnan(bestyval)
                    %nothing fits
                    break;
                end
                
                testinv = testinv / 2.0;
                scanamt = testinv * ((testamt-1)/2.0);
                ymin = bestyval - scanamt;
            end
            
            %Prepare return value
            line_left  = bestline_l;
            line_right = bestline_r;
        end
        
        function score = rescoreWeighted(data_l, data_r, fit_l, fit_r)
            n_left = size(data_l,1);
            n_right = size(data_r,1);
            n_total = n_left + n_right - 1;
            weighted_left = n_left * fit_l.rsquared;
            weighted_right = n_right * fit_r.rsquared;
            score = (weighted_left + weighted_right)/n_total;
        end
        
        function fit_line = fitAnchored(data, x_anchor, y_anchor, nitr)
            xokay = find(data(:,1) ~= x_anchor);
            
            %Gen a line between anchor and every data point
            [mtest, btest] = Seglr2.lineeqs(x_anchor, y_anchor, transpose(data(xokay,1)), transpose(data(xokay,2)));
            rsqtest = Seglr2.rsquared(data, mtest, btest);
            test_ct = size(rsqtest,2);
            
            if nargin < 4
                nitr = 3;
            end
            
            ivaldiv = 32;
            
            for i = 1:nitr
                %Find max rsq to test.
                [bestrsq, maxidx] = max(rsqtest);
                bestyintr = btest(maxidx);
                bestslope = mtest(maxidx);
            
                %Sort y intercepts
                btest_sorted = sort(btest);
                best_index = find(btest_sorted == bestyintr,1);
                if best_index > 1
                    b_min = btest_sorted(best_index-1);
                else
                    b_min = NaN;
                end
                
                if best_index < test_ct
                    b_max = btest_sorted(best_index+1);
                else
                    b_max = NaN;
                end
                
                if isnan(b_min) & isnan(b_max)
                    fit_line = Seglr2.genLineStruct(bestslope, bestyintr, bestrsq);
                    return;
                end
                
                if isnan(b_min)
                    b_min = bestyintr - (b_max - bestyintr);
                end
                
                if isnan(b_max)
                    b_max = bestyintr + (bestyintr - b_min);
                end
                
                %Vectorize
                intrv = (b_max - b_min) / (ivaldiv-1);
                if intrv == 0
                    fit_line = Seglr2.genLineStruct(bestslope, bestyintr, bestrsq);
                    return;
                end
                
                test_yintr = [b_min:intrv:b_max];
                zerox = zeros(1,size(test_yintr,2));
                
                %Generate new set of lines and test.
                [mtest, btest] = Seglr2.lineeqs(x_anchor, y_anchor, zerox, test_yintr);
                rsqtest = Seglr2.rsquared(data, mtest, btest);
                test_ct = size(rsqtest,2);
            end
            
            %Adjust best.
            [bestrsq, maxidx] = max(rsqtest);
            bestyintr = btest(maxidx);
            bestslope = mtest(maxidx);
            
            fit_line = Seglr2.genLineStruct(bestslope, bestyintr, bestrsq);
        end
        
        function [slope, yinter] = lineeqs(xa, ya, xtest, ytest)
            slope = (ytest - ya) ./ (xtest - xa);
            yinter = ya - (xa.*slope);
        end
        
        function rsq = rsquared(data, slope, yinter)
            %Same for all lines
            yavg = mean(data(:,2), 'all', 'omitnan');
            dvec = data(:,2) - yavg;
            
            evec = data(:,1) * slope;
            evec = yinter + evec;
            nvec = data(:,2) - evec;
            
            nvec = nvec.*nvec;
            dvec = dvec.*dvec;
            
            nsum = sum(nvec,1,'omitnan');
            dsum = sum(dvec,1,'omitnan');
            
            rsq = 1.0 - (nsum./dsum);
        end
        
        function rsq = rsquared_seg(data, breakx_idx, line_left, line_right)
            ptct = size(data,1);
            yavg = mean(data(:,2), 'all', 'omitnan');
            dvec = data(:,2) - yavg;
            
            evec = NaN(ptct,1);
            evec(1:breakx_idx,1) = line_left.yintr + (line_left.slope .* data(1:breakx_idx,1));
            evec(breakx_idx+1:ptct,1) = line_right.yintr + (line_right.slope .* data(breakx_idx+1:ptct,1));
            nvec = data(:,2) - evec;
            
            nvec = nvec.*nvec;
            dvec = dvec.*dvec;
            
            nsum = sum(nvec,1,'omitnan');
            dsum = sum(dvec,1,'omitnan');
            
            rsq = 1.0 - (nsum./dsum);
        end
        
        function fig_handle = renderFit(data, linfit, show_scan_ranges, figno)
            if nargin < 4
                figno = round(rand() * 10000);
            end
            
            color_grey = [0.290, 0.290, 0.290];
            color_red = [0.929, 0.229, 0.229];
            color_blue = [0.229, 0.229, 0.929];
            
            fig_handle = figure(figno);
            clf;
            ax = axes;
            plot(data(:,1), data(:,2),'LineWidth',2,'Color',color_grey);
            hold on;
            
            line_x = data(1:linfit.break_index,1);
            line_y = linfit.left.yintr + (line_x.*linfit.left.slope);
            plot(line_x, line_y,'LineWidth',1,'Color',color_red,'LineStyle','--');
            
            line_x = data(linfit.break_index:size(data,1),1);
            line_y = linfit.right.yintr + (line_x.*linfit.right.slope);
            plot(line_x, line_y,'LineWidth',1,'Color',color_blue,'LineStyle','--');
            
            if show_scan_ranges
                %TODO
            end
            
        end
        
    end
end