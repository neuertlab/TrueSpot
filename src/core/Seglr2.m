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
            linfit = struct("left", [], "right", [], "break_index", -1, "segrsq", 0.0, "rsqwavg", 0.0);
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
        
    end
end