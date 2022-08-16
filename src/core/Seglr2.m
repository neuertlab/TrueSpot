%
%%

classdef Seglr2
    
    methods (Static)
        
        function linestruct = genEmptyLineStruct()
            linestruct = struct("slope", 0.0, "yintr", 0.0, "rsquared", 0.0, "segrsq", 0.0);
        end
        
        function linestruct = genLineStruct(slope, yintr, rsq)
            linestruct = struct("slope", slope, "yintr", yintr, "rsquared", rsq, "segrsq", 0.0);
        end
        
        function linfit = genEmptySegfitStruct()
            linfit = struct("left", [], "right", [], "break_index", -1, "segrsq", 0.0);
        end
        
        function linfit = fitTo(data, itr_per_break, verbosity)
            min_points = 3;
            data_points = size(data,1);
            if data_points < ((min_points * 2) - 1)
                return;
            end
            
            if nargin < 2
                itr_per_break = 10;
            end
            if nargin < 3
                verbosity = 1;
            end
            
            %Preallocate result array.
            linearr(data_points) = Seglr2.genEmptySegfitStruct();
            dbg_srsq = NaN(data_points,1);
            
            for i = min_points:(data_points - min_points)
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
                [linearr(i).left, linearr(i).right] = Seglr2.tryBreakpoint(data, i, yscan_min, yscan_max, itr_per_break);
                linearr(i).segrsq = linearr(i).left.segrsq;
                dbg_srsq(i,1) = linearr(i).segrsq;
            end
            
            %Find segmented rsq maximum.
            best_srsq = 0.0;
            best_break = -1;
            
            for i = min_points:(data_points - min_points)
                if linearr(i).segrsq > best_srsq
                    best_srsq = linearr(i).segrsq;
                    best_break = i;
                end
            end
            
            linfit = [];
            if best_break > 0
                linfit = linearr(best_break);
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
                nitr = 10;
            end
            
            %Prepare bookkeeping
            map_l = containers.Map({0.0},{Seglr2.genEmptyLineStruct()});
            map_r = containers.Map({0.0},{Seglr2.genEmptyLineStruct()});
            
            %Set up some variables
            testcvr = ymax - ymin;
            testamt = 64;
            testinv = testcvr/(testamt-1);
            
            bestsegrsq = 0.000;
            bestyval = NaN;
            bestline_l = [];
            bestline_r = [];
            
            %Force lines to intersect at breakx. Try various y values.
            for i = 1:nitr
                testyval = ymin;
                for j = 1:testamt
                    %See if this value has already been tested...
                    if (testyval ~= 0.0) & (isKey(map_l, testyval))
                        fit_l = map_l(testyval);
                    else
                        fit_l = Seglr2.fitAnchored(data_l, x_anchor, testyval);
                        map_l(testyval) = fit_l;
                    end
                    
                    if (testyval ~= 0.0) & (isKey(map_r, testyval))
                        fit_r = map_r(testyval);
                    else
                        fit_r = Seglr2.fitAnchored(data_r, x_anchor, testyval);
                        map_r(testyval) = fit_r;
                    end
                    
                    %Find overall seg rsq
                    segrsq = Seglr2.rsquared_seg(data, break_idx, fit_l, fit_r);
                    fit_l.segrsq = segrsq;
                    fit_r.segrsq = segrsq;
                    if segrsq > bestsegrsq
                        bestsegrsq = segrsq;
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
        
        function index = optRsqValTwoPiece(rsq_left, rsq_right)
            avg = (rsq_left + rsq_right) ./ 2;
            diffsq = (rsq_left - rsq_right).^2;
            weighted = avg - diffsq;
            [~,index] = max(weighted,'all','omitnan');
        end
        
        function fit_line = fitAnchored(data, x_anchor, y_anchor, nitr)
            xokay = find(data(:,1) ~= x_anchor);
            
            %Gen a line between anchor and every data point
            [mtest, btest] = Seglr2.lineeqs(x_anchor, y_anchor, transpose(data(xokay,1)), transpose(data(xokay,2)));
            rsqtest = Seglr2.rsquared(data, mtest, btest);
            test_ct = size(rsqtest,2);
            
            if nargin < 4
                nitr = 10;
            end
            
            ivaldiv = 64;
            
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
            evec(1:breakx_idx) = line_left.yintr + (line_left.slope .* data(1:breakx_idx,1));
            evec(breakx_idx+1:ptct) = line_right.yintr + (line_right.slope .* data(breakx_idx+1:ptct,1));
            nvec = data(:,2) - evec;
            
            nvec = nvec.*nvec;
            dvec = dvec.*dvec;
            
            nsum = sum(nvec,1,'omitnan');
            dsum = sum(dvec,1,'omitnan');
            
            rsq = 1.0 - (nsum./dsum);
        end
        
    end
end