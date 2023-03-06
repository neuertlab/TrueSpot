%
%%

function spotanno = SpotAnnoSnap_DeepBlink(spotanno)
    maxdist_3 = 4;
    maxdist_z = 2;
    minth_add = 0.95;
    minth_check = 0.85;
    
    T = size(spotanno.positives,1);
    thresh_x = transpose(spotanno.threshold_table(:,1));
    minthidx_add = RNAUtils.findThresholdIndex(minth_add, thresh_x);
    minthidx_check = RNAUtils.findThresholdIndex(minth_check, thresh_x);
    
    for t = 1:T
        spotanno = spotanno.clearAtThreshold(t);
    end
    
    add_alloc = size(spotanno.positives{1,1}, 1);
    refspots = size(spotanno.ref_coords,1);
    temp_ref_tbl = zeros(refspots + add_alloc, 4);
    temp_ref_tbl(1:refspots,1:3) = spotanno.ref_coords(:,:);
    temp_tbl_sz = refspots;
    
    z_added_tbl = false(refspots, (maxdist_z * 2) + 1); %Or just use bookkeeping :P
    
    for t = T:-1:minthidx_check
        if (t < minthidx_add) & (nnz(temp_ref_tbl(:,4)) == temp_tbl_sz)
            %All snapped. No more to add.
            break;
        end
        fprintf('DEBUG -- SpotAnnoSnap_DeepBlink -- Looking at spots >= prob %.2f\n', thresh_x(t));
        
        thpos = spotanno.positives{t,1};
        spotcount = size(thpos,1);
        pos_tbl = NaN(spotcount,6);
        pos_tbl(:,1:3) = thpos(:,1:3);
        
        for r = 1:refspots
            if (t < minthidx_add) & (temp_ref_tbl(r,4))
                continue;
            end
           % fprintf('DEBUG -- SpotAnnoSnap_DeepBlink -- Analyzing reference spot [%d of %d] @(%d,%d,%d)\n', ...
           %     r, refspots, temp_ref_tbl(r,1), temp_ref_tbl(r,2), temp_ref_tbl(r,3));
            
            x_dist = pos_tbl(:,1) - temp_ref_tbl(r,1);
            y_dist = pos_tbl(:,2) - temp_ref_tbl(r,2);
            z_dist = pos_tbl(:,3) - temp_ref_tbl(r,3);
            
            pos_tbl(:,4) = abs(z_dist);
            pos_tbl(:,5) = sqrt(x_dist.^2 + y_dist.^2);
            pos_tbl(:,6) = sqrt(x_dist.^2 + y_dist.^2 + z_dist.^2);
            
            match_bool = (pos_tbl(:,6) <= maxdist_3);
            match_bool = and(match_bool, (pos_tbl(:,4) <= maxdist_z));
            
            if nnz(match_bool) > 0
                %Only run ismember if there are still matches remaining! Maybe
                    %even run on subset of matches that have already been found?
                [match_rows, ~] = find(match_bool);
                rmatches = pos_tbl(match_rows,:);
                match_bool = ~ismember(rmatches(:,1:3), temp_ref_tbl(:,1:3),'rows');
                
                if nnz(match_bool) <= 0
                    %fprintf('DEBUG -- SpotAnnoSnap_DeepBlink -- No remaining matches found!\n');
                    continue;
                end
                
                [match_rows, ~] = find(match_bool);
                rmatches = rmatches(match_rows,:);
                %fprintf('DEBUG -- SpotAnnoSnap_DeepBlink -- Matches found: %d\n', nnz(match_bool));
                
                if ~temp_ref_tbl(r,4)
                    %Pick closest spot to snap to.
                    [~,I] = min(rmatches(:,6));
                    old_pt = temp_ref_tbl(r,1:3);
                    temp_ref_tbl(r,1:3) = rmatches(I,1:3);
                    fprintf('DEBUG -- SpotAnnoSnap_DeepBlink -- (%d,%d,%d) snapped to (%d,%d,%d)\n', ...
                        old_pt(1,1), old_pt(1,2), old_pt(1,3),...
                        temp_ref_tbl(r,1), temp_ref_tbl(r,2), temp_ref_tbl(r,3));
                    temp_ref_tbl(r,4) = 1;
                    
                    %Remove snap target from rmatches
                    rmatches(I,:) = [];
                end
                
                if t >= minthidx_add & ~isempty(rmatches)
                    %Add sufficiently close z spots to ref table, if not
                    %already added.
                    rz = temp_ref_tbl(r,3);
                    ctr = 0;
                    for zz = (rz-maxdist_z):1:(rz+maxdist_z)
                        ctr = ctr + 1;
                        if zz == rz; continue; end
                        
                        %See if spot has already been added for this z...
                        if ~z_added_tbl(r,ctr)
                            zmatch = find(rmatches(:,3) == zz);
                            if ~isempty(zmatch)
                                zset = rmatches(zmatch,:);
                                [~,min_i] = min(zset(:,6));
                                fprintf('DEBUG -- SpotAnnoSnap_DeepBlink -- adding (%d,%d,%d)\n', ...
                                    zset(min_i,1), zset(min_i,2), zset(min_i,3));
                                temp_tbl_sz = temp_tbl_sz+1;
                                temp_ref_tbl(temp_tbl_sz,1:3) = zset(min_i,1:3);
                                temp_ref_tbl(temp_tbl_sz,4) = 1;
                                z_added_tbl(r,ctr) = true;
                            end
                        end
                    end
                end
            end
        end
        
    end

    %Save new ref table
    spotanno.ref_coords = temp_ref_tbl(1:temp_tbl_sz,1:3);
    
end