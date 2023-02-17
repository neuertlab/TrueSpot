%
%%

function spotanno = SpotAnnoSnap_DeepBlink(spotanno)
    maxdist_3 = 5;
    maxdist_z = 3;
    minthidx_add = 6;
    
    T = size(spotanno.positives,1);
    
    for t = 1:T
        spotanno = spotanno.clearAtThreshold(t);
    end
    
    add_alloc = size(spotanno.positives{1,1}, 1);
    refspots = size(spotanno.ref_coords,1);
    temp_ref_tbl = zeros(refspots + add_alloc, 4);
    temp_ref_tbl(1:refspots,1:3) = spotanno.ref_coords(:,:);
    temp_tbl_sz = refspots;
    
    for t = T:-1:1
        if (t < minthidx_add) & (nnz(temp_ref_tbl(:,4)) == temp_tbl_sz)
            %All snapped. No more to add.
            break;
        end
        
        thpos = spotanno.positives{t,1};
        spotcount = size(thpos,1);
        pos_tbl = NaN(spotcount,6);
        pos_tbl(:,1:3) = thpos(:,1:3);
        
        for r = 1:refspots
            if (t < minthidx_add) & (temp_ref_tbl(r,4))
                continue;
            end
            
            x_dist = pos_tbl(:,1) - temp_ref_tbl(r,1);
            y_dist = pos_tbl(:,2) - temp_ref_tbl(r,2);
            z_dist = pos_tbl(:,3) - temp_ref_tbl(r,3);
            
            pos_tbl(:,4) = abs(z_dist);
            pos_tbl(:,5) = sqrt(x_dist.^2 + y_dist.^2);
            pos_tbl(:,6) = sqrt(x_dist.^2 + y_dist.^2 + z_dist.^2);
            
            match_bool = (pos_tbl(:,6) <= maxdist_3);
            match_bool = and(match_bool, (pos_tbl(:,4) <= maxdist_z));
            match_bool = and(match_bool, ~ismember(pos_tbl(:,1:3), temp_ref_tbl(:,1:3),'rows'));
            
            if nnz(match_bool) > 0
                [match_rows, ~] = find(match_bool);
                rmatches = pos_tbl(match_rows,:);
                
                if ~temp_ref_tbl(r,4)
                    %Pick closest spot to snap to.
                    [~,I] = min(rmatches(:,6));
                    old_pt = temp_ref_tbl(r,1:3);
                    temp_ref_tbl(r,1:3) = rmatches(I,1:3);
                    fprintf('DEBUG -- SpotAnnoSnap_DeepBlink -- (%d,%d,%d) snapped to (%d,%d,%d)\n', ...
                        old_pt(1,1), old_pt(1,2), old_pt(1,3),...
                        temp_ref_tbl(r,1), temp_ref_tbl(r,2), temp_ref_tbl(r,3));
                    temp_ref_tbl(r,4) = 1;
                end
                
                if t >= minthidx_add
                    %Add sufficiently close z spots to ref table, if not
                    %already added.
                    rz = temp_ref_tbl(r,3);
                    not_added = ~ismember(rmatches(:,1:3), temp_ref_tbl(:,1:3),'rows');
                    %not_added = not_added & (rmatches(:,3) ~= rz);
                    add_count = nnz(not_added);
                    if add_count > 0
                        for zz = (rz-maxdist_z):1:(rz+maxdist_z)
                            if zz == rz; continue; end
                            zmatch = find(and(not_added, rmatches(:,3) == zz));
                            if ~isempty(zmatch)
                                zset = rmatches(zmatch,:);
                                [~,min_i] = min(zset(:,6));
                                fprintf('DEBUG -- SpotAnnoSnap_DeepBlink -- adding (%d,%d,%d)\n', ...
                                    zset(min_i,1), zset(min_i,2), zset(min_i,3));
                                temp_tbl_sz = temp_tbl_sz+1;
                                temp_ref_tbl(temp_tbl_sz,1:3) = zset(min_i,1:3);
                                temp_ref_tbl(temp_tbl_sz,4) = 1;
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