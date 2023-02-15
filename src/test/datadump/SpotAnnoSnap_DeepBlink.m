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
        
        spotcount = size(spotanno.positives{t,1},1);
        pos_tbl = NaN(spotcount,6);
        pos_tbl(:,1:3) = spotanno.positives{t,1};
        
        for r = 1:temp_tbl_size
            if (t < minthidx_add) & (temp_ref_tbl(r,4))
                continue;
            end
            
            x_dist = pos_tbl(:,1) - temp_ref_tbl(r,1);
            y_dist = pos_tbl(:,2) - temp_ref_tbl(r,2);
            z_dist = pos_tbl(:,3) - temp_ref_tbl(r,3);
            
            pos_tbl(:,4) = z_dist;
            pos_tbl(:,5) = sqrt(x_dist.^2 + y_dist.^2);
            pos_tbl(:,6) = sqrt(x_dist.^2 + y_dist.^2 + z_dist.^2);
            
            match_bool = (pos_tbl(:,6) <= maxdist_3);
            match_bool = and(match_bool, (pos_tbl(:,4) <= maxdist_z));
            
            if nnz(match_bool) > 0
                [match_rows, ~] = find(match_bool);
                rmatches = pos_tbl(match_rows,:);
                
                if ~temp_ref_tbl(r,4)
                    %Pick closest spot to snap to.
                    [~,I] = min(rmatches(:,6));
                    old_pt = temp_ref_tbl(r,1:3);
                    temp_ref_tbl(r,1:3) = rmatches(I,1:3);
                    fprintf('DEBUG -- (%d,%d,%d) snapped to (%d,%d,%d)\n', ...
                        old_pt(1,1), old_pt(1,2), old_pt(1,3),...
                        temp_ref_tbl(r,1), temp_ref_tbl(r,2), temp_ref_tbl(r,3));
                    temp_ref_tbl(r,4) = 1;
                end
                
                if t >= minthidx_add
                    %Add sufficiently close z spots to ref table, if not
                    %already added.
                    not_added = ~ismember(rmatches(1:3), temp_ref_tbl(1:3));
                    add_count = nnz(not_added);
                    if add_count > 0
                        mcount = size(rmatches,1);
                        for j = 1:mcount
                            if ~not_added(j,1); continue; end
                            fprintf('DEBUG -- adding (%d,%d,%d)\n', ...
                                rmatches(j,1), rmatches(j,2), rmatches(j,3));
                            temp_tbl_sz = temp_tbl_sz+1;
                            temp_ref_tbl(temp_tbl_sz,1:3) = rmatches(j,1:3);
                            temp_ref_tbl(temp_tbl_sz,4) = 1;
                        end
                    end
                end
            end
        end
        
    end

end