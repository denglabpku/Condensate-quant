function [dist2center, radius1, radius2, radius3, spots_resize] = getDist2CenterRadius(rc_index, spots_resize, CDcenter, labels)

    % filter spots_resize in ROI
    h = size(labels, 1); w = size(labels, 2); numberOfPages = size(labels, 3);
    for frame_iter = 1:numberOfPages
        temp = spots_resize{frame_iter, 1};
        filter_idx = temp(:, 1)>0.5 & temp(:, 1)<h-0.5 & temp(:, 2)>0.5 & temp(:, 2)<w-0.5;
        spots_resize{frame_iter, 1} = temp(filter_idx, :);
    end

    for frame_iter = 1:numberOfPages
        
        tempCDcenter = CDcenter(:, :, frame_iter);
        temp_labels = labels(:, :, frame_iter);
        tempLocs = find(temp_labels > 0);

        [row,col] = ind2sub([h, w],tempLocs);
        [k, dist] = dsearchn([row,col], rc_index(frame_iter, :));
        bw = temp_labels==temp_labels(row(k), col(k));
        labelprop = regionprops(bw, 'EquivDiameter');

        idx = bw(sub2ind([h, w], round(spots_resize{frame_iter,1}(:, 1)), round(spots_resize{frame_iter,1}(:, 2))));
        spots_select = spots_resize{frame_iter,1}(idx, :);

        tempCDcenter = tempCDcenter&bw;
        tempLocs = find(tempCDcenter > 0);
        bwprops = regionprops(bw, 'Centroid');
        tempcenter = round(bwprops.Centroid);
        if ~isempty(tempLocs) && ~isempty(spots_select)
            % [r, c] = ind2sub([h, w],tempLocs);
            D = vecnorm(spots_select - tempcenter, 2, 2);
            [~, idc] = min(D);
            center = spots_select(idc,:);  
        else
            center = tempcenter;
        end
    
        B = bwboundaries(bw);        % 边界
        boundary = B{1};               % [row col]
        p = rc_index(frame_iter, :);
        
        v = p - center;
        v = v / norm(v);
        
        t = (boundary - center) * v';  % 投影到方向上
        proj = center + t.*v;
        
        dist_perp = vecnorm(boundary - proj,2,2);
        
        [~, id] = min(dist_perp + 1e6*(t<0)); % 只取前向
        intersect = boundary(id,:);

        % ===== 新增：rc_index 最近的 boundary 点 =====
        d_bp = vecnorm(boundary - p, 2, 2);      % boundary 到 rc_index 的距离
        [~, id2] = min(d_bp);
        nearest_boundary = boundary(id2, :);

        
        radius1(frame_iter) = norm(intersect - center);
        radius2(frame_iter) = labelprop.EquivDiameter/2;
        radius3(frame_iter) = norm(nearest_boundary - center);
        dist2center(frame_iter) = norm(p - center);

    end

end
