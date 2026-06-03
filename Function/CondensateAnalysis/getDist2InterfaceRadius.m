function [dist_p2boundary, dist_p2centroid, dist_boundary2centroid] = getDist2InterfaceRadius(rc_index, spots_resize, CDcenter, labels)

    % filter spots_resize in ROI
    h = size(labels, 1); w = size(labels, 2); numberOfPages = size(labels, 3);
    for frame_iter = 1:numberOfPages
        temp = spots_resize{frame_iter, 1};
        filter_idx = temp(:, 1)>0.5 & temp(:, 1)<h-0.5 & temp(:, 2)>0.5 & temp(:, 2)<w-0.5;
        spots_resize{frame_iter, 1} = temp(filter_idx, :);
    end

    for frame_iter = 1:numberOfPages
        temp_labels = labels(:, :, frame_iter);
        bw = temp_labels > 0;              % 所有非零区域
        cc = bwconncomp(bw, 4);       % 4邻域连通域
        
        newLabel = zeros(size(temp_labels), 'like', temp_labels);
        
        for i = 1:cc.NumObjects
            newLabel(cc.PixelIdxList{i}) = i;
        end
        labels(:, :, frame_iter) = newLabel;
    end

    for frame_iter = 1:numberOfPages
        
        tempCDcenter = CDcenter(:, :, frame_iter);
        temp_labels = labels(:, :, frame_iter);
        tempLocs = find(temp_labels > 0);

        [row,col] = ind2sub([h, w],tempLocs);
        [k, dist] = dsearchn([row,col], rc_index(frame_iter, :));
        bw = temp_labels==temp_labels(row(k), col(k));

        idx = bw(sub2ind([h, w], round(spots_resize{frame_iter,1}(:, 1)), round(spots_resize{frame_iter,1}(:, 2))));
        spots_select = spots_resize{frame_iter,1}(idx, :);

        tempCDcenter = tempCDcenter&bw;
        tempLocs = find(tempCDcenter > 0);
        if ~isempty(tempLocs)
            [r, c] = ind2sub([h, w],tempLocs);
            % center = [r, c];
            center = spots_select; %[r, c];
        else
            bwprops = regionprops(bw, 'Centroid');
            center = round(bwprops.Centroid);
        end
    
        B = bwboundaries(bw);        % 边界
        boundary = B{1};               % [row col]
        p = rc_index(frame_iter, :);

        % ---------- p -> nearest boundary ----------
        D_pb = vecnorm(boundary - p, 2, 2);
        [min_pb, ib] = min(D_pb);
        nearest_boundary = boundary(ib,:);

        % ---------- p -> nearest centroid ----------
        D_pc = vecnorm(center - p, 2, 2);
        [min_pc, ic] = min(D_pc);
        nearest_centroid = center(ic,:);
    
        dist_p2centroid(frame_iter) = min_pc;
    
        % ---------- nearest boundary <-> nearest centroid ----------
        dist_boundary2centroid(frame_iter) = norm(nearest_boundary - nearest_centroid);

        mask = temp_labels>0;
        if mask(ceil(rc_index(frame_iter, 1)), ceil(rc_index(frame_iter, 2)))
            dist_p2boundary(frame_iter) = min_pb;
        else
            dist_p2boundary(frame_iter) = -1*min_pb;
        end

    end

end
