function dist2bound = getDist2Bound(rc_index, mask, CDboundary)

% rc_index = rc_index*resize_factor;
% mask = logical(zeros(size(CDboundary)));
% for frame_iter = 1:numberOfPages
%     temp_img = imresize(img_processed(:, :, frame_iter), resize_factor);  
%     mask(:, :, frame_iter) = temp_img>cutoff_new(frame_iter);
% end

    h = size(mask, 1); w = size(mask, 2); numberOfPages = size(mask, 3);
    dist2bound = zeros(1, numberOfPages);

    for frame_iter = 1:numberOfPages
        tempBoundary = CDboundary(:, :, frame_iter);
        tempLocs = find(tempBoundary > 0);

        if isempty(tempLocs)
            dist2bound(frame_iter) = nan;
        else
            [row,col] = ind2sub([h, w],tempLocs);
            [~, dist] = dsearchn([row,col], rc_index(frame_iter, :));
            if mask(round(rc_index(frame_iter, 1)), round(rc_index(frame_iter, 2)), frame_iter)
                dist2bound(frame_iter) = dist;
            else
                dist2bound(frame_iter) = -1*dist;
            end
        end
    end

end