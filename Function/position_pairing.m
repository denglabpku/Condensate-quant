function [paired_pos1, paired_pos2] = position_pairing(pos1,pos2,dist_thresh,use3D)

    if use3D
        [k1, distToPos1] = dsearchn(pos2, pos1);
        [k2, ~] = dsearchn(pos1, pos2);
    else
        [k1, distToPos1] = dsearchn(pos2(:, 1:2), pos1(:, 1:2));
        [k2, ~] = dsearchn(pos1(:, 1:2), pos2(:, 1:2));
    end

    pos1_pair = [(1:length(k1))', k1];
    pos2_pair = [(1:length(k2))', k2];

    pos1_pair(:, 3) = pos2_pair(pos1_pair(:, 2), 2);
    pos1_pair(:, 4) = distToPos1;

    idx_filter = (pos1_pair(:, 1) == pos1_pair(:, 3)) & (pos1_pair(:, 4)<=dist_thresh);

    pos1_pair = pos1_pair(idx_filter, :);

    paired_pos1 = pos1(pos1_pair(:, 1), :);
    paired_pos2 = pos2(pos1_pair(:, 2), :);

    

end