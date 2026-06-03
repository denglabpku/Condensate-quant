function [labels_filtered, stats] = filter_condensates_Area_PC(labels, img, min_pixel_num, min_partition_coefficient, resize_factor)

% ===== 初始化 =====
labels_filtered = labels;
label_ids = unique(labels);
label_ids(label_ids == 0) = [];

stats = struct('label', {}, 'area', {}, 'PC_max', {}, 'PC_mean', {});

% ===== 遍历每个 condensate =====
for i = 1:length(label_ids)

    id = label_ids(i);
    mask = (labels == id);

    % ---- 面积 ----
    area_pixel = sum(mask(:));

    % ---- 构建 inside / outside (ring) ----
    mask_erode  = mask; %imerode(mask, strel('disk', resize_factor));
    mask_dilate1 = imdilate(mask, strel('disk', 2*resize_factor));
    mask_dilate2 = imdilate(mask, strel('disk', 3*resize_factor));

    inside  = mask_erode;
    outside = mask_dilate2 & ~mask_dilate1 & (labels == 0);

    % 避免空集合
    if sum(inside(:)) < 1 || sum(outside(:)) < 5
        labels_filtered(mask) = 0;
        continue;
    end

    inside_val  = img(inside);
    outside_val = img(outside);

    % ---- Partition coefficient（robust）----
    PC_max = max(inside_val) / mean(outside_val);
    PC_mean = mean(inside_val) / mean(outside_val);

    % ---- 保存统计 ----
    stats(end+1).label = id;
    stats(end).area = area_pixel;
    stats(end).PC_max = PC_max;
    stats(end).PC_mean = PC_mean;

    % ---- 筛选 ----
    if area_pixel < min_pixel_num || PC_max < min_partition_coefficient
        labels_filtered(mask) = 0;
    end

end

end