function full_img = reconstructFromPatches(patches, coordinates, patch_size, img_size, gap)
% RECONSTRUCTFROMPATCHESHALFOVERLAP 将多个 patch 拼接为整图，只保留每个 patch 的有效中心区域
%
% 输入:
%   patches      : cell 数组，patch{i} 是大小为 [pt, ph, pw] 的 patch
%   coordinates  : cell 数组，对应每个 patch 的原图位置
%   patch_size   : [pt, ph, pw] 每个 patch 的大小
%   img_size     : [T, H, W] 原图大小
%   gap          : [gt, gh, gw] 滑动步长（小于 patch_size 表示有 overlap）
%
% 输出:
%   full_img     : 重建后的整图（无缝、无平均，直接拼）

% 计算每个方向需要裁剪多少（左右各一半 overlap）
cut_t1 = floor((patch_size(1) - gap(1)) / 2);
cut_t2 = ceil((patch_size(1) - gap(1)) / 2);

cut_h1 = floor((patch_size(2) - gap(2)) / 2);
cut_h2 = ceil((patch_size(2) - gap(2)) / 2);

cut_w1 = floor((patch_size(3) - gap(3)) / 2);
cut_w2 = ceil((patch_size(3) - gap(3)) / 2);

full_img = zeros(img_size, 'like', patches{1});

for i = 1:length(patches)
    patch = patches{i};
    coord = coordinates{i};

    % 计算 patch 内保留区域（去掉边缘 overlap）
    if coord.t_start == 1
        t1 = 1;
    else
        t1 = cut_t1 + 1;
    end
    if coord.t_end == img_size(1)
        t2 = patch_size(1);
    else
        t2 = patch_size(1) - cut_t2;
    end

    if coord.h_start == 1
        h1 = 1;
    else
        h1 = cut_h1 + 1;
    end
    if coord.h_end == img_size(2)
        h2 = patch_size(2);
    else
        h2 = patch_size(2) - cut_h2;
    end

    if coord.w_start == 1
        w1 = 1;
    else
        w1 = cut_w1 + 1;
    end
    if coord.w_end == img_size(3)
        w2 = patch_size(3);
    else
        w2 = patch_size(3) - cut_w2;
    end

    % patch 中裁剪后保留区域
    patch_crop = patch(t1:t2, h1:h2, w1:w2);

    % 将裁剪后的 patch 放回整图
    ts = coord.t_start + (t1 - 1);
    te = coord.t_start + (t2 - 1);
    hs = coord.h_start + (h1 - 1);
    he = coord.h_start + (h2 - 1);
    ws = coord.w_start + (w1 - 1);
    we = coord.w_start + (w2 - 1);

    full_img(ts:te, hs:he, ws:we) = patch_crop;
end
end

