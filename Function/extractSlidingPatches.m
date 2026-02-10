function [patches, coordinates] = extractSlidingPatches(img, patch_size, gap)
% img: 输入图像，大小为 [T, H, W]
% patch_size: 三元素数组 [pt, ph, pw]
% gap: 三元素数组 [gt, gh, gw]

[T, H, W] = size(img);
pt = patch_size(1); ph = patch_size(2); pw = patch_size(3);
gt = gap(1);        gh = gap(2);        gw = gap(3);

% 计算 patch 数量
num_t = ceil((T - pt + gt) / gt);
num_h = ceil((H - ph + gh) / gh);
num_w = ceil((W - pw + gw) / gw);

patches = {};  % 存放切块数据
coordinates = {};  % 存放坐标信息

count = 1;
for t = 0:(num_t - 1)
    for h = 0:(num_h - 1)
        for w = 0:(num_w - 1)
            % 时间方向
            if t < num_t - 1
                ts = t * gt + 1;
                te = ts + pt - 1;
            else
                te = T;
                ts = T - pt + 1;
            end

            % 高度方向
            if h < num_h - 1
                hs = h * gh + 1;
                he = hs + ph - 1;
            else
                he = H;
                hs = H - ph + 1;
            end

            % 宽度方向
            if w < num_w - 1
                ws = w * gw + 1;
                we = ws + pw - 1;
            else
                we = W;
                ws = W - pw + 1;
            end

            % 提取 patch
            patch = img(ts:te, hs:he, ws:we);
            patches{count} = patch;

            % 保存坐标
            coordinates{count} = struct( ...
                't_start', ts, 't_end', te, ...
                'h_start', hs, 'h_end', he, ...
                'w_start', ws, 'w_end', we ...
            );
            count = count + 1;
        end
    end
end
end
