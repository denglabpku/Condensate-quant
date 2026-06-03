function [spots, quality] = log_detector_fft(I, blob_diameter, threshold_value, mask)

    if nargin<3
        threshold_value = 0;
    end

    % Convert to double
    I = double(I);
    [H, W] = size(I);

    if nargin<4
        mask = true(H, W);
    end

    % Estimate sigma (TrackMate convention)
    sigma = blob_diameter / (2*sqrt(2));

    % Fourier coordinates
    [kx, ky] = meshgrid( ...
        ifftshift((-floor(W/2)):(ceil(W/2)-1))/W, ...
        ifftshift((-floor(H/2)):(ceil(H/2)-1))/H );

    k2 = kx.^2 + ky.^2;

    % Build LoG filter in Fourier space
    LoG_kernel = - (2*pi)^2 * k2 .* exp(-2*(pi^2)*sigma^2*k2);

    % Apply filter via FFT
    F = fft2(I);
    I_log = real(ifft2(F .* LoG_kernel));

    % --- TrackMate uses inverted LoG for bright spots ---
    % (bright blobs -> negative response -> take negative to make peaks positive)
    I_log = -I_log;

    % Find local maxima (bright blobs)
    bw = imregionalmax(I_log);

    bw = bw & (I_log > threshold_value);

    % Threshold based on LoG response (optional)
    % threshold = 0.1 * max(I_log(:));
    % bw = bw & (I_log > threshold);

    % Extract coordinates
    [y, x] = find(bw);
    pts = [x, y];

    % Non-maximum suppression (too close points)
    minDist = sigma;
    keep = true(size(pts,1),1);
    for i = 1:size(pts,1)
        if ~keep(i), continue; end
        d = sqrt(sum((pts - pts(i,:)).^2, 2));
        tooClose = d < minDist & (1:size(pts,1))' > i;
        keep(tooClose) = false;
    end
    pts = pts(keep,:);

    % Subpixel quadratic refinement
    offsets = zeros(size(pts));
    for i = 1:size(pts,1)
        x0 = pts(i,1); y0 = pts(i,2);
        if x0>1 && y0>1 && x0<W && y0<H
            patch = I_log(y0-1:y0+1, x0-1:x0+1);
            [dx, dy] = subpixel_quadratic_fit(patch);
            if abs(dx)<1 && abs(dy)<1
                offsets(i,:) = [dx, dy];
            end
        end
    end

    % Final coordinates
    spots = [pts(:,1) + offsets(:,1), pts(:,2) + offsets(:,2)];

    % Compute quality (LoG response intensity)
    % TrackMate's definition: quality = LoG response at spot position (positive)
    quality = zeros(size(spots,1),1);
    for i = 1:size(spots,1)
        x = round(spots(i,1));
        y = round(spots(i,2));
        quality(i) = I_log(y,x);
    end

    % apply mask
    inside = mask(sub2ind(size(mask), round(spots(:, 2)), round(spots(:, 1))));
    spots = spots(inside, :);
    quality = quality(inside, :);
end

function [dx, dy] = subpixel_quadratic_fit(patch)
    % 3x3 quadratic fitting
    [X, Y] = meshgrid(-1:1, -1:1);
    Z = patch(:);
    A = [X(:).^2, Y(:).^2, X(:).*Y(:), X(:), Y(:), ones(9,1)];
    coeff = A\Z;
    a=coeff(1); b=coeff(2); c=coeff(3); d=coeff(4); e=coeff(5);
    H = [2*a, c; c, 2*b];
    g = [d; e];
    offset = -H\g;
    dx = offset(1);
    dy = offset(2);
end