function [spots, quality, I_log] = log_detector_3d_Sage(I, blob_diameter_xy, blob_diameter_z, threshold_value, mask)
    % I: 3D image stack (Y x X x Z)
    % sigmaX, sigmaY, sigmaZ: expected blob scale in each dimension (pixels)
    % Output: threshold_value and quality is just the I_log local maximal
    % value

    % Zuhui Wang based on Sage, D., F.R. Neumann, F. Hediger, S.M. Gasser, and M. Unser. 2005. “Automatic Tracking of Individual Fluorescence Particles: Application to the Study of Chromosome Dynamics.” IEEE Transactions on Image Processing 14 (9): 1372–83. https://doi.org/10.1109/TIP.2005.852787.

    I = double(I);
    [H, W, D] = size(I);

    if nargin<5
        mask = true(H, W, D);
    end
    if length(size(mask))==2
        mask = repmat(mask, [1, 1, D]);
    end

    sigmaX = blob_diameter_xy / (2*sqrt(2));
    sigmaY = sigmaX;
    sigmaZ = blob_diameter_z / (2*sqrt(2));

    % --- Build 1D Gaussian and 2nd derivative kernels for each axis ---
    % X axis
    szX = ceil(3*sigmaX);
    x = -szX:szX;
    gX = exp(-x.^2/(2*sigmaX^2)); gX = gX/sum(gX);
    g2X = (x.^2/sigmaX^4 - 1/sigmaX^2).*exp(-x.^2/(2*sigmaX^2));
    % g2X = g2X - mean(g2X);

    % Y axis
    szY = ceil(3*sigmaY);
    y = -szY:szY;
    gY = exp(-y.^2/(2*sigmaY^2)); gY = gY/sum(gY);
    g2Y = (y.^2/sigmaY^4 - 1/sigmaY^2).*exp(-y.^2/(2*sigmaY^2));
    % g2Y = g2Y - mean(g2Y);

    % Z axis
    szZ = ceil(3*sigmaZ);
    z = -szZ:szZ;
    gZ = exp(-z.^2/(2*sigmaZ^2)); gZ = gZ/sum(gZ);
    g2Z = (z.^2/sigmaZ^4 - 1/sigmaZ^2).*exp(-z.^2/(2*sigmaZ^2));
    % g2Z = g2Z - mean(g2Z);

    % --- Apply separable LoG ---
    % Term 1: second derivative in x, Gaussian in y,z
    I1 = imfilter(I, reshape(g2X,[numel(x),1,1]), 'replicate');
    I1 = imfilter(I1, reshape(gY,[1,numel(y),1]), 'replicate');
    I1 = imfilter(I1, reshape(gZ,[1,1,numel(z)]), 'replicate');

    % Term 2: second derivative in y, Gaussian in x,z
    I2 = imfilter(I, reshape(gX,[numel(x),1,1]), 'replicate');
    I2 = imfilter(I2, reshape(g2Y,[1,numel(y),1]), 'replicate');
    I2 = imfilter(I2, reshape(gZ,[1,1,numel(z)]), 'replicate');

    % Term 3: second derivative in z, Gaussian in x,y
    I3 = imfilter(I, reshape(gX,[numel(x),1,1]), 'replicate');
    I3 = imfilter(I3, reshape(gY,[1,numel(y),1]), 'replicate');
    I3 = imfilter(I3, reshape(g2Z,[1,1,numel(z)]), 'replicate');

    % Combine
    I_log = I1 + I2 + I3;

    % Invert for bright blobs
    I_log = -I_log; 

    % --- Detection pipeline (same as before) ---
    bw = imregionalmax(I_log,26);

    % Threshold based on LoG response (optional)
    % threshold = threshold_value * max(I_log(:));
    bw = bw & (I_log > threshold_value) & mask;
    
    % Non-maximum suppression using morphological operation
    minDist = mean([sigmaX,sigmaY,sigmaZ]);
    se = strel('sphere', round(minDist));
    I_supp = I_log .* bw;
    I_supp = I_supp .* (I_log == imdilate(I_log, se));

    [y,x,z] = ind2sub(size(I_supp), find(I_supp));
    pts = [x,y,z];    

    % Subpixel refinement
    offsets = zeros(size(pts));
    for i = 1:size(pts,1)
        x0 = pts(i,1); y0 = pts(i,2); z0 = pts(i,3);
        if x0>1 && y0>1 && z0>1 && x0<W && y0<H && z0<D
            patch = I_log(y0-1:y0+1, x0-1:x0+1, z0-1:z0+1);
            [dx, dy, dz] = subpixel_trackmate_3d(patch);
            if all(abs([dx,dy,dz])<1)
                offsets(i,:) = [dx, dy, dz];
            end
        end
    end

    spots = pts + offsets;

    % Quality = LoG response
    quality = zeros(size(spots,1),1);
    for i = 1:size(spots,1)
        x = round(spots(i,1));
        y = round(spots(i,2));
        z = round(spots(i,3));
        quality(i) = I_log(y,x,z);
    end
end

function [dx, dy, dz] = subpixel_quadratic_fit_3d(patch)
    % Fit quadratic surface to 3x3x3 patch
    [X,Y,Z] = ndgrid(-1:1,-1:1,-1:1);
    A = [X(:).^2, Y(:).^2, Z(:).^2, ...
         X(:).*Y(:), X(:).*Z(:), Y(:).*Z(:), ...
         X(:), Y(:), Z(:), ones(numel(X),1)];
    coeff = A\patch(:);

    % Coefficients: aX^2 + bY^2 + cZ^2 + dXY + eXZ + fYZ + gX + hY + iZ + j
    H = [2*coeff(1), coeff(4), coeff(5); ...
         coeff(4), 2*coeff(2), coeff(6); ...
         coeff(5), coeff(6), 2*coeff(3)];
    g = [coeff(7); coeff(8); coeff(9)];
    offset = -H\g;
    dx = offset(1); dy = offset(2); dz = offset(3);
end

function [dx,dy,dz] = subpixel_trackmate_3d(patch)

    % patch is 3x3x3, center = patch(2,2,2)

    % --- X direction ---
    a = patch(2,1,2);
    b = patch(2,2,2);
    c = patch(2,3,2);
    dx = (a - c) / (2*(a - 2*b + c));

    % --- Y direction ---
    a = patch(1,2,2);
    b = patch(2,2,2);
    c = patch(3,2,2);
    dy = (a - c) / (2*(a - 2*b + c));

    % --- Z direction ---
    a = patch(2,2,1);
    b = patch(2,2,2);
    c = patch(2,2,3);
    dz = (a - c) / (2*(a - 2*b + c));
end