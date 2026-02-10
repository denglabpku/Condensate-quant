%% IMAGE PROCESSING PIPELINE FOR BIOMOLECULAR CONDENSATE DYNAMICS
%  Project: Multi-color Analysis of RNA, OCT4 & BRD4 Condensates
%  Contact: Bo Wang, Peking University
%  
%  Description: 
%  This script processes multi-channel 3D/time-lapse TIFF images to quantify
%  statistics between RNA location biomolecular condensate (CD) morphology. The pipeline includes:
%    1. Deep-learning denoising via Noise2Void.
%    2. Richardson-Lucy Deconvolution.
%    3. Nucleus segmentation and HMRF-based condensate identification.
%    4. Morphometric statistics (Distance).

clc;close all;clear;
%% Denoising and Deconvolution parameter
% Pre-trained N2V network for OCT4/BRD4 live-SR data
onnx_path = '..\onnx\N2V_2D_LiveSR_125_1_E10_xy3z0.onnx';
is_GPU_avaliable = true;

% Initialize Network and extract patching requirements
net = importNetworkFromONNX(onnx_path);
inputsize = net.Layers(1).InputSize;
patch_h = inputsize(2);
patch_w = inputsize(3);
patch_t = inputsize(1);
disp('Load deep-learning denoising model complete!');
disp(['Input image size: H=', num2str(patch_h), '; W=', num2str(patch_w), '; T=', num2str(patch_t), '.']);

% Load Point Spread Function (PSF) for multi-channel deconvolution
psf_SR_488 = TIFFreader('..\PSF\psf_SR_channel488_2D.tif', 'double');
psf_SR_560 = TIFFreader('..\PSF\psf_SR_channel561_2D.tif', 'double');
psf_SR_640 = TIFFreader('..\PSF\psf_SR_channel642_2D.tif', 'double');

%% data loading & pre-processing
filepath_list = {''};

for filepath_iter = 1:length(filepath_list)

filepath = filepath_list{filepath_iter};
output_path = filepath;

filename_list = dir([filepath, '*.tif']);

pixelSize = 95;      % Physical pixel size in nanometers (nm)
resize_factor = 10;  % Sub-pixel interpolation factor for morphological precision
roi_width = 31;

for file_iter = 1:length(filename_list)

    filename = filename_list(file_iter).name;
    mkdir([output_path, filename(1:(end-4))]);
    disp(['Processing ', filename, ' ...']);

    %%%%%%%%%%%%%%%%%%%%%%%%%% read img sequence %%%%%%%%%%%%%%%%%%%%%%%%%%

    % Use Bio-Formats to import OME-TIFF metadata
    r = bfGetReader([filepath, filename]);
    omeMeta = r.getMetadataStore();
    
    sizeX = r.getSizeX();
    sizeY = r.getSizeY();
    sizeZ = r.getSizeZ();
    sizeC = r.getSizeC();
    sizeT = r.getSizeT();

    img_stack = zeros(sizeY, sizeX, sizeZ, sizeC, sizeT, 'uint16');
    for t = 1:sizeT
        for z = 1:sizeZ
            for c = 1:sizeC
                index = r.getIndex(z-1, c-1, t-1) + 1;
                img_stack(:,:,z,c,t) = bfGetPlane(r, index);
            end
        end
    end

    r.close(); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% denoising %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Adding overlap to eliminate edge artifacts
    img_stack_denoised = img_stack;

    if sizeY<patch_h || sizeX<patch_w || sizeT<patch_t
        warning('ImageTooSmall:SizeCheck', ...
            'The input image is too small (%dx%dx%d). Minimum required is %dx%dx%d.', ...
            sizeY, sizeX, sizeT, patch_h, patch_w, patch_t);
    end
    overlap_factor=0.25;
    gap_h = round(patch_h*(1-overlap_factor));% Patch gap in height
    gap_w = round(patch_w*(1-overlap_factor));% Patch gap in width
    gap_t = round(patch_t*(1-overlap_factor));% Patch gap in slice

    for z = 1:sizeZ
        for c = 1:sizeC
            % Intensity normalization (mean subtraction)
            img = single(squeeze(img_stack(:, :, z, c, :))); img_mean = mean(img(:));
            img_permute = permute(img, [3, 1, 2])-img_mean;

            % Sliding window patch extraction
            [patches, coordinates] = extractSlidingPatches(img_permute, [patch_t, patch_h, patch_w], [gap_t, gap_h, gap_w]);
            
            denoised_patches = {};
            h = waitbar(0, ['Processing Channel ', num2str(c), ' with Deep Learning ...']);
            for patch_iter = 1:length(patches)
                waitbar(patch_iter/length(patches), h, ['Processing Channel ', num2str(c), ' with Deep Learning ...']);
                if is_GPU_avaliable
                    I_dlarray = gpuArray(single(reshape(patches{patch_iter}, [patch_t, patch_h, patch_w, 1])));
                else
                    I_dlarray = single(reshape(patches{patch_iter}, [patch_t, patch_h, patch_w, 1]));
                end
                % Neural network inference
                denoised_patches{patch_iter} = predict(net, I_dlarray);
            end
            close(h); 

            % Reassemble patches and restore original intensity scale
            img_denoised_permute = reconstructFromPatches(denoised_patches, coordinates, [patch_t, patch_h, patch_w], size(img_permute), [gap_t, gap_h, gap_w]);
            img_denoised = permute(img_denoised_permute, [2, 3, 1])+img_mean;
            img_stack_denoised(:, :, z, c, :) = gather(img_denoised);
    
        end
    end
    disp('Image denoising complete!');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% deconvolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Richardson-Lucy algorithm
    pad = 20; % padding to remove edge strip
    img_stack_deconv = double(img_stack_denoised);
    img_stack_reconv = double(img_stack_denoised);

    % 488-channel deconvolution
    img_denoised = double(img_stack_denoised(:, :, :, 1, :));
    [img_deconv, img_reconv] = deconv_fixed_iter(img_denoised, psf_SR_488, pad, 20);
    img_stack_deconv(:, :, :, 1, :) = img_deconv;
    img_stack_reconv(:, :, :, 1, :) = img_reconv;

    % 560-channel deconvolution
    img_denoised = double(img_stack_denoised(:, :, :, 2, :));
    [img_deconv, img_reconv] = deconv_fixed_iter(img_denoised, psf_SR_560, pad, 20);
    img_stack_deconv(:, :, :, 2, :) = img_deconv;
    img_stack_reconv(:, :, :, 2, :) = img_reconv;

    % 640-channel deconvolution
    img_denoised = double(img_stack_denoised(:, :, :, 3, :));
    [img_deconv, img_reconv] = deconv_fixed_iter(img_denoised, psf_SR_640, pad, 20);
    img_stack_deconv(:, :, :, 3, :) = img_deconv;
    img_stack_reconv(:, :, :, 3, :) = img_reconv;

    % export denoised and deconvolved image stacks
    bfsave(uint16(img_stack_reconv), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-denoised.ome.tif']);
    img_stack_deconv(:, :, :, 1, :) = double(img_stack_denoised(:, :, :, 1, :));

    channel_488 = squeeze(img_stack_deconv(:, :, :, 1, :));
    channel_560 = squeeze(img_stack_deconv(:, :, :, 2, :));
    channel_640 = squeeze(img_stack_deconv(:, :, :, 3, :));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% nucleus mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate global nucleus mask using Channel 560
    nucleus_mask = NucleusMask2(max(channel_560, [], 3), 15);
    
    % Detect 3D spots using Laplacian of Gaussian (LoG)
    [spots_488, quality, I_log] = log_detector_3d_Sage(channel_488, 3, 3, mean(channel_488(:))/5);
    id = spots_488(:, 1)>5 & spots_488(:, 1)<(sizeX-5);
    spots_488 = spots_488(id, :); quality = quality(id);
    [~, i] = max(quality); spots_488 = spots_488(i, :);

    z_layer = round(spots_488(3));

    img_series_max = img_stack_deconv(:, :, z_layer, :, :);
    z_stack = 1; c_channel = sizeC; 
    h = sizeY; w = sizeX;
    numberOfPages = sizeT;
    channel_labels = {'RNA', 'OCT4', 'BRD4'};

    spots_488_3D = spots_488;

    save([output_path, filename(1:(end-4)), '.mat'], "img_series_max", "nucleus_mask", "resize_factor", "channel_labels", "pixelSize");
    nucleus_mask = repmat(nucleus_mask, 1, 1, numberOfPages);
    
    %% dealing with foci channel
    foci_result = struct();

    img = img_series_max(:, :, 1, 1);
    [spots_488, quality] = log_detector_fft(img, 3, 0, nucleus_mask);
    id = abs(spots_488(:, 1)-spots_488_3D(1))<3 & abs(spots_488(:, 2)-spots_488_3D(2))<3;
    spots_488 = spots_488(id, :); quality = quality(id);
    [~, i] = max(quality); spots_488 = spots_488(i, :);

    spots_488_2D = spots_488;

    foci_result(1).name = channel_labels{1};
    foci_result(1).spots_2D = spots_488_2D;
    foci_result(1).spots_3D = spots_488_3D;
    foci_result(1).rc_index = [spots_488_2D(:, 2), spots_488_2D(:, 1)]-0.5;

    for c_iter = 1 % RNA channel
    
        disp(['Processing ', channel_labels{c_iter}, ' channel ...']);       
        img_series = squeeze(img_series_max(:, :, :, c_iter, :));
        rc_index = foci_result(c_iter).rc_index;      
        h = size(img_series, 1); w = size(img_series, 2);
        numberOfPages = size(img_series, 3);

        % calculate foci intensity relative to background
        [rna_bkg,base_bkg] = IntensityCalculation(img_series,nucleus_mask, rc_index, 4, 5);
        norm_base_bkg = movmean(base_bkg, 10)-500;
        intensity = zeros(size(norm_base_bkg));
        for frame_iter = 1:numberOfPages
            intensity(frame_iter) = rna_bkg(frame_iter, 1)/norm_base_bkg(frame_iter)*norm_base_bkg(1);
        end
    
        % define roi_window and roi_resize for roi selection
        if c_iter==1  %  ROI selection
            roi_window = logical(zeros(size(nucleus_mask)));
            for frame_iter = 1:numberOfPages
                roi_window(ceil(rc_index(frame_iter, 1)), ceil(rc_index(frame_iter, 2)), frame_iter) = 1;
            end
            roi_resize = imresize(imdilate(roi_window, true(roi_width, roi_width)), resize_factor, "nearest");
        end
    
        % export processed DNA and RNA channel images
        rc_index_rescale = rc_index*resize_factor;
        img_processed_roi = zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages);
        img_center = logical(zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages));
        for frame_iter = 1:numberOfPages
            disp(['Processing Frame ', num2str(frame_iter), ' ...']);
            temp_img = imresize(img_series(:, :, frame_iter), resize_factor, 'nearest');

            refined_bw = false(size(temp_img));
            refined_bw(ceil(rc_index_rescale(frame_iter, 1)), ceil(rc_index_rescale(frame_iter, 2))) = 1;
            temp_mask = imdilate(refined_bw, strel('disk', 3));
            temp_roi_resize = roi_resize(:, :, frame_iter);
        
            [row1, row2, col1, col2] = getROIboundary(temp_roi_resize, roi_width*resize_factor);
            img_processed_roi(row1:row2, col1:col2, frame_iter) = reshape(temp_img(temp_roi_resize), [row2-row1+1, col2-col1+1]);
            img_center(row1:row2, col1:col2, frame_iter) = reshape(temp_mask(temp_roi_resize), [row2-row1+1, col2-col1+1]);
        end
        img_processed_roi_export = imresize(img_processed_roi, 1/resize_factor, "nearest");
        TIFwriter(uint16(img_processed_roi_export), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '.tif']);
        TIFwriter(uint8(img_center), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-Center-roi.tif'], 'lzw');
    
        % important parameter: rc_index, bw, intensity
        foci_result(c_iter).rna_bkg = rna_bkg;
        foci_result(c_iter).base_bkg = base_bkg;
        foci_result(c_iter).intensity = intensity;
    end
    
    save([output_path, filename(1:(end-4)), '.mat'], "foci_result", "roi_window", '-append');

%% dealing with condensate channel
condensate_result = struct; 
seg_points = [5, 5];% Threshold for condensate class, above this will be considered as condensates

for c_iter = 2:3 % OCT4 and BRD4 channel

    img_series = squeeze(img_series_max(:, :, :, c_iter, :));

    %%%%%%%%%%%%%%%% get CD center, interface and boundary %%%%%%%%%%%%%%%%
    
    CDcenter = logical(zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages));
    CDinterface = CDcenter; CDboundary = CDcenter; CDmask = CDboundary;

    img_processed_roi = zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages);
    img_processed_bicubic_roi = zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages);
    labels = uint16(zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages));
    spots = cell(numberOfPages, 1);

    for frame_iter = 1:numberOfPages
        disp(['Processing Frame ', num2str(frame_iter), ' ...']);
        temp_img = imresize(img_series(:, :, frame_iter), resize_factor);
        temp_roi_resize = roi_resize(:, :, frame_iter);
        [row1, row2, col1, col2] = getROIboundary(temp_roi_resize, roi_width*resize_factor);
        img_processed_bicubic_roi(row1:row2, col1:col2, frame_iter) = reshape(temp_img(temp_roi_resize), [row2-row1+1, col2-col1+1]);

        % Hidden Markov Random Field (HMRF) to include spatial dependency
        temp_nucleus_mask = imresize(nucleus_mask(:, :, frame_iter), resize_factor, "nearest");
        temp_nucleus_mask = reshape(temp_nucleus_mask(temp_roi_resize), [row2-row1+1, col2-col1+1]);
        nclust = 7;% Number of intensity clusters;
        seg_point = seg_points(c_iter-1);% Threshold for condensate class, above this will be considered as condensates
        temp_img_roi = img_processed_bicubic_roi(:, :, frame_iter);
        [HMRFseg, ~] = HMRFseg4img(temp_img_roi, temp_nucleus_mask, nclust, 0.1, 10^(-8));
        bw_HMRF = HMRFseg.img_class>=seg_point;
        local_thresh = min(temp_img_roi(bw_HMRF));
        bw_local = bw_HMRF;

        % CD center, interface, boundary and labels identification
        [temp_spots, quality] = log_detector_fft(img_series(:, :, frame_iter), 3, 0, nucleus_mask & img_series(:, :, frame_iter)>=local_thresh);
        spots{frame_iter, 1} = [temp_spots(:, 2), temp_spots(:, 1)]-0.5;
        temp_CDcenter = false(size(temp_img));
        spots_resize = ceil((temp_spots-0.5)*resize_factor);
        for i = 1:size(temp_spots, 1)
            temp_CDcenter(spots_resize(i, 2), spots_resize(i, 1)) = 1;
        end
        CDcenter(row1:row2, col1:col2, frame_iter) = reshape(temp_CDcenter(temp_roi_resize), [row2-row1+1, col2-col1+1]);    
        CDinterface(:, :, frame_iter) = getCDinterface(img_processed_bicubic_roi(:, :, frame_iter), bw_local, 0);
        [temp_CDboundary, temp_labels] = getCDboundary(img_processed_bicubic_roi(:, :, frame_iter), bw_local, CDcenter(:, :, frame_iter), CDinterface(:, :, frame_iter), 0);
        CDboundary(:, :, frame_iter) = temp_CDboundary;
        labels(:, :, frame_iter) = temp_labels;
        CDmask(:, :, frame_iter) = bw_local;
    
        temp_img = imresize(img_series(:, :, frame_iter), resize_factor, "nearest");
        img_processed_roi(row1:row2, col1:col2, frame_iter) = reshape(temp_img(temp_roi_resize), [row2-row1+1, col2-col1+1]);

        % export HMRF segmentation and visualization
        if frame_iter == 1
            fig1 = figure;
            fig1.Units = "inches";
            fig1.Position = [7.4,3.9,9.8,6.4];
            subplot(1, 2, 1);
            imagesc(img_processed_roi(:, :, frame_iter));
            colormap(gca, "gray")
            daspect([1, 1, 1]);
            axis off
            subplot(1, 2, 2);
            cmap = uint8([4, 0, 0; 56, 46, 142; 137, 48, 141; 215, 31, 40; 239, 127, 25; 244, 191, 27; 244, 237, 70; 255, 255, 255]);
            nucleus_partition = HMRFseg.img_class;
            imagesc(nucleus_partition);
            if min(min(nucleus_partition)) == 1
                colormap(gca, cmap(2:end, :));
            else
                colormap(gca, cmap);
            end
            axis off
            daspect([1, 1, 1]);
            print(fig1, [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-HMRFseg.png'], '-dpng');
            close;
        end
    end
    CDmask_export = uint8(CDmask)*50;
    CDmask_export(CDboundary)=255;
    
    % Export morphological results as TIFF
    img_processed_roi_export = imresize(img_processed_roi, 1/resize_factor, "nearest");
    TIFwriter(uint16(img_processed_roi_export), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-roi.tif']);
    TIFwriter(uint16(img_processed_bicubic_roi), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-roi-bicubic.tif']);
    TIFwriter(uint8(CDcenter), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-CDcenter.tif'], 'lzw');
    TIFwriter(uint8(CDinterface), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-CDinterface.tif'], 'lzw');
    TIFwriter(uint8(CDboundary), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-CDboundary.tif'], 'lzw');
    TIFwriter(uint8(CDmask_export), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-CDmask.tif']);

    % Store results in struct
    condensate_result(c_iter).name = channel_labels{c_iter};
    condensate_result(c_iter).spots_2d = spots;
    condensate_result(c_iter).img_processed_roi = img_processed_roi;
    condensate_result(c_iter).img_processed_bicubic_roi = img_processed_bicubic_roi;
    condensate_result(c_iter).CDcenter = CDcenter;
    condensate_result(c_iter).CDinterface = CDinterface;
    condensate_result(c_iter).CDboundary = CDboundary;
    condensate_result(c_iter).CDmask = CDmask;
    condensate_result(c_iter).labels = labels;

    %%% calculate the distance to centroid and boundary within roi mask %%%

    spots_resize = spots;
    dist2bound = zeros(1, numberOfPages);
    dist2center = zeros(1, numberOfPages);
    condensate_radius = zeros(1, numberOfPages);
    equiv_radius = zeros(1, numberOfPages);
    boundary2center = zeros(1, numberOfPages);
    dist2boundary = zeros(1, numberOfPages);
    dist2centroid = zeros(1, numberOfPages);
    boundary2centroid = zeros(1, numberOfPages);
    for f_iter = 1 % RNA channel
        rc_index = foci_result(f_iter).rc_index*resize_factor;
        for frame_iter = 1:numberOfPages
            [row, col] = find(roi_resize(:, :, frame_iter));
            min_row = min(row); min_col = min(col);
            rc_index(frame_iter, :) = rc_index(frame_iter, :)-[min_row-1, min_col-1]+0.5;
            spots_resize{frame_iter, 1} = spots{frame_iter, 1}*resize_factor-[min_row-1, min_col-1]+0.5;
        end
        % Metric calculations: distance from RNA to condensate boundary/center
        dist2bound(f_iter, :) = getDist2Bound(rc_index, CDmask, CDboundary)*pixelSize/resize_factor;

        [dist, radius1, radius2, radius3, spots_resize] = getDist2CenterRadius(rc_index, spots_resize, CDcenter, labels);
        dist2center(f_iter, :) = dist*pixelSize/resize_factor;
        condensate_radius(f_iter, :) = radius1*pixelSize/resize_factor;
        equiv_radius(f_iter, :) = radius2*pixelSize/resize_factor;
        boundary2center(f_iter, :) = radius3*pixelSize/resize_factor;
        [dist_2boundary, dist_2centroid, boundary_2centroid] = getDist2InterfaceRadius(rc_index, spots_resize, CDcenter, labels);
        dist2boundary(f_iter, :) = dist_2boundary*pixelSize/resize_factor;
        dist2centroid(f_iter, :) = dist_2centroid*pixelSize/resize_factor;
        boundary2centroid(f_iter, :) = boundary_2centroid*pixelSize/resize_factor;

        condensate_result(f_iter).spots_resize = rc_index;
    end

    condensate_result(c_iter).spots_resize = spots_resize;
    condensate_result(c_iter).dist2bound = dist2bound;
    condensate_result(c_iter).dist2center = dist2center;
    condensate_result(c_iter).condensate_radius = condensate_radius;
    condensate_result(c_iter).equiv_radius = equiv_radius;
    condensate_result(c_iter).boundary2center = boundary2center;
    condensate_result(c_iter).dist2boundary = dist2boundary;
    condensate_result(c_iter).dist2centroid = dist2centroid;
    condensate_result(c_iter).boundary2centroid = boundary2centroid;

end

% export condensate segmentation result and centroids summary plot
save([output_path, filename(1:(end-4)), '.mat'], "condensate_result", '-append');

% visualization
frame_iter = 1;
fig = figure('Color','w','Position',[200 200 1350 680]);
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
imagesc(condensate_result(2).img_processed_roi(:, :, frame_iter))
axis image off
colormap gray
hold on
visboundaries(condensate_result(2).CDmask(:, :, frame_iter), 'Color', 'g', 'LineWidth', 1);
plot(condensate_result(1).spots_resize(2), condensate_result(1).spots_resize(1), 'bo','MarkerSize',7,'LineWidth',1.2);
spots = condensate_result(2).spots_resize{frame_iter, 1};
plot(spots(:,2), spots(:,1), 'go','MarkerSize',7,'LineWidth',1.2);
title('OCT4','FontSize',11,'FontWeight','bold')
hold off

nexttile
imagesc(condensate_result(3).img_processed_roi(:, :, frame_iter))
axis image off
colormap gray
hold on
visboundaries(condensate_result(3).CDmask(:, :, frame_iter), 'Color', 'm', 'LineWidth', 1);
plot(condensate_result(1).spots_resize(2), condensate_result(1).spots_resize(1), 'bo', 'MarkerSize',7,'LineWidth',1.2);
spots = condensate_result(3).spots_resize{frame_iter, 1};
plot(spots(:,2), spots(:,1), 'mo', 'MarkerSize',7,'LineWidth',1.2);
title('BRD4','FontSize',11,'FontWeight','bold')
hold off

print(fig, [output_path, filesep, filename(1:(end-4)), '.png'], '-dpng');

%%
close all;
end

end

%%

%% calculate distance of DNA and RNA

filepath = '';
filename_list = dir([filepath, '*.mat']);
mkdir([filepath, 'result']);

for file_iter = 6%1:length(filename_list)
    filename = filename_list(file_iter).name;

    load([filepath, filename], "foci_result", "condensate_result", "channel_labels");

    pixelSize = 160;

% plot distance to condensates boundary

time_used = 294;

fig1 = figure;
fig1.Units = 'inches'; fig1.Position = [10.3,5.3,5.8,8.3];
% RNA intensity
subplot(2, 1, 1);
plot(movmean(foci_result(1).intensity, 5));
ylabel('Intensity (A.U.)');
title('RNA intensity (A.U.)');
xlabel('Frame(2s)');
% xlim([1, length(foci_result(2).intensity)]);
xlim([1, time_used]);
ylim([0, 5e4]);

subplot(2, 1, 2);
y = movmean(condensate_result(2).dist2bound(1, :), 5);
% y = condensate_result(3).dist2bound(f_iter, :);
x = 1:length(y);
plot(x, y, 'Color', 'green');
hold on;
a1=fill([x, fliplr(x)], [y, zeros(size(y))], 'green', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
y = movmean(condensate_result(3).dist2bound(1, :), 5);
% y = condensate_result(4).dist2bound(f_iter, :);
x = 1:length(y);
plot(x, y, 'Color', 'magenta');
a2=fill([x, fliplr(x)], [y, zeros(size(y))], 'magenta', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
line([1, length(y)], [0, 0], 'Color','red','LineStyle','--');
hold off
xlabel('Frame(2s)');
ylabel('Distance(nm)');
title([channel_labels{1}, ' distance to condensates boundary(nm)']);
legend([a1, a2], {'OCT4', 'BRD4'});
% xlim([1, length(y)]);
xlim([1, time_used]);
ylim([-200, 400]);


% print(fig1, [filepath, 'result', filesep, filename(1:(end-4)), '-smooth.png'], '-dpng');
% close;
end