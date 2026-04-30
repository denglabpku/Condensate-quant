%% IMAGE PROCESSING PIPELINE FOR BIOMOLECULAR CONDENSATE DYNAMICS
%  Project: Multi-color Analysis of DNA, RNA, OCT4 & BRD4 Condensates
%  Contact: Bo Wang, Peking University
%  
%  Description: 
%  This script processes multi-channel 3D/time-lapse TIFF images to quantify
%  statistics between RNA location biomolecular condensate (CD) morphology. The pipeline includes:
%    1. Deep-learning denoising via Noise2Void.
%    2. Richardson-Lucy Deconvolution.
%    3. Nucleus segmentation and HMRF-based condensate identification.
%    4. Morphometric statistics (Distance).
%    5. Visualization (scatter plot, contour plot)
clc;close all;clear;rng(42)
%% Denoising and Deconvolution parameter
% Pre-trained N2V network for OCT4/BRD4 live-SR data
onnx_path = '..\onnx\N2V_2D_OCT4_BRD4_liveSR_125_1_E10_xy3z0.onnx';
is_GPU_avaliable = true;

% Initialize Network and extract patching requirements
inputsize = net.Layers(1).InputSize;
patch_h = inputsize(2);
patch_w = inputsize(3);
patch_t = inputsize(1);
disp('Load deep-learning denoising model complete!');
disp(['Input image size: H=', num2str(patch_h), '; W=', num2str(patch_w), '; T=', num2str(patch_t), '.']);

% Load Point Spread Function (PSF) for multi-channel deconvolution
psf_SR_405 = TIFFreader('..\PSF\psf_SR_channel405_2D.tif', 'double');
psf_SR_488 = TIFFreader('..\PSF\psf_SR_channel488_2D.tif', 'double');
psf_SR_560 = TIFFreader('..\PSF\psf_SR_channel561_2D.tif', 'double');
psf_SR_640 = TIFFreader('..\PSF\psf_SR_channel642_2D.tif', 'double');

%% data loading and pre-processing
filepath_list = {'.\Promoter_RNA_SCR_BRD4', ...
                 '.\Promoter_RNA_SCR_OCT4'};
condensate_name_list = {'BRD4', 'OCT4'};
channel_labels = {'Promoter', 'RNA', 'Enhancer', 'Condensate'};
for filepath_iter = 1:length(filepath_list)

filepath = filepath_list{filepath_iter};

condensate_name = condensate_name_list{filepath_iter};

for with_RNA = [1, 0]

if with_RNA
    output_path = [filepath, 'filter_result_with_RNA_averaged_boundary', filesep];
    mkdir(output_path)

    filename_list = readlines([filepath, filesep, 'Cell_with_RNA.txt']);
    filename_list = strrep(filename_list, "'", "");
else
    output_path = [filepath, 'filter_result_without_RNA_averaged_boundary', filesep];
    mkdir(output_path)

    filename_list = readlines([filepath, filesep, 'Cell_without_RNA.txt']);
    filename_list = strrep(filename_list, "'", "");
end

filename_list = filename_list(strlength(filename_list) > 0); % remove empty string

pixelSize = 95;      % Physical pixel size in nanometers (nm)
resize_factor = 10;  % Sub-pixel interpolation factor for morphological precision
roi_width = 31;
%%
for file_iter = 1:length(filename_list)

    filename = filename_list{file_iter};
    filename = [filename(1:end-3),'tif'];
    % filename = filename_list(file_iter).name;
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

    if sizeY<patch_h || sizeX<patch_w
        warning('ImageTooSmall:SizeCheck', ...
            'The input image is too small (%dx%d). Minimum required is %dx%d.', ...
            sizeY, sizeX, patch_h, patch_w);
    end
    overlap_factor=0.25;
    gap_h = round(patch_h*(1-overlap_factor));% Patch gap in height
    gap_w = round(patch_w*(1-overlap_factor));% Patch gap in width
    gap_t = round(patch_t*(1-overlap_factor));% Patch gap in slice

    for c = 1:sizeC
        % Intensity normalization (mean subtraction)
        img = single(img_stack(:, :, :, c)); img_mean = mean(img(:));
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
        img_stack_denoised(:, :, :, c) = gather(img_denoised);

    end
    disp('Image denoising complete!');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% deconvolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Richardson-Lucy algorithm
    % Fixed at 20 iterations based on previous test, should be adjusted according to your data typ
    pad = 20; % padding to remove edge strip
    img_stack_deconv = double(img_stack_denoised);
    img_stack_reconv = double(img_stack_denoised);

    layer_select = [];
    for layer_iter = 1:sizeZ
        img = img_stack_denoised(:, :, layer_iter, 1);
        if mean(img(:))>=500
            layer_select = [layer_select, layer_iter];
        end
    end

    % 405-channel deconvolution
    img_denoised = double(img_stack_denoised(:, :, layer_select, 1));
    [img_deconv, img_reconv] = deconv_fixed_iter(img_denoised, psf_SR_405, pad, 20);
    img_stack_deconv(:, :, layer_select, 1) = img_deconv;
    img_stack_reconv(:, :, layer_select, 1) = img_reconv;

    % 488-channel deconvolution
    img_denoised = double(img_stack_denoised(:, :, :, 2));
    [img_deconv, img_reconv] = deconv_fixed_iter(img_denoised, psf_SR_488, pad, 20);
    img_stack_deconv(:, :, :, 2) = img_deconv;
    img_stack_reconv(:, :, :, 2) = img_reconv;

    % 560-channel deconvolution
    img_denoised = double(img_stack_denoised(:, :, 1:(sizeZ-1), 3));
    [img_deconv, img_reconv] = deconv_fixed_iter(img_denoised, psf_SR_560, pad, 20);
    img_stack_deconv(:, :, 1:(sizeZ-1), 3) = img_deconv;
    img_stack_reconv(:, :, 1:(sizeZ-1), 3) = img_reconv;

    % 640-channel deconvolution
    img_denoised = double(img_stack_denoised(:, :, 1:(sizeZ-1), 4));
    [img_deconv, img_reconv] = deconv_fixed_iter(img_denoised, psf_SR_640, pad, 20);
    img_stack_deconv(:, :, 1:(sizeZ-1), 4) = img_deconv;
    img_stack_reconv(:, :, 1:(sizeZ-1), 4) = img_reconv;

    % export denoised and deconvolved image stacks
    bfsave(uint16(img_stack_deconv), [output_path, filename(1:(end-4)), '-denoised-deconv.ome.tif']);
    bfsave(uint16(img_stack_reconv), [output_path, filename(1:(end-4)), '-denoised.ome.tif']);
    img_stack_deconv(:, :, :, 1) = double(img_stack_denoised(:, :, :, 1));
    img_stack_deconv(:, :, :, 2) = double(img_stack_denoised(:, :, :, 2));
    img_stack_deconv(:, :, :, 3) = double(img_stack_denoised(:, :, :, 3));

    channel_405 = img_stack_deconv(:, :, :, 1);
    channel_488 = img_stack_deconv(:, :, :, 2);
    channel_560 = img_stack_deconv(:, :, 1:(sizeZ-1), 3);
    channel_640 = img_stack_deconv(:, :, 1:(sizeZ-1), 4);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% nucleus mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate global nucleus mask using Channel 560
    nucleus_mask = NucleusMask2(mean(img_stack_reconv(:, :, layer_select, 4), 3), 15);

    temp_channel_405 = channel_405(:, :, layer_select);
    % Detect 3D spots using Laplacian of Gaussian (LoG)
    [spots_405, quality, I_log] = log_detector_3d_Sage(temp_channel_405, 3, 3, mean(channel_405(:))/5);
    id = spots_405(:, 1)>5 & spots_405(:, 1)<(sizeX-5);
    spots_405 = spots_405(id, :); quality = quality(id);
    [~, i] = max(quality); spots_405 = spots_405(i, :);
    spots_405(3) = spots_405(3)+layer_select(1)-1;

    [spots_560, quality, I_log] = log_detector_3d_Sage(channel_560, 3, 3, 0);
    threshold = 2*median(quality);
    spots_560 = spots_560(quality>threshold, :);
    [k2, dist] = dsearchn(spots_560, spots_405);
    spots_560 = spots_560(k2, :);
    
    if with_RNA
        [spots_488, quality, I_log] = log_detector_3d_Sage(channel_488, 3, 3, 0);
        % threshold = trackmate_auto_threshold(quality, 10);
        threshold = 2*median(quality);
        spots_488 = spots_488(quality>threshold, :);
        if size(spots_488, 1)>=1
            [k1, dist] = dsearchn(spots_488(:, 1:2), spots_405(:, 1:2));
            % 
            % spots_488_temp = spots_488;
            % spots_488_temp(:,3) = spots_488(:,3)*3;
            % spots_405_temp = spots_405;
            % spots_405_temp(:,3) = spots_405(:,3)*3;
            % [k1, dist] = dsearchn(spots_488_temp(:, 1:3), spots_405_temp(:, 1:3));
            
            if dist<5 %pixel
                spots_488 = spots_488(k1, :);
            else
                spots_488 = spots_405; % no RNA
            end
        else
            spots_488 = spots_405; % no RNA
        end
    else
        spots_488 = spots_405; % no RNA
    end

    z_layer_405 = round(spots_405(3));
    z_layer_488 = round(spots_488(3));
    z_layer_560 = round(spots_560(3));

    z_layer_used = max(min([z_layer_405, z_layer_488, z_layer_560]), layer_select(1)):min(max([z_layer_405, z_layer_488, z_layer_560]), sizeZ-1);
    
    if isempty(z_layer_used)
        if z_layer_405 == sizeZ
            z_layer_used = sizeZ-1;
        elseif z_layer_405 < layer_select(1)
            z_layer_used = layer_select(1);
        end
    end
    
    % define 'final' img_series
    img_series_max = zeros(sizeY, sizeX, 1, sizeC);
    img_series_max(:, :, :, 1) = channel_405(:, :, z_layer_405);
    img_series_max(:, :, :, 2) = channel_488(:, :, z_layer_488);
    img_series_max(:, :, :, 3) = channel_560(:, :, z_layer_560);
    img_series_max(:, :, :, 4) = mean(channel_640(:, :, z_layer_used), 3);

    z_stack = 1; c_channel = sizeC; 
    h = sizeY; w = sizeX;
    numberOfPages = sizeT;

    spots_405_3D = spots_405;
    spots_488_3D = spots_488;
    spots_560_3D = spots_560;

    save([output_path, filename(1:(end-4)), '.mat'], "img_series_max", "nucleus_mask", "resize_factor", "channel_labels", "pixelSize");

    %%
    foci_result = struct();
    %% dealing with foci channel

    img = img_series_max(:, :, :, 1);
    [spots_405, quality] = log_detector_fft(img, 3, 0, nucleus_mask);
    id = abs(spots_405(:, 1)-spots_405_3D(1))<3 & abs(spots_405(:, 2)-spots_405_3D(2))<3;
    spots_405 = spots_405(id, :); quality = quality(id);
    [~, i] = max(quality); spots_405 = spots_405(i, :);

    img = img_series_max(:, :, :, 3);
    [spots_560, quality] = log_detector_fft(img, 3, 0, nucleus_mask);
    threshold = 2*median(quality);
    spots_560 = spots_560(quality>threshold, :);
    [k2, dist] = dsearchn(spots_560, spots_405);
    spots_560 = spots_560(k2, :);

    img = img_series_max(:, :, :, 2);
    [spots_488, quality] = log_detector_fft(img, 3, 0, nucleus_mask);
    threshold = 2*median(quality);
    spots_488 = spots_488(quality>threshold, :);
    if with_RNA
        if size(spots_488, 1)>=1
            [k1, dist] = dsearchn(spots_488, spots_405);
            if dist<5 %pixel
                spots_488 = spots_488(k1, :);
            else
                spots_488 = spots_405; % no RNA
            end
        else
            spots_488 = spots_405; % no RNA
        end
    else
        spots_488 = spots_405; % no RNA
    end

    spots_405_2D = spots_405;
    spots_488_2D = spots_488;
    spots_560_2D = spots_560;

    foci_result(1).name = channel_labels{1};
    foci_result(1).spots_3D = spots_405_3D;
    foci_result(1).spots_2D = spots_405_2D;
    foci_result(1).rc_index = [spots_405_2D(2), spots_405_2D(1)]-0.5;
    foci_result(2).name = channel_labels{2};
    foci_result(2).spots_3D = spots_488_3D;
    foci_result(2).spots_2D = spots_488_2D;
    foci_result(2).rc_index = [spots_488_2D(2), spots_488_2D(1)]-0.5;
    foci_result(3).name = channel_labels{3};
    foci_result(3).spots_3D = spots_560_3D;
    foci_result(3).spots_2D = spots_560_2D;
    foci_result(3).rc_index = [spots_560_2D(2), spots_560_2D(1)]-0.5;

    for c_iter = 1:3 % Promoter RNA and SCR channel
    
        disp(['Processing ', channel_labels{c_iter}, ' channel ...']);
        
        img_series = img_series_max(:, :, :, c_iter);
        rc_index = foci_result(c_iter).rc_index;
        
        h = size(img_series, 1); w = size(img_series, 2);
        numberOfPages = size(img_series, 3);

        % calculte foci location
        rc_index_rescale = rc_index*resize_factor;
        refined_bw = logical(imresize(zeros(size(nucleus_mask)), resize_factor));
        for frame_iter = 1:numberOfPages
            refined_bw(round(rc_index_rescale(frame_iter, 1)), round(rc_index_rescale(frame_iter, 2)), frame_iter) = 1;
        end
        % TIFwriter(uint8(imdilate(refined_bw, strel('disk', 5))), [filepath, filename(1:(end-4)), '-', channel_labels{c_iter}, '-Center.tif'], 'lzw');
        
        % calculate foci intensity relative to background
        [rna_bkg,base_bkg] = IntensityCalculation(img_series,nucleus_mask, rc_index, 4, 5);
        norm_base_bkg = movmean(base_bkg, 10)-500;
        intensity = zeros(size(norm_base_bkg));
        for frame_iter = 1:numberOfPages
            % intensity(frame_iter) = rna_bkg(frame_iter, 1)/norm_base_bkg(frame_iter)*norm_base_bkg(1);
            intensity(frame_iter) = rna_bkg(frame_iter, 1)/base_bkg;
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
        img_processed_roi = zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages);
        img_center = logical(zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages));
        temp_refined_bw = imdilate(refined_bw, strel('disk', 3));
        for frame_iter = 1:numberOfPages
            disp(['Processing Frame ', num2str(frame_iter), ' ...']);
            temp_img = imresize(img_series(:, :, frame_iter), resize_factor, 'nearest');
            temp_mask = temp_refined_bw(:, :, frame_iter);
            temp_roi_resize = roi_resize(:, :, frame_iter);
        
            [row1, row2, col1, col2] = getROIboundary(temp_roi_resize, roi_width*resize_factor);
            img_processed_roi(row1:row2, col1:col2, frame_iter) = reshape(temp_img(temp_roi_resize), [row2-row1+1, col2-col1+1]);
            img_center(row1:row2, col1:col2, frame_iter) = reshape(temp_mask(temp_roi_resize), [row2-row1+1, col2-col1+1]);
        end
        TIFwriter(uint16(img_processed_roi), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '.tif']);
        TIFwriter(uint8(img_center), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-Center-roi.tif'], 'lzw');
    
        % important parameter: rc_index, bw, intensity
        foci_result(c_iter).rna_bkg = rna_bkg;
        foci_result(c_iter).base_bkg = base_bkg;
        foci_result(c_iter).intensity = intensity;
        foci_result(c_iter).intensity2 = intensity2;
    end
    
    save([output_path, filename(1:(end-4)), '.mat'], "foci_result", "roi_window", '-append');

%% dealing with condensate channel
condensate_result = struct; 
seg_points = 5; % Threshold for condensate class, above this will be considered as condensates


for c_iter = 4 % OCT4 or BRD4 channel

    img_series = img_series_max(:, :, :, c_iter);

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
        nclust = 7;% Number of intensity clusters;
        seg_point = seg_points;% Threshold for condensate class, above this will be considered as condensates
        temp_img_roi = img_processed_bicubic_roi(:, :, frame_iter);
        [HMRFseg, ~] = HMRFseg4img(temp_img_roi, true(size(temp_img_roi)), nclust, 0.1, 10^(-8));
        bw_HMRF = HMRFseg.img_class>=seg_point;
        local_thresh = min(temp_img_roi(bw_HMRF));
        bw_local = bw_HMRF;

        % local_thresh = local_otsu_rank(img_processed_bicubic_roi(:, :, frame_iter),neighbor_radius*resize_factor, 256);
        % bw_local_otsu = img_processed_bicubic_roi(:, :, frame_iter)>=local_thresh;

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
        % [temp_CDboundary, temp_labels] = getCDboundaryWS(img_processed_bicubic_roi(:, :, frame_iter), bw_local, CDcenter(:, :, frame_iter), CDinterface(:, :, frame_iter), 0);
        [temp_CDboundary, temp_labels] = getCDboundary(img_processed_bicubic_roi(:, :, frame_iter), bw_local, CDcenter(:, :, frame_iter), CDinterface(:, :, frame_iter), 0);
        CDboundary(:, :, frame_iter) = temp_CDboundary;
        labels(:, :, frame_iter) = temp_labels;
        CDmask(:, :, frame_iter) = bw_local;

        temp_img = imresize(img_series(:, :, frame_iter), resize_factor, "nearest");
        img_processed_roi(row1:row2, col1:col2, frame_iter) = reshape(temp_img(temp_roi_resize), [row2-row1+1, col2-col1+1]);

        % export HMRF segmentation
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
    TIFwriter(uint16(img_processed_roi), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-roi.tif']);
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
    dist2bound = zeros(3, numberOfPages);
    dist2center = zeros(3, numberOfPages);
    condensate_radius = zeros(3, numberOfPages);
    equiv_radius = zeros(3, numberOfPages);
    boundary2center = zeros(3, numberOfPages);
    dist2boundary = zeros(3, numberOfPages);
    dist2centroid = zeros(3, numberOfPages);
    boundary2centroid = zeros(3, numberOfPages);
    for f_iter = 1:3 % Promoter RNA and SCR channel
        rc_index = foci_result(f_iter).rc_index*resize_factor;
        for frame_iter = 1:numberOfPages
            [row, col] = find(roi_resize(:, :, frame_iter));
            min_row = min(row); min_col = min(col);
            rc_index(frame_iter, :) = rc_index(frame_iter, :)-[min_row-1, min_col-1]+0.5;
            spots_resize{frame_iter, 1} = spots{frame_iter, 1}*resize_factor-[min_row-1, min_col-1]+0.5;
        end
        % Metric calculations: distance from DNA/RNA to condensate boundary/center
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
fig = figure('Color','w','Position',[200 200 680 680]);
t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile
imagesc(condensate_result(4).img_processed_roi)
axis image off
colormap gray
hold on
visboundaries(condensate_result(4).CDmask, 'Color', 'm', 'LineWidth', 1);
plot(condensate_result(1).spots_resize(2), condensate_result(1).spots_resize(1), 'bo','MarkerSize',7,'LineWidth',1.2);
plot(condensate_result(2).spots_resize(2), condensate_result(2).spots_resize(1), 'ro','MarkerSize',7,'LineWidth',1.2);
plot(condensate_result(3).spots_resize(2), condensate_result(3).spots_resize(1), 'go','MarkerSize',7,'LineWidth',1.2);
spots = condensate_result(4).spots_resize{frame_iter, 1};
plot(spots(:,2), spots(:,1), 'mo','MarkerSize',7,'LineWidth',1.2);
title([condensate_name, ', layer:', num2str(z_layer_used)],'FontSize',11,'FontWeight','bold') %gai cheng used
hold off

print(fig, [output_path, filesep, filename(1:(end-4)), '.png'], '-dpng');

close all;
end
end
end

%% calculate distance from random sites to boundary, centroids and condensate radius
% Analyzation and visualization from pre-processed data start here
% No need to re-run previous steps

roi_width = 31;
resize_factor = 10;
pixelSize = 95;

filepath = '.\Promoter_RNA_SCR_OCT4'; condensate = 'OCT4';
% filepath = '.\Promoter_RNA_SCR_BRD4'; condensate = 'BRD4';

filepath_withRNA = [filepath, 'filter_result_with_RNA_averaged_boundary', filesep];
filepath_withoutRNA = [filepath, 'filter_result_without_RNA_averaged_boundary', filesep];

%Calculate cell in on state
filename_withRNA_list = dir([filepath_withRNA, '*.mat']);
for file_iter = 1:length(filename_withRNA_list)
    filename = filename_withRNA_list(file_iter).name;
    load([filepath_withRNA, filename], "foci_result", "condensate_result", "channel_labels");

    %generate random localizations
    rand_index = [(1 + (roi_width-2)*rand)*resize_factor, (1 + (roi_width-2)*rand)*resize_factor];

    c_iter = 4;
    spots_resize = condensate_result(c_iter).spots_resize;
    CDcenter = condensate_result(c_iter).CDcenter;
    labels = condensate_result(c_iter).labels;
    [dist_2boundary, dist_2centroid, boundary_2centroid] = getDist2InterfaceRadius(rand_index, spots_resize, CDcenter, labels);
    
    condensate_result(c_iter).rand_dist2centroid = min(vecnorm(rand_index-center_list,2,2))*pixelSize/resize_factor;
    condensate_result(c_iter).rand2boundary = dist_2boundary*pixelSize/resize_factor;
    condensate_result(c_iter).rand2centroid = dist_2centroid*pixelSize/resize_factor;
    condensate_result(c_iter).randradius = boundary_2centroid*pixelSize/resize_factor;

    %save in append mode to previous data struct
    save([filepath_withRNA, filename], "condensate_result", '-append');
end

%Calculate cell in off state
filename_woRNA_list = dir([filepath_withoutRNA, '*.mat']);
for file_iter = 1:length(filename_woRNA_list)
    filename = filename_woRNA_list(file_iter).name;
    load([filepath_withoutRNA, filename], "foci_result", "condensate_result", "channel_labels");
    
    rand_index = [(1 + (roi_width-2)*rand)*resize_factor, (1 + (roi_width-2)*rand)*resize_factor];

    c_iter = 4;
    spots_resize = condensate_result(c_iter).spots_resize;
    CDcenter = condensate_result(c_iter).CDcenter;
    labels = condensate_result(c_iter).labels;
    [dist_2boundary, dist_2centroid, boundary_2centroid] = getDist2InterfaceRadius(rand_index, spots_resize, CDcenter, labels);

    condensate_result(c_iter).rand2boundary = dist_2boundary*pixelSize/resize_factor;
    condensate_result(c_iter).rand2centroid = dist_2centroid*pixelSize/resize_factor;
    condensate_result(c_iter).randradius = boundary_2centroid*pixelSize/resize_factor;

    %save in append mode to previous data struct
    save([filepath_withoutRNA, filename], "condensate_result", '-append');
end
%% calculate distance between DNA/RNA sites

filepath = '.\Promoter_RNA_SCR_OCT4'; condensate = 'OCT4';
% filepath = '.\Promoter_RNA_SCR_BRD4'; condensate = 'BRD4';
infotable_path = fullfile(filepath,'OCT4_label_z.xlsx');
% infotable_path = fullfile(filepath,'BRD4_label_z.xlsx');

infotable = readtable(infotable_path);
infotable = table2struct(infotable);
for i = 1:size(infotable,1)
    tempcell = infotable(i).Var1;
    % strrep(filename_list, "'", "")
    tempcell = strrep(tempcell, '''', '');
    infotable(i).name = [tempcell(1:end-3),'mat'];
end

filepath_withRNA = [filepath, 'filter_result_with_RNA_averaged_boundary', filesep];
filepath_withoutRNA = [filepath, 'filter_result_without_RNA_averaged_boundary', filesep];
pixelSize = 95;

%Calculate cell in on state
filename_withRNA_list = dir([filepath_withRNA, '*.mat']);
dist_summary_withRNA = struct();
epdist_rna_on = zeros(length(filename_withRNA_list),3);
for file_iter = 1:length(filename_withRNA_list)
    filename = filename_withRNA_list(file_iter).name;
    load([filepath_withRNA, filename], "foci_result", "condensate_result");

    dist_summary_withRNA(file_iter).filename = filename;

    zstep = getZstep(infotable,filename);

    dist_summary_withRNA(file_iter).prom2cdbound = condensate_result(4).dist2boundary(1);
    dist_summary_withRNA(file_iter).rna2cdbound = condensate_result(4).dist2boundary(2);
    dist_summary_withRNA(file_iter).scr2cdbound = condensate_result(4).dist2boundary(3);
    dist_summary_withRNA(file_iter).rand2cdbound = condensate_result(4).rand2boundary;
    dist_summary_withRNA(file_iter).prom2rna = norm(foci_result(1).rc_index-foci_result(2).rc_index)*pixelSize;
    dist_summary_withRNA(file_iter).scr2rna = norm(foci_result(3).rc_index-foci_result(2).rc_index)*pixelSize;
    dist_summary_withRNA(file_iter).prom2scr = norm(foci_result(1).rc_index-foci_result(3).rc_index)*pixelSize;
    
    prom_3d_pix = [foci_result(1).rc_index, foci_result(1).spots_3D(3)*zstep/pixelSize];
    rna_3d_pix = [foci_result(2).rc_index, foci_result(2).spots_3D(3)*zstep/pixelSize];
    scr_3d_pix = [foci_result(3).rc_index, foci_result(3).spots_3D(3)*zstep/pixelSize];

    dist_summary_withRNA(file_iter).prom2rna3d = norm(prom_3d_pix-rna_3d_pix)*pixelSize;
    dist_summary_withRNA(file_iter).scr2rna3d = norm(scr_3d_pix-rna_3d_pix)*pixelSize;
    dist_summary_withRNA(file_iter).prom2scr3d = norm(prom_3d_pix-scr_3d_pix)*pixelSize;

    epdist_rna_on(file_iter,2) = norm(prom_3d_pix-scr_3d_pix)*pixelSize;

    dist_summary_withRNA(file_iter).rnaintensity = foci_result(2).intensity;

    epdist_rna_on(file_iter,1) = foci_result(2).intensity;

    epdist_rna_on(file_iter,3) = foci_result(2).base_bkg;
end

%Calculate cell in off state
filename_woRNA_list = dir([filepath_withoutRNA, '*.mat']);
dist_summary_woRNA = struct();
epdist_rna_off = zeros(length(filename_withRNA_list),2);
for file_iter = 1:length(filename_woRNA_list)
    filename = filename_woRNA_list(file_iter).name;
    load([filepath_withoutRNA, filename], "foci_result", "condensate_result");

    zstep = getZstep(infotable,filename);

    dist_summary_woRNA(file_iter).filename = filename;

    dist_summary_woRNA(file_iter).prom2cdbound = condensate_result(4).dist2boundary(1);
    dist_summary_woRNA(file_iter).scr2cdbound = condensate_result(4).dist2boundary(3);
    dist_summary_woRNA(file_iter).rand2cdbound = condensate_result(4).rand2boundary;
    dist_summary_woRNA(file_iter).prom2scr = norm(foci_result(1).rc_index-foci_result(3).rc_index)*pixelSize;
    
    prom_3d_pix = [foci_result(1).rc_index, foci_result(1).spots_3D(3)*zstep/pixelSize];
    scr_3d_pix = [foci_result(3).rc_index, foci_result(3).spots_3D(3)*zstep/pixelSize];
    
    dist_summary_woRNA(file_iter).prom2scr3d = norm(prom_3d_pix-scr_3d_pix)*pixelSize;

    epdist_rna_off(file_iter,2) = norm(prom_3d_pix-scr_3d_pix)*pixelSize;

    dist_summary_woRNA(file_iter).rnaintensity = foci_result(2).intensity;

    epdist_rna_off(file_iter,1) = foci_result(2).intensity;
    epdist_rna_off(file_iter,3) = foci_result(2).base_bkg;
end


%% calculate distance from DNA(enhancer)/RNA sites to boundary, centroids and condensate radius

filepath = '.\Promoter_RNA_SCR_OCT4'; condensate = 'OCT4';
% filepath = '.\Promoter_RNA_SCR_BRD4'; condensate = 'BRD4';

filepath_withRNA = [filepath, 'filter_result_with_RNA_averaged_boundary', filesep];
filepath_withoutRNA = [filepath, 'filter_result_without_RNA_averaged_boundary', filesep];

pixelSize = 95;

filename_withRNA_list = dir([filepath_withRNA, '*.mat']);
PRE_centroid_radius = zeros(length(filename_withRNA_list), 9); 
PRE_rand_centroid_radius = zeros(length(filename_withRNA_list), 2);

PRE_each_nearest_centroid_radius = zeros(length(filename_withRNA_list), 4); 
PRE_all_nearest_centroid_radius = zeros(length(filename_withRNA_list), 3);

%Calculate cell in on state
for file_iter = 1:length(filename_withRNA_list)
    filename = filename_withRNA_list(file_iter).name;
    load([filepath_withRNA, filename], "foci_result", "condensate_result", "channel_labels");

    PRE_centroid_radius(file_iter, :) = [condensate_result(4).dist2centroid(1), condensate_result(4).boundary2centroid(1), ...
                                         condensate_result(4).dist2centroid(2), condensate_result(4).boundary2centroid(2), ...
                                         condensate_result(4).dist2centroid(3), condensate_result(4).boundary2centroid(3), ...
                                         norm(foci_result(1).rc_index-foci_result(2).rc_index)*pixelSize, ...
                                         norm(foci_result(1).rc_index-foci_result(3).rc_index)*pixelSize, ...
                                         norm(foci_result(2).rc_index-foci_result(3).rc_index)*pixelSize];
    PRE_rand_centroid_radius(file_iter, :) = [condensate_result(4).rand2centroid, condensate_result(4).randradius];
    
    center_list = condensate_result(4).spots_resize;
    center_list = center_list{1};
    PRE_each_nearest_centroid_radius(file_iter,:) = [min(vecnorm(condensate_result(1).spots_resize-center_list,2,2))*pixelSize/resize_factor, ...
                                                     min(vecnorm(condensate_result(2).spots_resize-center_list,2,2))*pixelSize/resize_factor, ...
                                                     min(vecnorm(condensate_result(3).spots_resize-center_list,2,2))*pixelSize/resize_factor, ...
                                                     condensate_result(4).rand_dist2centroid];
    [d,center_id] = min(vecnorm(mean([condensate_result(1).spots_resize; condensate_result(2).spots_resize; condensate_result(3).spots_resize],1)-center_list,2,2));
    center_loc = center_list(center_id,:);
    PRE_all_nearest_centroid_radius(file_iter,:) = [norm(condensate_result(1).spots_resize-center_loc)*pixelSize/resize_factor, ...
                                                    norm(condensate_result(2).spots_resize-center_loc)*pixelSize/resize_factor, ...
                                                    norm(condensate_result(3).spots_resize-center_loc)*pixelSize/resize_factor];

end

%Calculate cell in off state
filename_woRNA_list = dir([filepath_withoutRNA, '*.mat']);
PE_centroid_radius = zeros(length(filename_woRNA_list), 5); 
PE_rand_centroid_radius = zeros(length(filename_woRNA_list), 2); 
for file_iter = 1:length(filename_woRNA_list)
    filename = filename_woRNA_list(file_iter).name;
    load([filepath_withoutRNA, filename], "foci_result", "condensate_result", "channel_labels");

    PE_centroid_radius(file_iter, :) = [condensate_result(4).dist2centroid(1), condensate_result(4).boundary2centroid(1), ...
                                         condensate_result(4).dist2centroid(3), condensate_result(4).boundary2centroid(3), ...
                                         norm(foci_result(1).rc_index-foci_result(3).rc_index)*pixelSize];
    PE_rand_centroid_radius(file_iter, :) = [condensate_result(4).rand2centroid, condensate_result(4).randradius];
end

%% Scatter plot

range_x_limit = 700;
range_y_limit = 500;

figure('Color','w','Position',[100 100 1100 550])
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% MATLAB default colormap
c_blue  = [0 0.4470 0.7410];
c_red   = [0.8500 0.3250 0.0980];
c_black = [0.15 0.15 0.15];
c_purple = [0.45 0.15 0.60];

% marker parameters
ms = 12;
alpha = 0.35;

% ================= on state row =================

nexttile
plot_scatter_panel(PRE_centroid_radius(:,1), ...
                   PRE_centroid_radius(:,2), ...
                   c_blue, ms, alpha, range_x_limit, range_y_limit)
title(['Promoter(ON) – ', condensate]);

nexttile
plot_scatter_panel(PRE_centroid_radius(:,5), ...
                   PRE_centroid_radius(:,6), ...
                   c_purple, ms, alpha, range_x_limit, range_y_limit)
title(['SCR(ON) – ', condensate]);

nexttile
plot_scatter_panel(PRE_centroid_radius(:,3), ...
                   PRE_centroid_radius(:,4), ...
                   c_red, ms, alpha, range_x_limit, range_y_limit)
title(['RNA(ON) – ', condensate]);

% ================= off state row =================

nexttile
plot_scatter_panel(PE_centroid_radius(:,1), ...
                   PE_centroid_radius(:,2), ...
                   c_blue, ms, alpha, range_x_limit, range_y_limit)
title(['Promoter(OFF) – ', condensate]);

nexttile
plot_scatter_panel(PE_centroid_radius(:,3), ...
                   PE_centroid_radius(:,4), ...
                   c_purple, ms, alpha, range_x_limit, range_y_limit)
title(['SCR(OFF) – ', condensate]);

nexttile
plot_scatter_panel([PE_rand_centroid_radius(:,1); PRE_rand_centroid_radius(:,1)], ...
                   [PE_rand_centroid_radius(:,2); PRE_rand_centroid_radius(:,2)], ...
                   c_black, ms, alpha, range_x_limit, range_y_limit)
title(['Random – ', condensate])

% ========= unified format =========
allAxes = findall(gcf,'Type','axes');
set(allAxes,'Box','off','TickDir','out','LineWidth',1.2,'FontSize',11)
xlabel(t,'DNA/RNA foci to condensate centroid (nm)','FontSize',13,'FontWeight','bold')
ylabel(t,'Condensate boundary to centroid (nm)','FontSize',13,'FontWeight','bold')

set(gcf,'Renderer','painters')

%% Contour plot

range_x_limit = 700;
range_y_limit = 500;

figure('Color','w','Position',[100 100 1100 550])
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% MATLAB default colormap
c_blue  = [0 0.4470 0.7410];
c_red   = [0.8500 0.3250 0.0980];
c_black = [0.15 0.15 0.15];
c_purple = [0.45 0.15 0.60];

nx = 100; ny = 100;
[xg, yg] = meshgrid(linspace(0,range_x_limit,nx), ...
                    linspace(0,range_y_limit,ny));
gridpts = [xg(:) yg(:)];

% ================= on state row =================

nexttile
plot_density_contour(PRE_centroid_radius(:,1), ...
                   PRE_centroid_radius(:,2), ...
                   gridpts, xg, yg, c_blue, range_x_limit, range_y_limit)
title(['Promoter(ON) – ', condensate]);

nexttile
plot_density_contour(PRE_centroid_radius(:,5), ...
                   PRE_centroid_radius(:,6), ...
                   gridpts, xg, yg, c_purple, range_x_limit, range_y_limit)
title(['SCR(ON) – ', condensate]);

nexttile
plot_density_contour(PRE_centroid_radius(:,3), ...
                   PRE_centroid_radius(:,4), ...
                   gridpts, xg, yg, c_red, range_x_limit, range_y_limit)
title(['RNA(ON) – ', condensate]);

% ================= off state row =================

nexttile
plot_density_contour(PE_centroid_radius(:,1), ...
                   PE_centroid_radius(:,2), ...
                   gridpts, xg, yg, c_blue, range_x_limit, range_y_limit)
title(['Promoter(OFF) – ', condensate]);

nexttile
plot_density_contour(PE_centroid_radius(:,3), ...
                   PE_centroid_radius(:,4), ...
                   gridpts, xg, yg, c_purple, range_x_limit, range_y_limit)
title(['SCR(OFF) – ', condensate]);

nexttile
plot_density_contour([PE_rand_centroid_radius(:,1); PRE_rand_centroid_radius(:,1)], ...
                   [PE_rand_centroid_radius(:,2); PRE_rand_centroid_radius(:,2)], ...
                     gridpts, xg, yg, c_black, range_x_limit, range_y_limit)
title(['Random – ', condensate])

% ========= unified format =========
allAxes = findall(gcf,'Type','axes');
set(allAxes,'Box','off','TickDir','out','LineWidth',1.2,'FontSize',11)
xlabel(t,'DNA/RNA foci to condensate centroid (nm)','FontSize',13,'FontWeight','bold')
ylabel(t,'Condensate boundary to centroid (nm)','FontSize',13,'FontWeight','bold')

set(gcf,'Renderer','painters')