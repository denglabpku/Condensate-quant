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
script_dir = pwd;
repo_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(repo_root, 'Function')));

% Pre-trained N2V network for OCT4/BRD4 live-SR data
onnx_path = fullfile(repo_root, 'onnx', 'N2V_2D_OCT4_BRD4_liveSR_125_1_E10_xy3z0.onnx');
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
psf_SR_405 = TIFFreader(fullfile(repo_root, 'PSF', 'psf_SR_channel405_2D.tif'), 'double');
psf_SR_488 = TIFFreader(fullfile(repo_root, 'PSF', 'psf_SR_channel488_2D.tif'), 'double');
psf_SR_560 = TIFFreader(fullfile(repo_root, 'PSF', 'psf_SR_channel561_2D.tif'), 'double');

%% data loading and pre-processing
% These folders map to the promoter-RNA-SCR datasets with either BRD4 or OCT4 condensate labeling.

filepath_list = {fullfile('D:\ImageData\Supplementary_Imaging_Data\ExampleData', 'Promoter_RNA_SCR')};
channel_labels = {'Promoter', 'RNA', 'Enhancer'};

for filepath_iter = 1:length(filepath_list)

filepath = filepath_list{filepath_iter};

for with_RNA = [1, 0]

if with_RNA
    output_path = [filepath, filesep, 'filter_result_with_RNA_averaged_boundary_intensity', filesep];
    mkdir(output_path)
    
    filename_list = readlines([filepath, filesep, 'Cell_with_RNA.txt']);
    filename_list = strrep(filename_list, "'", "");
else
    output_path = [filepath, filesep, 'filter_result_without_RNA_averaged_boundary_intensity', filesep];
    mkdir(output_path)
    
    filename_list = readlines([filepath, filesep, 'Cell_without_RNA.txt']);
    filename_list = strrep(filename_list, "'", "");
end

filename_list = filename_list(strlength(filename_list) > 0); % remove empty string
error_log_path = fullfile(output_path, 'Cell_analysis_errors.txt');
initializeCellProcessingErrorLog(error_log_path);

pixelSize = 95; %nm
roi_width = 31;

for file_iter = 1:length(filename_list)

    filename = filename_list{file_iter};
    filename = [filename(1:end-3),'tif'];
    mkdir([output_path, filename(1:(end-4))]);
    disp(['Processing ', filename, ' ...']);
    r = [];
    try

    %%%%%%%%%%%%%%%%%%%%%%%%%% read img sequence %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    r = bfGetReader([filepath, filesep, filename]);

    % r = bfGetReader([file_path, filename]);


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
    img_stack_denoised = img_stack;

    if sizeY<patch_h || sizeX<patch_w
        warning('ImageTooSmall:SizeCheck', ...
            'The input image is too small (%dx%d). Minimum required is %dx%d.', ...
            sizeY, sizeX, patch_h, patch_w);
        continue
    end
    overlap_factor=0.25;
    gap_h = round(patch_h*(1-overlap_factor));% Patch gap in height
    gap_w = round(patch_w*(1-overlap_factor));% Patch gap in width
    gap_t = round(patch_t*(1-overlap_factor));% Patch gap in slice

    for c = 1:min(sizeC, 3)

        img = single(img_stack(:, :, :, c)); img_mean = mean(img(:));
        img_permute = permute(img, [3, 1, 2])-img_mean;
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
            denoised_patches{patch_iter} = predict(net, I_dlarray);
        end
        close(h); 
        
        img_denoised_permute = reconstructFromPatches(denoised_patches, coordinates, [patch_t, patch_h, patch_w], size(img_permute), [gap_t, gap_h, gap_w]);
        img_denoised = permute(img_denoised_permute, [2, 3, 1])+img_mean;
        img_stack_denoised(:, :, :, c) = gather(img_denoised);

    end
    disp('Image denoising complete!');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% deconvolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    % export denoised and deconvolved image stacks
    bfsave(uint16(img_stack_deconv), fullfile(output_path, [filename(1:(end-4)), '-denoised-deconv.ome.tif']));
    bfsave(uint16(img_stack_reconv), fullfile(output_path, [filename(1:(end-4)), '-denoised.ome.tif']));

    img_stack_deconv(:, :, :, 1) = double(img_stack_denoised(:, :, :, 1));
    img_stack_deconv(:, :, :, 2) = double(img_stack_denoised(:, :, :, 2));
    img_stack_deconv(:, :, :, 3) = double(img_stack_denoised(:, :, :, 3));

    channel_405 = img_stack_deconv(:, :, :, 1);
    channel_488 = img_stack_deconv(:, :, :, 2);
    channel_560 = img_stack_deconv(:, :, 1:(sizeZ-1), 3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% nucleus mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % nucleus_mask = NucleusMask2(mean(img_stack_reconv(:, :, layer_select, 4), 3), 15);
    nucleus_mask = NucleusMask2(mean(img_stack_reconv(:, :, layer_select, 1), 3), 15);

    temp_channel_405 = channel_405(:, :, layer_select);
    [spots_405, quality, I_log, err_xy_405, err_xyz_405, err_2D_405] = log_detector_3d_Sage_v3(temp_channel_405, 3, 3, mean(channel_405(:))/5);
    
    
    id = spots_405(:, 1)>5 & spots_405(:, 1)<(sizeX-5);
    spots_405 = spots_405(id, :); quality = quality(id);
    err_xy_405 = err_xy_405(id, :);
    err_xyz_405 = err_xyz_405(id, :);
    err_2D_405 = err_2D_405(id,:);
    [~, i] = max(quality); spots_405 = spots_405(i, :);
    err_xy_405 = err_xy_405(i, :);
    err_xyz_405 = err_xyz_405(i, :);
    err_2D_405 = err_2D_405(i, :);
    spots_405(3) = spots_405(3)+layer_select(1)-1;

    [spots_560, quality, I_log, err_xy_560, err_xyz_560, err_2D_560] = log_detector_3d_Sage_v3(channel_560, 3, 3, 0);
    threshold = 2*median(quality);
    spots_560 = spots_560(quality>threshold, :);
    err_xy_560 = err_xy_560(quality>threshold, :);
    err_xyz_560 = err_xyz_560(quality>threshold, :);
    err_2D_560 = err_2D_560(quality>threshold, :);
    [k2, dist] = dsearchn(spots_560, spots_405);
    spots_560 = spots_560(k2, :);
    err_xy_560 = err_xy_560(k2,:);
    err_xyz_560 = err_xyz_560(k2,:);
    err_2D_560 = err_2D_560(k2,:);

    if with_RNA
        [spots_488, quality, I_log, err_xy_488, err_xyz_488, err_2D_488] = log_detector_3d_Sage_v3(channel_488, 3, 3, 0);
        threshold = 2*median(quality);
        spots_488 = spots_488(quality>threshold, :);
        err_xy_488 = err_xy_488(quality>threshold, :);
        err_xyz_488 = err_xyz_488(quality>threshold, :);
        err_2D_488 = err_2D_488(quality>threshold, :);

        if size(spots_488, 1)>=1
            [k1, dist] = dsearchn(spots_488(:, 1:2), spots_405(:, 1:2));
            
            if dist<5 %pixel
                spots_488 = spots_488(k1, :);
                err_xy_488 = err_xy_488(k1, :);
                err_xyz_488 = err_xyz_488(k1, :);
                err_2D_488 = err_2D_488(k1, :);
            else
                spots_488 = spots_405; % no RNA
            end
        else
            spots_488 = spots_405; % no RNA
        end
    else
        spots_488 = spots_405; % no RNA
    end

    z_layer_405 = max(min(round(spots_405(3)),size(channel_405, 3)),1);
    z_layer_488 = max(min(round(spots_488(3)),size(channel_488, 3)),1);
    z_layer_560 = max(min(round(spots_560(3)),size(channel_560, 3)),1);

    z_layer_used = max(min([z_layer_405, z_layer_488, z_layer_560]), layer_select(1)):min(max([z_layer_405, z_layer_488, z_layer_560]), sizeZ-1);
    
    if isempty(z_layer_used)
        if z_layer_405 == sizeZ
            z_layer_used = sizeZ-1;
        elseif z_layer_405 < layer_select(1)
            z_layer_used = layer_select(1);
        end
    end
    
    img_series_max = zeros(sizeY, sizeX, 1, min(sizeC, 3)); 
    img_series_max(:, :, :, 1) = channel_405(:, :, z_layer_405);
    img_series_max(:, :, :, 2) = channel_488(:, :, z_layer_488);
    img_series_max(:, :, :, 3) = channel_560(:, :, z_layer_560);

    z_stack = 1; c_channel = sizeC; 
    h = sizeY; w = sizeX;
    numberOfPages = sizeT;

    resize_factor = 10;
    spots_405_3D = spots_405;
    spots_488_3D = spots_488;
    spots_560_3D = spots_560;

    save([output_path, filename(1:(end-4)), '.mat'], "img_series_max", "nucleus_mask", "resize_factor", "channel_labels", "pixelSize");

    %% dealing with foci channel
    foci_result = struct();

    img = img_series_max(:, :, :, 1);
    [spots_405, quality] = log_detector_fft(img, 3, 0, nucleus_mask);
    id = abs(spots_405(:, 1)-spots_405_3D(1))<3 & abs(spots_405(:, 2)-spots_405_3D(2))<3;
    spots_405 = spots_405(id, :); quality = quality(id);
    [~, i] = max(quality); spots_405 = spots_405(i, :);

    if isempty(spots_405)
        continue
    end

    img = img_series_max(:, :, :, 3);
    [spots_560, quality] = log_detector_fft(img, 3, 0, nucleus_mask);
    threshold = 2*median(quality);
    spots_560 = spots_560(quality>threshold, :);
    [k2, dist] = dsearchn(spots_560, spots_405);
    spots_560 = spots_560(k2, :);

    if isempty(spots_560)
        continue
    end

    img = img_series_max(:, :, :, 2);
    [spots_488, quality] = log_detector_fft(img, 3, 0, nucleus_mask);
    threshold = 2*median(quality);
    spots_488 = spots_488(quality>threshold, :);

    if isempty(spots_488)
        continue
    end

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
    foci_result(1).err_xy = err_xy_405;
    foci_result(1).err_xyz = err_xyz_405;
    foci_result(1).err_2D = err_2D_405;

    foci_result(2).name = channel_labels{2};
    foci_result(2).spots_3D = spots_488_3D;
    foci_result(2).spots_2D = spots_488_2D;
    foci_result(2).rc_index = [spots_488_2D(2), spots_488_2D(1)]-0.5;
    foci_result(2).err_xy = err_xy_488;
    foci_result(2).err_xyz = err_xyz_488;
    foci_result(2).err_2D = err_2D_488;

    foci_result(3).name = channel_labels{3};
    foci_result(3).spots_3D = spots_560_3D;
    foci_result(3).spots_2D = spots_560_2D;
    foci_result(3).rc_index = [spots_560_2D(2), spots_560_2D(1)]-0.5;
    foci_result(3).err_xy = err_xy_560;
    foci_result(3).err_xyz = err_xyz_560;
    foci_result(3).err_2D = err_2D_560;

    for c_iter = 1:3 % DNA and RNA channel
    
        disp(['Processing ', channel_labels{c_iter}, ' channel ...']);
        
        img_series = img_series_max(:, :, :, c_iter);
        rc_index = foci_result(c_iter).rc_index;
        
        h = size(img_series, 1); w = size(img_series, 2);
        numberOfPages = size(img_series, 3);

        % calculte intensity
        rc_index_rescale = rc_index*resize_factor;
        refined_bw = logical(imresize(zeros(size(nucleus_mask)), resize_factor));
        for frame_iter = 1:numberOfPages
            refined_bw(round(rc_index_rescale(frame_iter, 1)), round(rc_index_rescale(frame_iter, 2)), frame_iter) = 1;
        end
    
        [rna_bkg,base_bkg] = IntensityCalculation(img_series,nucleus_mask, rc_index, 4, 5);
        norm_base_bkg = movmean(base_bkg, 10)-500;
        intensity = zeros(size(norm_base_bkg));
        for frame_iter = 1:numberOfPages
            intensity(frame_iter) = rna_bkg(frame_iter, 1)/base_bkg;
        end
        
        [rna_bkg2,base_bkg2] = IntensityCalculation(img_series,nucleus_mask, rc_index, 3, 4);
        norm_base_bkg2 = movmean(base_bkg2, 10)-500;
        intensity2 = zeros(size(norm_base_bkg2));
        for frame_iter = 1:numberOfPages
            intensity2(frame_iter) = rna_bkg2(frame_iter, 1)/norm_base_bkg2(frame_iter)*norm_base_bkg2(1);
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

    close all;
    catch ME
        cleanupCellProcessingResources(r);
        logCellProcessingError(error_log_path, filename, ME);
        warning('CellProcessing:Failed', ...
            'Failed processing %s. Skipping to next cell. See %s.', filename, error_log_path);
        continue
    end
end
end
end

%% Calculate distances between DNA, RNA and enhancer sites.
% This block processes one unified folder directly and uses the original getZstep function.
% Required file:
%   zstep_lookup_all_files.csv
% Required columns in zstep_lookup_all_files.csv:
%   name, Var2
% where name is the .mat filename and Var2 is the z-step in nm.

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

filepath = fullfile(script_dir, 'rawdata');

pixelSize = 95; % nm

zstep_csv_path = fullfile(filepath, 'zstep_sheet.csv');
zstep_table = readtable(zstep_csv_path);
zstep_info = table2struct(zstep_table);

for i = 1:numel(zstep_info)
    zstep_info(i).name = char(string(zstep_info(i).name));
end

filepath_withRNA = fullfile(filepath, 'filter_result_with_RNA_averaged_boundary_intensity');
filepath_withoutRNA = fullfile(filepath, 'filter_result_without_RNA_averaged_boundary_intensity');

dist_summary_withRNA = struct();
dist_summary_woRNA = struct();

idx_withRNA = 1;
idx_woRNA = 1;

epdist_rna_on = [];
epdist_rna_off = [];

if exist(filepath_withRNA, 'dir')
    filename_withRNA_list = dir(fullfile(filepath_withRNA, '*.mat'));

    fprintf('Processing with-RNA group: %d files\n', numel(filename_withRNA_list));

    for file_iter = 1:length(filename_withRNA_list)

        filename = filename_withRNA_list(file_iter).name;

        data = load(fullfile(filepath_withRNA, filename), "foci_result");

        if ~isfield(data, 'foci_result')
            warning('foci_result not found in file: %s', filename);
            continue;
        end

        foci_result = data.foci_result;

        zstep = getZstep(zstep_info, filename);

        dist_summary_withRNA(idx_withRNA).filename = filename;
        dist_summary_withRNA(idx_withRNA).zstep_nm = zstep;

        prom_xy = foci_result(1).spots_3D(1:2);
        rna_xy  = foci_result(2).spots_3D(1:2);
        scr_xy  = foci_result(3).spots_3D(1:2);

        dist_summary_withRNA(idx_withRNA).prom2rna = norm(prom_xy - rna_xy) * pixelSize;
        dist_summary_withRNA(idx_withRNA).scr2rna  = norm(scr_xy  - rna_xy) * pixelSize;
        dist_summary_withRNA(idx_withRNA).prom2scr = norm(prom_xy - scr_xy) * pixelSize;

        dist_summary_withRNA(idx_withRNA).dz12 = abs(foci_result(1).spots_3D(3) - foci_result(2).spots_3D(3));
        dist_summary_withRNA(idx_withRNA).dz13 = abs(foci_result(1).spots_3D(3) - foci_result(3).spots_3D(3));
        dist_summary_withRNA(idx_withRNA).dz23 = abs(foci_result(2).spots_3D(3) - foci_result(3).spots_3D(3));

        prom_3d_pix = [foci_result(1).rc_index, foci_result(1).spots_3D(3) * zstep / pixelSize];
        rna_3d_pix  = [foci_result(2).rc_index, foci_result(2).spots_3D(3) * zstep / pixelSize];
        scr_3d_pix  = [foci_result(3).rc_index, foci_result(3).spots_3D(3) * zstep / pixelSize];

        dist_summary_withRNA(idx_withRNA).prom2rna3d = norm(prom_3d_pix - rna_3d_pix) * pixelSize;
        dist_summary_withRNA(idx_withRNA).scr2rna3d  = norm(scr_3d_pix  - rna_3d_pix) * pixelSize;
        dist_summary_withRNA(idx_withRNA).prom2scr3d = norm(prom_3d_pix - scr_3d_pix) * pixelSize;

        dist_summary_withRNA(idx_withRNA).err_xy_prom = foci_result(1).err_xy * pixelSize;
        dist_summary_withRNA(idx_withRNA).err_xyz_prom = sqrt(((foci_result(1).err_xyz)^2 - (foci_result(1).err_xy)^2) * zstep^2 + ...
                                                               (foci_result(1).err_xy)^2 * pixelSize^2);
        dist_summary_withRNA(idx_withRNA).err_2D_prom = foci_result(1).err_2D * pixelSize;

        dist_summary_withRNA(idx_withRNA).err_xy_rna = foci_result(2).err_xy * pixelSize;
        dist_summary_withRNA(idx_withRNA).err_xyz_rna = sqrt(((foci_result(2).err_xyz)^2 - (foci_result(2).err_xy)^2) * zstep^2 + ...
                                                              (foci_result(2).err_xy)^2 * pixelSize^2);
        dist_summary_withRNA(idx_withRNA).err_2D_rna = foci_result(2).err_2D * pixelSize;

        dist_summary_withRNA(idx_withRNA).err_xy_scr = foci_result(3).err_xy * pixelSize;
        dist_summary_withRNA(idx_withRNA).err_xyz_scr = sqrt(((foci_result(3).err_xyz)^2 - (foci_result(3).err_xy)^2) * zstep^2 + ...
                                                              (foci_result(3).err_xy)^2 * pixelSize^2);
        dist_summary_withRNA(idx_withRNA).err_2D_scr = foci_result(3).err_2D * pixelSize;

        dist_summary_withRNA(idx_withRNA).rnaintensity = foci_result(2).intensity;

        epdist_rna_on(idx_withRNA, 1) = foci_result(2).intensity;
        epdist_rna_on(idx_withRNA, 2) = norm(prom_3d_pix - scr_3d_pix) * pixelSize;
        epdist_rna_on(idx_withRNA, 3) = foci_result(2).base_bkg;

        idx_withRNA = idx_withRNA + 1;
    end
else
    warning('with-RNA result folder not found: %s', filepath_withRNA);
end

if exist(filepath_withoutRNA, 'dir')
    filename_woRNA_list = dir(fullfile(filepath_withoutRNA, '*.mat'));

    fprintf('Processing without-RNA group: %d files\n', numel(filename_woRNA_list));

    for file_iter = 1:length(filename_woRNA_list)

        filename = filename_woRNA_list(file_iter).name;

        data = load(fullfile(filepath_withoutRNA, filename), "foci_result");

        if ~isfield(data, 'foci_result')
            warning('foci_result not found in file: %s', filename);
            continue;
        end

        foci_result = data.foci_result;

        zstep = getZstep(zstep_info, filename);

        dist_summary_woRNA(idx_woRNA).filename = filename;
        dist_summary_woRNA(idx_woRNA).zstep_nm = zstep;

        prom_xy = foci_result(1).spots_3D(1:2);
        scr_xy  = foci_result(3).spots_3D(1:2);

        dist_summary_woRNA(idx_woRNA).prom2scr = norm(prom_xy - scr_xy) * pixelSize;

        prom_3d_pix = [foci_result(1).rc_index, foci_result(1).spots_3D(3) * zstep / pixelSize];
        scr_3d_pix  = [foci_result(3).rc_index, foci_result(3).spots_3D(3) * zstep / pixelSize];

        dist_summary_woRNA(idx_woRNA).prom2scr3d = norm(prom_3d_pix - scr_3d_pix) * pixelSize;

        dist_summary_woRNA(idx_woRNA).err_xy_prom = foci_result(1).err_xy * pixelSize;
        dist_summary_woRNA(idx_woRNA).err_xyz_prom = sqrt(((foci_result(1).err_xyz)^2 - (foci_result(1).err_xy)^2) * zstep^2 + ...
                                                            (foci_result(1).err_xy)^2 * pixelSize^2);
        dist_summary_woRNA(idx_woRNA).err_2D_prom = foci_result(1).err_2D * pixelSize;

        dist_summary_woRNA(idx_woRNA).err_xy_scr = foci_result(3).err_xy * pixelSize;
        dist_summary_woRNA(idx_woRNA).err_xyz_scr = sqrt(((foci_result(3).err_xyz)^2 - (foci_result(3).err_xy)^2) * zstep^2 + ...
                                                           (foci_result(3).err_xy)^2 * pixelSize^2);
        dist_summary_woRNA(idx_woRNA).err_2D_scr = foci_result(3).err_2D * pixelSize;

        dist_summary_woRNA(idx_woRNA).rnaintensity = foci_result(2).intensity;

        epdist_rna_off(idx_woRNA, 1) = foci_result(2).intensity;
        epdist_rna_off(idx_woRNA, 2) = norm(prom_3d_pix - scr_3d_pix) * pixelSize;
        epdist_rna_off(idx_woRNA, 3) = foci_result(2).base_bkg;

        idx_woRNA = idx_woRNA + 1;
    end
else
    warning('without-RNA result folder not found: %s', filepath_withoutRNA);
end

fprintf('\nMerge complete!\n');
fprintf('with-RNA group: %d cells.\n', idx_withRNA - 1);
fprintf('without-RNA group: %d cells.\n', idx_woRNA - 1);