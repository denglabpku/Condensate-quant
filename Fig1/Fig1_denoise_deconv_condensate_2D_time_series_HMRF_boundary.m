%% IMAGE PROCESSING PIPELINE FOR BIOMOLECULAR CONDENSATE DYNAMICS
%  Project: Multi-color Analysis of OCT4 & BRD4 Condensates
%  Contact: Bo Wang, Peking University
%  
%  Description: 
%  This script processes multi-channel 3D/time-lapse TIFF images to quantify
%  biomolecular condensate (CD) morphology. The pipeline includes:
%    1. Deep-learning denoising via Noise2Void.
%    2. Richardson-Lucy Deconvolution.
%    3. Nucleus segmentation and HMRF-based condensate identification.
%    4. Morphometric statistics (Area, Radius, Density).

clc;close all;clear;

%% 1. DENOISING & DECONVOLUTION CONFIGURATION
% Pre-trained N2V network for OCT4/BRD4 live-SR data
% onnx_path = '..\onnx\N2V_2D_LiveSR_125_1_E10_xy3z0.onnx';
onnx_path = '../onnx/N2V_2D_LiveSR_125_1_E10_xy3z0.onnx';
is_GPU_avaliable = true;

% Initialize Network and extract patching requirements
net = importNetworkFromONNX(onnx_path);
inputsize = net.Layers(1).InputSize; % Format: [T, H, W, C]
patch_h = inputsize(2);
patch_w = inputsize(3);
patch_t = inputsize(1);
disp('Deep-learning model loaded. Initializing GPU-accelerated denoising...');

% Load Point Spread Function (PSF) for 561nm channel deconvolution
% psf_SR_560 = TIFFreader('..\PSF\psf_SR_channel561_2D.tif', 'double');
psf_SR_560 = TIFFreader('../PSF/psf_SR_channel561_2D.tif', 'double');


%% 2. DATA LOADING & PRE-PROCESSING
filepath_list = {''};
pixelSize = 95;      % Physical pixel size in nanometers (nm)
resize_factor = 10;  % Sub-pixel interpolation factor for morphological precision

for filepath_iter = 1:length(filepath_list)

filepath = filepath_list{filepath_iter};
output_path = filepath;

filename_list = dir([filepath, '*.tif']);

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
    %% 3. NOISE2VOID IMAGE DENOISING
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
    %% 4. RL-DECONVOLUTION
    % Richardson-Lucy algorithm
    pad = 20; % padding to remove edge strip
    img_stack_deconv = double(img_stack_denoised);

    img_denoised = double(img_stack_denoised(:, :, :, 1, :));
    [img_deconv, ~] = deconv_fixed_iter(img_denoised, psf_SR_560, pad, 20);
    img_stack_deconv(:, :, :, 1, :) = img_deconv;

    % export denoised and deconvolved image stacks
    img_stack_denoised = double(squeeze(img_stack_denoised));
    img_stack_deconv = squeeze(img_stack_deconv);
    TIFwriter(uint16(img_stack_deconv), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-denoised-deconv.tif']);
    TIFwriter(uint16(img_stack_denoised), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-denoised.tif']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% nucleus mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 5. HMRF-BASED CONDENSATE SEGMENTATION
    % Generate global nucleus mask
    nucleus_mask = NucleusMask2(max(img_stack_denoised, [], 3), 15);

    z_stack = 1; c_channel = sizeC; 
    h = sizeY; w = sizeX;
    numberOfPages = sizeT;

    save([output_path, filename(1:(end-4)), '.mat'], "img_stack_deconv", "nucleus_mask", "resize_factor", "pixelSize");

%% dealing with condensate channel
condensate_result = struct; 
nclust = 7;% Number of intensity clusters;
seg_point = 5;% Threshold for condensate class, above this will be considered as condensates

img_series = img_stack_deconv;

%%%%%%%%%%%%%%%% get CD center, interface and boundary %%%%%%%%%%%%%%%%

CDcenter = logical(zeros(h*resize_factor, w*resize_factor, numberOfPages));
CDinterface = CDcenter; CDboundary = CDcenter; CDmask = CDboundary;
labels = uint16(CDcenter);
area = cell(numberOfPages, 1);

nucleus_mask_resize = imresize(nucleus_mask, resize_factor, "nearest");
for frame_iter = 1:numberOfPages
    disp(['Processing Frame ', num2str(frame_iter), ' ...']);
    temp_img = imresize(img_series(:, :, frame_iter), resize_factor);

    % Hidden Markov Random Field (HMRF) to include spatial dependency
    temp_resize_factor = 2;
    temp_img_for_HMRF = imresize(img_series(:, :, frame_iter), temp_resize_factor);
    [HMRFseg, ~] = HMRFseg4img(temp_img_for_HMRF, imresize(nucleus_mask, temp_resize_factor, "nearest"), nclust, 0.1, 10^(-8));
    bw_HMRF = HMRFseg.img_class>=seg_point;
    local_thresh = min(temp_img_for_HMRF(bw_HMRF));

    bw_local = temp_img>local_thresh & nucleus_mask_resize; 

    CDcenter(:, :, frame_iter) = getCDcenter(img_series(:, :, frame_iter), bw_local, resize_factor, 0);
    CDinterface(:, :, frame_iter) = getCDinterface(temp_img, bw_local, 0);
    [temp_CDboundary, temp_labels] = getCDboundary(temp_img, bw_local, CDcenter(:, :, frame_iter), CDinterface(:, :, frame_iter), 0);
    CDboundary(:, :, frame_iter) = temp_CDboundary;
    labels(:, :, frame_iter) = temp_labels;
    CDmask(:, :, frame_iter) = bw_local;

    %%%%%%%%%%%%%%%%%%% calculate the condensate area %%%%%%%%%%%%%%%%%%%%%

    [C,ia,ic] = unique(temp_labels(:));
    a_counts = accumarray(ic,1);
    value_counts = [C, a_counts];
    value_counts = value_counts(C>0, :);
    area{frame_iter, 1} = value_counts(:, 2);

    % export HMRF segmentation
    if frame_iter == 1
        fig1 = figure;
        fig1.Units = "inches";
        fig1.Position = [7.4,3.9,9.8,6.4];
        subplot(1, 2, 1);
        imagesc(temp_img_for_HMRF);
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
        print(fig1, [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-HMRFseg.png'], '-dpng');
        close;
    end
end
CDmask_export = uint8(CDmask)*50;
CDmask_export(CDboundary)=255;

% Export morphological results as TIFF
TIFwriter(uint8(CDcenter), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-CDcenter.tif'], 'lzw');
TIFwriter(uint8(CDinterface), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-CDinterface.tif'], 'lzw');
TIFwriter(uint8(CDboundary), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-CDboundary.tif'], 'lzw');
TIFwriter(uint8(CDmask_export), [output_path, filename(1:(end-4)), filesep, filename(1:(end-4)), '-CDmask.tif']);

% Store results in struct
condensate_result.CDcenter = CDcenter;
condensate_result.CDinterface = CDinterface;
condensate_result.CDboundary = CDboundary;
condensate_result.CDmask = CDmask;
condensate_result.labels = labels;
condensate_result.area = area;

% important parameter: condensate_result
save([output_path, filename(1:(end-4)), '.mat'], "condensate_result", '-append');

end

close all;

end

%% Calculate morphometric statistics (Area, Radius, Density)
% Analyzation and visualization from pre-processed data start here
% No need to re-run previous steps

filepath = '';
output_path = filepath;

condensate_statistics = struct;
condensate_label = {'BRD4', 'OCT4'};

for filepath_iter = 1:2

    condensate_name = condensate_label{filepath_iter};

    filename_list = dir([filepath, '*', condensate_name, '*.mat']);
    condensate_statistics(filepath_iter).name = condensate_name;

    pixelSize = 95; %nm
    resize_factor = 10;
    numberOfPages = 61;
    
    radius_thresh = 0;
    area_cutoff = pi*radius_thresh^2/(pixelSize/resize_factor)^2;
    
    condensate_num_over_time = zeros(numberOfPages, length(filename_list));
    condensate_area_over_time = cell(numberOfPages, length(filename_list));
    mean_condensate_area_over_time = zeros(numberOfPages, length(filename_list));
    area_in_total = [];
    
    for file_iter = 1:length(filename_list)
    
        filename = filename_list(file_iter).name;
        disp(['Processing ', filename, ' ...']);
        load([filepath, filename], "condensate_result");
    
        for frame_iter = 1:numberOfPages
    
            temp_area = condensate_result.area{frame_iter, 1};
            temp_area = temp_area(temp_area>=area_cutoff);
            condensate_num_over_time(frame_iter, file_iter) = length(temp_area);
            condensate_area_over_time{frame_iter, file_iter} = temp_area*(pixelSize/resize_factor)^2;
            if frame_iter == 20
                area_in_total = [area_in_total; temp_area*(pixelSize/resize_factor)^2];
            end
            mean_condensate_area_over_time(frame_iter, file_iter) = mean(temp_area)*(pixelSize/resize_factor)^2;
    
        end
    
    end
    
    %save in struct
    condensate_statistics(filepath_iter).num_over_time = condensate_num_over_time;
    condensate_statistics(filepath_iter).mean_num = mean(condensate_statistics(filepath_iter).num_over_time, 1);
    condensate_statistics(filepath_iter).area_over_time = condensate_area_over_time;
    condensate_statistics(filepath_iter).mean_area_over_time = mean_condensate_area_over_time;
    condensate_statistics(filepath_iter).area_in_total = area_in_total;
    condensate_statistics(filepath_iter).radius = sqrt(condensate_statistics(filepath_iter).area_in_total/pi);

end

save([filepath, 'condensate_statistics.mat'], "condensate_statistics");


%% Visualization of morphometric statistics

figure;
shadedErrorBar(1:numberOfPages, mean(condensate_statistics(1).num_over_time, 2), std(condensate_statistics(1).num_over_time, [], 2), ...
    'lineProps',{'-', 'color', '#0072BD', 'linewidth', 2});
hold on
shadedErrorBar(1:numberOfPages, mean(condensate_statistics(2).num_over_time, 2), std(condensate_statistics(2).num_over_time, [], 2), ...
    'lineProps',{'-', 'color', '#D95319', 'linewidth', 2});
hold off
xlim([1, numberOfPages]);
ylim([0, 500]);
ylabel('average CD number per cell');

figure;
shadedErrorBar(1:numberOfPages, mean(sqrt(condensate_statistics(1).mean_area_over_time/pi), 2), std(sqrt(condensate_statistics(1).mean_area_over_time/pi), [], 2), ...
    'lineProps',{'-', 'color', '#0072BD', 'linewidth', 2});
hold on
shadedErrorBar(1:numberOfPages, mean(sqrt(condensate_statistics(2).mean_area_over_time/pi), 2), std(sqrt(condensate_statistics(2).mean_area_over_time/pi), [], 2), ...
    'lineProps',{'-', 'color', '#D95319', 'linewidth', 2});
hold off
xlim([1, numberOfPages]);
ylabel('average CD radius per cell');
ylim([0, 300])

%% Loading SIM Data

oct4_cell_number = [];
oct4_cell_radius = [];
brd4_cell_number = [];
brd4_cell_radius = [];

for iter = 1:22

    oct4_cell_number = [oct4_cell_number; length(merge_tab_stats3D{1, iter})];
    oct4_cell_radius = [oct4_cell_radius; merge_tab_stats3D{1, iter}];

    brd4_cell_number = [brd4_cell_number; length(merge_tab_stats3D{2, iter})];
    brd4_cell_radius = [brd4_cell_radius; merge_tab_stats3D{2, iter}];

end
