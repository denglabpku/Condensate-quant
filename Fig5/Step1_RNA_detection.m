%  Step1_RNA_detection.m
%  Bo Wang, July 2025
close all; clear; clc;

%  DESCRIPTION
%  This script processes time-lapse RNA bursting data. It takes a .tif file 
%  containing the RNA signal and an Ilastik-generated probability map (.h5 file) 
%  as inputs, and outputs the detected transcriptional spots and their tracks 
%  over time.
%%

% define input image sequence filename
filepath_list = {'D:/ImageData/Untreated/max_projection/20241112_B27E24-NC-Z-MULTI-ABSOLUTEZ-30S-006-01-Max.tif', ...
                 'D:/ImageData/JQ1-treat/max_projection/20241112_B27E27-JQ1-1uM-15min-Z-MULTI-ABSOLUTEZ-30S-005-15-Max.tif'};
 
for cell_iter = 1:length(filepath_list)

    filepath = filepath_list{cell_iter};
    
    % read Ilastik-generated probability map (.h5 file) 
    dataset = h5read([filepath(1:(end-4)), '_Probabilities.h5'], '/exported_data');
    % read time-lapse RNA bursting data （.tif file）
    img_series = TIFFreader(filepath);
    
    h = size(img_series, 1); w = size(img_series, 2); numberOfPages = size(img_series, 3);
    
    reshaped_img = reshape(dataset(1, :, :, :), [w, h, numberOfPages]);
    
    threshold = 0.5;
    bw = reshaped_img>=threshold;
    
    temp = logical(zeros(h, w, numberOfPages));
    for i = 1:numberOfPages
        temp(:, :, i) = bw(:, :, i)';
    end
    bw = temp;
    
    roitable = readtable([filepath(1:(end-4)), '-ROI.txt']);
    roilist = zeros(size(roitable, 1)/4, 4); % xmin, ymin, xmax, ymax
    for i = 1:size(roitable, 1)/4
        roilist(i, 1:2)=roitable{4*(i-1)+1, 1:2};
        roilist(i, 3:4)=roitable{4*(i-1)+3, 1:2};
    end
    
    spots = cell(0);
    for frame_iter = 1:numberOfPages
        img = img_series(:, :, frame_iter);
        ccomp = bwconncomp(bw(:, :, frame_iter));
    
        locs_intensity = zeros(ccomp.NumObjects, 3);
        for i = 1:ccomp.NumObjects
            [m, idx] = max(img(ccomp.PixelIdxList{1, i}));
            locs = ccomp.PixelIdxList{1, i}(idx);
            locs_intensity(i, 3) = m;
            [a, b] = ind2sub(ccomp.ImageSize,locs);
            locs_intensity(i, 1:2) = [a, b];
        end
    
        spots{frame_iter} = locs_intensity;
    end
    
    tracks = cell(0);
    for track_iter = 1:size(roilist, 1)
        xmin = roilist(track_iter, 1); ymin = roilist(track_iter, 2);
        xmax = roilist(track_iter, 3); ymax = roilist(track_iter, 4);
    
        temp = zeros(numberOfPages, 4);
        for frame_iter = 1:numberOfPages
            idx = find((xmin<spots{frame_iter}(:, 2))&(spots{frame_iter}(:, 2)<=xmax)&(ymin<spots{frame_iter}(:, 1))&(spots{frame_iter}(:, 1)<=ymax));
            temp(frame_iter, 1) = frame_iter;
            if ~isempty(idx)
                res = spots{frame_iter}(idx, :);
                [m, i] = max(res(:, 3));
                temp(frame_iter, 2:4) = res(i, :);
            end
        end
        tracks{track_iter, 1} = temp;
    
        % fill zeros
        temp_i = temp(:, 4);
        % begin zeros
        non_zero_indices = find(temp_i ~= 0);
        if temp(1, 4) == 0
            temp(1:(non_zero_indices(1)-1), 2) = temp(non_zero_indices(1), 2);
            temp(1:(non_zero_indices(1)-1), 3) = temp(non_zero_indices(1), 3);
            for i = 1:(non_zero_indices(1)-1)
                temp(i, 4) = img_series(temp(i, 2), temp(i, 3), i);
            end
        end
        % end zeros
        if temp(end, 4) == 0
            temp((non_zero_indices(end)+1):end, 2) = temp(non_zero_indices(end), 2);
            temp((non_zero_indices(end)+1):end, 3) = temp(non_zero_indices(end), 3);
            for i = (non_zero_indices(end)+1):numberOfPages
                temp(i, 4) = img_series(temp(i, 2), temp(i, 3), i);
            end
        end
        % middle zeros
        % linear padding
        temp(temp(:, 2) == 0, 2) = nan;
        temp(temp(:, 3) == 0, 3) = nan;
        
        temp = fillmissing(temp, 'linear', 1);
        temp = round(temp);
        for frame_iter = 1:numberOfPages
            if temp(frame_iter, 4)==0
                temp(frame_iter, 4) = img_series(temp(frame_iter, 2), temp(frame_iter, 3), frame_iter);
            end
        end
    
        tracks{track_iter, 2} = temp; % frame, x, y, intensity
    end
    
     % calculate RNA intensity: background(7*7), RNA(5*5)
     rna_width = 2; bkg_width = 3;
     for track_iter = 1:size(tracks, 1)
        temp = tracks{track_iter, 2};
        rna_bkg = zeros(numberOfPages, 2);
        rna_bkg_cor = zeros(numberOfPages, 2);
        for frame_iter = 1:numberOfPages
            temp_img = img_series(:, :, frame_iter);
            rna_mask = logical(zeros(size(temp_img))); bkg_mask = rna_mask;
            rna_mask(max(temp(frame_iter, 2)-rna_width, 1):min(temp(frame_iter, 2)+rna_width, h), max(temp(frame_iter, 3)-rna_width, 1):min(temp(frame_iter, 3)+rna_width, w)) = 1;
            bkg_mask(max(temp(frame_iter, 2)-bkg_width, 1):min(temp(frame_iter, 2)+bkg_width, h), max(temp(frame_iter, 3)-bkg_width, 1):min(temp(frame_iter, 3)+bkg_width, w)) = 1;
            bkg_mask = bkg_mask&~rna_mask;
    
            bkg = mean(double(temp_img(bkg_mask)));
            rna = sum(double(temp_img(rna_mask))-bkg);
            rna_bkg(frame_iter, :) = [rna, bkg];
            rna_bkg_cor(frame_iter, :) = [rna, bkg];
        end
        tracks{track_iter, 3} = rna_bkg;
        tracks{track_iter, 4} = rna_bkg_cor;
        tracks{track_iter, 5} = tracks{track_iter, 2}(:, 4) - rna_bkg_cor(:, 2);
     end
    
    save([filepath(1:(end-4)), '-tracks.mat'], "spots", "tracks");

end
