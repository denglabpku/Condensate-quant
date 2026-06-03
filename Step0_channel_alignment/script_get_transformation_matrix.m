%  Project: Multi-channel Bead-based Chromatic Registration
%  Contact: Bo Wang, Peking University
%
%  Description:
%  This script performs multi-channel spatial registration using 3D fluorescent
%  bead calibration datasets acquired in both super-resolution (SR) and
%  confocal imaging modes.
%
%  The workflow includes:
%    1. Loading bead localization coordinates from four channels
%       (405 nm, 488 nm, 561 nm, 642 nm).
%    2. Pairwise bead matching across channels based on spatial proximity.
%    3. Identification of beads detected simultaneously in all four channels.
%    4. Removal of outlier bead pairs using distance-based filtering.
%    5. Estimation of 3D affine and 2D similarity geometric transformations.
%    6. Transformation of non-reference channels into the 488 nm reference frame.
%    7. Quantitative evaluation of registration accuracy by residual
%       displacement analysis in x, y, and z dimensions.
%
%  Output:
%  The script saves channel registration transformation matrices for
%  subsequent multi-color image alignment and chromatic aberration correction.

clc;close all;clear; rng(42);

%% ============================================================
%% Super-resolution (SR) mode
%% ============================================================

dist_thresh = 0.3;  % Maximum allowed pairing distance (μm)

% Containers for pooled bead coordinates across all datasets
total_pos_405 = []; total_pos_488 = []; total_pos_561 = []; total_pos_642 = [];

% List of SR calibration datasets
filename_list = {'D:\ImageData\BeadsCalibration\20260121\20260121_Beads_SR_4channel_004', ...
                 'D:\ImageData\BeadsCalibration\20260121\20260121_Beads_SR_4channel_005', ...
                 'D:\ImageData\BeadsCalibration\20260121\20260121_Beads_SR_4channel_006', ...
                 'D:\ImageData\BeadsCalibration\20260123\20260123-Beads-4channel-SR-004', ...
                 'D:\ImageData\BeadsCalibration\20260123\20260123-Beads-4channel-SR-005', ...
                 'D:\ImageData\BeadsCalibration\20260123\20260123-Beads-4channel-SR-006'};

for filename_iter = 1:length(filename_list)

    filename = filename_list{filename_iter};

    % --------------------------------------------------------
    % Load bead localization coordinates detected by TrackMate
    % Columns 5–7 correspond to x, y, z coordinates
    % The first row is skipped because it contains metadata
    % --------------------------------------------------------

    beads_pos_405 = readtable([filename, '_405channel.csv']);
    beads_pos_405 = table2array(beads_pos_405(2:end, 5:7));
    
    beads_pos_488 = readtable([filename, '_488channel.csv']);
    beads_pos_488 = table2array(beads_pos_488(2:end, 5:7));
    
    beads_pos_561 = readtable([filename, '_560channel.csv']);
    beads_pos_561 = table2array(beads_pos_561(2:end, 5:7));
    
    beads_pos_642 = readtable([filename, '_640channel.csv']);
    beads_pos_642 = table2array(beads_pos_642(2:end, 5:7));

    % --------------------------------------------------------
    % Pair beads between 488 nm and the other channels
    % 488 nm is used as the reference channel
    % --------------------------------------------------------
    
    [paired_pos_488, paired_pos_561] = position_pairing(beads_pos_488,beads_pos_561,dist_thresh,0);
    [paired_pos_488_2, paired_pos_642] = position_pairing(beads_pos_488,beads_pos_642,dist_thresh,0);
    [paired_pos_488_3, paired_pos_405] = position_pairing(beads_pos_488, beads_pos_405, dist_thresh, 0);

    % --------------------------------------------------------
    % Identify beads detected simultaneously in all four channels
    %
    % Because each pairing step may retain a different subset of beads,
    % coordinate strings are used as unique identifiers
    % --------------------------------------------------------
    
    key1 = string(paired_pos_488(:,1)) + "_" + string(paired_pos_488(:,2));
    if size(paired_pos_488,2) == 3
        key1 = key1 + "_" + string(paired_pos_488(:,3));
    end
    key2 = string(paired_pos_488_2(:,1)) + "_" + string(paired_pos_488_2(:,2));
    if size(paired_pos_488_2,2) == 3
        key2 = key2 + "_" + string(paired_pos_488_2(:,3));
    end
    key3 = string(paired_pos_488_3(:,1)) + "_" + string(paired_pos_488_3(:,2));
    if size(paired_pos_488_3,2) == 3
        key3 = key3 + "_" + string(paired_pos_488_3(:,3));
    end
    [key_tmp, idx1, idx2] = intersect(key1, key2);
    [key_all, idx3, idx4] = intersect(key_tmp, key3);
    
    idx_561 = idx1(idx3);
    idx_642 = idx2(idx3);
    idx_405 = idx4;

    % --------------------------------------------------------
    % Append matched bead coordinates to pooled arrays
    % --------------------------------------------------------
    
    total_pos_405 = [total_pos_405; paired_pos_405(idx_405, :)];
    total_pos_488 = [total_pos_488; paired_pos_488(idx_561, :)];
    total_pos_561 = [total_pos_561; paired_pos_561(idx_561, :)];
    total_pos_642 = [total_pos_642; paired_pos_642(idx_642, :)];

end

id = abs(total_pos_488(:, 3)-total_pos_405(:, 3))<0.3;
total_pos_405 = total_pos_405(id, :);
total_pos_488 = total_pos_488(id, :);
total_pos_561 = total_pos_561(id, :);
total_pos_642 = total_pos_642(id, :);

% ------------------------------------------------------------
% Add a shared random z-offset to all channels
%
% This improves numerical stability for 3D affine fitting when
% beads are concentrated within a narrow z-plane
% ------------------------------------------------------------

dz = -0.5 + 2 * rand(size(total_pos_405, 1), 1);
total_pos_405(:, 3) = total_pos_405(:, 3) + dz;
total_pos_488(:, 3) = total_pos_488(:, 3) + dz;
total_pos_561(:, 3) = total_pos_561(:, 3) + dz;
total_pos_642(:, 3) = total_pos_642(:, 3) + dz;

% ------------------------------------------------------------
% Estimate 3D affine transformations
%
% Naming convention:
% tform_SR_3D_A_B transforms channel B into channel A
% ------------------------------------------------------------

tform_SR_3D_488_405 = fitgeotform3d(total_pos_405, total_pos_488, 'affine');
tform_SR_3D_561_405 = fitgeotform3d(total_pos_405, total_pos_561, 'affine');
tform_SR_3D_488_561 = fitgeotform3d(total_pos_561, total_pos_488, 'affine');
tform_SR_3D_561_642 = fitgeotform3d(total_pos_642, total_pos_561, 'affine');
tform_SR_3D_488_642 = fitgeotform3d(total_pos_642, total_pos_488, 'affine');

% ------------------------------------------------------------
% Estimate 2D similarity transformations using x-y coordinates
% ------------------------------------------------------------

tform_SR_2D_488_405 = fitgeotform2d(total_pos_405(:, 1:2), total_pos_488(:, 1:2), 'similarity');
tform_SR_2D_561_405 = fitgeotform2d(total_pos_405(:, 1:2), total_pos_561(:, 1:2), 'similarity');
tform_SR_2D_488_561 = fitgeotform2d(total_pos_561(:, 1:2), total_pos_488(:, 1:2), 'similarity');
tform_SR_2D_561_642 = fitgeotform2d(total_pos_642(:, 1:2), total_pos_561(:, 1:2), 'similarity');
tform_SR_2D_488_642 = fitgeotform2d(total_pos_642(:, 1:2), total_pos_488(:, 1:2), 'similarity');

%% ============================================================
% Confocal mode
% ============================================================

% The confocal pipeline follows the same procedure as the SR mode.

dist_thresh = 0.3; %um

total_pos_405 = []; total_pos_488 = []; total_pos_561 = []; total_pos_642 = [];

filename_list = {'D:\ImageData\BeadsCalibration\20260121\20260121_Beads_noSR_4channel_007', ...
                 'D:\ImageData\BeadsCalibration\20260121\20260121_Beads_noSR_4channel_008', ...
                 'D:\ImageData\BeadsCalibration\20260121\20260121_Beads_noSR_4channel_009', ...
                 'D:\ImageData\BeadsCalibration\20260123\20260123-Beads-4channel-noSR-001', ...
                 'D:\ImageData\BeadsCalibration\20260123\20260123-Beads-4channel-noSR-002', ...
                 'D:\ImageData\BeadsCalibration\20260123\20260123-Beads-4channel-noSR-007', ...
                 'D:\ImageData\BeadsCalibration\20260123\20260123-Beads-4channel-noSR-008', ...
                 };

for filename_iter = 1:length(filename_list)

    filename = filename_list{filename_iter};

    beads_pos_405 = readtable([filename, '_405channel.csv']);
    beads_pos_405 = table2array(beads_pos_405(2:end, 5:7));
    
    beads_pos_488 = readtable([filename, '_488channel.csv']);
    beads_pos_488 = table2array(beads_pos_488(2:end, 5:7));
    
    beads_pos_561 = readtable([filename, '_560channel.csv']);
    beads_pos_561 = table2array(beads_pos_561(2:end, 5:7));
    
    beads_pos_642 = readtable([filename, '_640channel.csv']);
    beads_pos_642 = table2array(beads_pos_642(2:end, 5:7));

    [paired_pos_488, paired_pos_561] = position_pairing(beads_pos_488,beads_pos_561,dist_thresh,0);
    [paired_pos_488_2, paired_pos_642] = position_pairing(beads_pos_488,beads_pos_642,dist_thresh,0);
    [paired_pos_488_3, paired_pos_405] = position_pairing(beads_pos_488, beads_pos_405, dist_thresh, 0);
    
    key1 = string(paired_pos_488(:,1)) + "_" + string(paired_pos_488(:,2));
    if size(paired_pos_488,2) == 3
        key1 = key1 + "_" + string(paired_pos_488(:,3));
    end
    key2 = string(paired_pos_488_2(:,1)) + "_" + string(paired_pos_488_2(:,2));
    if size(paired_pos_488_2,2) == 3
        key2 = key2 + "_" + string(paired_pos_488_2(:,3));
    end
    key3 = string(paired_pos_488_3(:,1)) + "_" + string(paired_pos_488_3(:,2));
    if size(paired_pos_488_3,2) == 3
        key3 = key3 + "_" + string(paired_pos_488_3(:,3));
    end
    [key_tmp, idx1, idx2] = intersect(key1, key2);
    [key_all, idx3, idx4] = intersect(key_tmp, key3);
    
    idx_561 = idx1(idx3);
    idx_642 = idx2(idx3);
    idx_405 = idx4;
    
    total_pos_405 = [total_pos_405; paired_pos_405(idx_405, :)];
    total_pos_488 = [total_pos_488; paired_pos_488(idx_561, :)];
    total_pos_561 = [total_pos_561; paired_pos_561(idx_561, :)];
    total_pos_642 = [total_pos_642; paired_pos_642(idx_642, :)];

end

total_pos_405 = total_pos_405(id, :);
total_pos_488 = total_pos_488(id, :);
total_pos_561 = total_pos_561(id, :);
total_pos_642 = total_pos_642(id, :);

dz = -0.5 + 2 * rand(size(total_pos_405, 1), 1);
total_pos_405(:, 3) = total_pos_405(:, 3) + dz;
total_pos_488(:, 3) = total_pos_488(:, 3) + dz;
total_pos_561(:, 3) = total_pos_561(:, 3) + dz;
total_pos_642(:, 3) = total_pos_642(:, 3) + dz;

tform_confocal_3D_488_405 = fitgeotform3d(total_pos_405, total_pos_488, 'affine');
tform_confocal_3D_561_405 = fitgeotform3d(total_pos_405, total_pos_561, 'affine');
tform_confocal_3D_488_561 = fitgeotform3d(total_pos_561, total_pos_488, 'affine');
tform_confocal_3D_561_642 = fitgeotform3d(total_pos_642, total_pos_561, 'affine');
tform_confocal_3D_488_642 = fitgeotform3d(total_pos_642, total_pos_488, 'affine');

tform_confocal_2D_488_405 = fitgeotform2d(total_pos_405(:, 1:2), total_pos_488(:, 1:2), 'similarity');
tform_confocal_2D_561_405 = fitgeotform2d(total_pos_405(:, 1:2), total_pos_561(:, 1:2), 'similarity');
tform_confocal_2D_488_561 = fitgeotform2d(total_pos_561(:, 1:2), total_pos_488(:, 1:2), 'similarity');
tform_confocal_2D_561_642 = fitgeotform2d(total_pos_642(:, 1:2), total_pos_561(:, 1:2), 'similarity');
tform_confocal_2D_488_642 = fitgeotform2d(total_pos_642(:, 1:2), total_pos_488(:, 1:2), 'similarity');

% ============================================================
% Save transformation matrices
% ============================================================

save('./ImageRegistor/tform_20260123.mat', 'tform_SR_3D_488_405', 'tform_SR_3D_561_405', 'tform_SR_3D_488_561', 'tform_SR_3D_561_642', 'tform_SR_3D_488_642', ...
                                           'tform_SR_2D_488_405', 'tform_SR_2D_561_405', 'tform_SR_2D_488_561', 'tform_SR_2D_561_642', 'tform_SR_2D_488_642', ...
                                           'tform_confocal_3D_488_405', 'tform_confocal_3D_561_405', 'tform_confocal_3D_488_561', 'tform_confocal_3D_561_642', 'tform_confocal_3D_488_642', ...
                                           'tform_confocal_2D_488_405', 'tform_confocal_2D_561_405', 'tform_confocal_2D_488_561', 'tform_confocal_2D_561_642', 'tform_confocal_2D_488_642');
%% transform 3D beads position

total_pos_405_warped = transformPointsForward(tform_SR_3D_488_405, total_pos_405);
total_pos_561_warped = transformPointsForward(tform_SR_3D_488_561, total_pos_561);
total_pos_642_warped = transformPointsForward(tform_SR_3D_488_642, total_pos_642);

%% calculate 3-channel distance in x-y-z axis

d_488_405 = total_pos_488-total_pos_405; d_488_405_warped = total_pos_488-total_pos_405_warped;
d_488_561 = total_pos_488-total_pos_561; d_488_561_warped = total_pos_488-total_pos_561_warped;
d_488_642 = total_pos_488-total_pos_642; d_488_642_warped = total_pos_488-total_pos_642_warped;
d_561_642 = total_pos_561-total_pos_642; d_561_642_warped = total_pos_561_warped-total_pos_642_warped;

% ------------------------------------------------------------
% Plot residual displacement histograms along x, y, z
%
% Blue: before correction
% Orange: after correction
% ------------------------------------------------------------

title_label = {'x-axis', 'y-axis', 'z-axis'};
bin_width = [5e-3, 5e-3, 1e-2];
figure;
for i = 1:3
subplot(4, 3, i);
histogram(d_488_405(:, i), 'BinWidth',bin_width(i));
hold on
histogram(d_488_405_warped(:, i), 'BinWidth',bin_width(i));
hold off
title(title_label{i})
xlabel('Distance(um)');
end
for i = 1:3
subplot(4, 3, i+3);
histogram(d_488_561(:, i), 'BinWidth',bin_width(i));
hold on
histogram(d_488_561_warped(:, i), 'BinWidth',bin_width(i));
hold off
title(title_label{i})
xlabel('Distance(um)');
end
for i = 1:3
subplot(4, 3, i+6);
histogram(d_561_642(:, i), 'BinWidth',bin_width(i));
hold on
histogram(d_561_642_warped(:, i), 'BinWidth',bin_width(i));
hold off
title(title_label{i})
xlabel('Distance(um)');
end
for i = 1:3
subplot(4, 3, i+9);
histogram(d_488_642(:, i), 'BinWidth',bin_width(i));
hold on
histogram(d_488_642_warped(:, i), 'BinWidth',bin_width(i));
hold off
title(title_label{i})
xlabel('Distance(um)');
end
legend({'Before Correction', 'After Correction'});

%%

mean(d_488_405_warped)
std(d_488_405_warped)
mean(d_488_561_warped)
std(d_488_561_warped)
mean(d_488_642_warped)
std(d_488_642_warped)
mean(d_561_642_warped)
std(d_561_642_warped)

%%
figure;
subplot(1, 2, 1);
hold on
scatter(total_pos_405(:, 1), total_pos_405(:, 2), 10, 'b', 'filled');
scatter(total_pos_488(:, 1), total_pos_488(:, 2), 10, 'r', 'filled');
scatter(total_pos_561(:, 1), total_pos_561(:, 2), 10, 'g', 'filled');
scatter(total_pos_642(:, 1), total_pos_642(:, 2), 10, 'm', 'filled');
hold off
xlabel('X / μm'); ylabel('Y / μm');
subplot(1, 2, 2);
hold on
scatter(total_pos_405_warped(:, 1), total_pos_405_warped(:, 2), 10, 'b', 'filled');
scatter(total_pos_488(:, 1), total_pos_488(:, 2), 10, 'r', 'filled');
scatter(total_pos_561_warped(:, 1), total_pos_561_warped(:, 2), 10, 'g', 'filled');
scatter(total_pos_642_warped(:, 1), total_pos_642_warped(:, 2), 10, 'm', 'filled');
hold off
xlabel('X / μm'); ylabel('Y / μm');
legend({'405nm', '488nm', '561nm', '642nm'});

figure;
subplot(1, 2, 1);
scatter3(total_pos_405(:, 1), total_pos_405(:, 2), total_pos_405(:, 3), 10, 'b', 'filled');
hold on
scatter3(total_pos_488(:, 1), total_pos_488(:, 2), total_pos_488(:, 3), 10, 'r', 'filled');
scatter3(total_pos_561(:, 1), total_pos_561(:, 2), total_pos_561(:, 3), 10, 'g', 'filled');
scatter3(total_pos_642(:, 1), total_pos_642(:, 2), total_pos_642(:, 3), 10, 'm', 'filled');
hold off
xlabel('X / μm'); ylabel('Y / μm'); zlabel('Z / μm');
subplot(1, 2, 2);

scatter3(total_pos_405_warped(:, 1), total_pos_405_warped(:, 2), total_pos_405_warped(:, 3), 10, 'b', 'filled');
hold on
scatter3(total_pos_488(:, 1), total_pos_488(:, 2), total_pos_488(:, 3), 10, 'r', 'filled');
scatter3(total_pos_561_warped(:, 1), total_pos_561_warped(:, 2), total_pos_561_warped(:, 3), 10, 'g', 'filled');
scatter3(total_pos_642_warped(:, 1), total_pos_642_warped(:, 2), total_pos_642_warped(:, 3), 10, 'm', 'filled');
hold off
xlabel('X / μm'); ylabel('Y / μm'); zlabel('Z / μm');
% legend({'405nm', '488nm', '561nm', '642nm'});

