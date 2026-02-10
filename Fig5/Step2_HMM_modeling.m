%  Step2_HMM_modeling.m
%  Bo Wang, July 2025
clc;close all;clear;

%  DESCRIPTION
%  This script takes RNA bursting .mat files from different treatment conditions 
%  and performs Hidden Markov Model (HMM) analysis to segment transcriptional states. 
%  Subsequently, one- or two-exponential fitting is applied to the 1-CDF of the 
%  ON and OFF state durations.

%% merged data of different groups
roi_idx = [];
merged_tracks = cell(0);
dirpath = 'D:/ImageData/Untreated/max_projection/';
filelist = dir([dirpath, '*-tracks.mat']);
for file_iter = 1:size(filelist, 1)
    filename = filelist(file_iter).name;
    load([dirpath, filename], "tracks");
    merged_tracks = vertcat(merged_tracks, tracks);
    roi_idx = [roi_idx; size(tracks, 1)];
end
dirpath2 = 'D:/ImageData/JQ1-treat/max_projection/';
filelist2 = dir([dirpath2, '*-tracks.mat']);
for file_iter = 1:size(filelist2, 1)
    filename = filelist2(file_iter).name;
    load([dirpath2, filename], "tracks");
    merged_tracks = vertcat(merged_tracks, tracks);
end

num = [218, 196];
cumsum_num = cumsum(num);
numberOfPages = size(merged_tracks{1, 1}, 1);
roi_idx(:, 2) = cumsum(roi_idx);

%% Plot heatmap
% based on max intensity
rna_intensity_heatmap = zeros(size(merged_tracks, 1), numberOfPages);
for i = 1:size(merged_tracks, 1)
    l = length(merged_tracks{i, 5});
    rna_intensity_heatmap(i, 1:l) = merged_tracks{i, 5}';
end

figure;
h1 = heatmap(rna_intensity_heatmap);
colorused = slanCM(97); colorused = colorused(end:-1:1, :);
colorused = [imresize(colorused(1:128, :), [85, 3], 'bilinear'); colorused(129:256, :)];
h1.Colormap = colorused;
h1.ColorLimits = [0, 6e3];
h1.GridVisible = 'off';
xlabel('Frame');
ylabel('Single Cell');
for i = 1:numberOfPages
    h1.XDisplayLabels{i, 1} = '';
end
for i = 1:size(merged_tracks, 1)
    h1.YDisplayLabels{i, 1} = '';
end

%% hidden Markov model
observations = rna_intensity_heatmap(:, :)';
numStates = 2;
numObservations = length(observations);

% initial parameters estimation
A_est = [0.8, 0.2; 0.2, 0.8];
mu_est = [0; 2000];
sigma_est = [400; 1000];
pi_est = [0.5 0.5];

max_iter = 100;
tol = 1e-6; % Convergence threshold
loglik_prev = -inf;

for iter = 1:max_iter
    % E step: forward-backward algorithm
    % calculate alpha
    alpha = zeros(numStates, numObservations);
    c = zeros(1, numObservations); % scaling factor
    
    % B_obs: The probability of observing a specific observation given a particular state
    % for Gaussian distribution，P(O_t | S_t=s_i) = N(O_t | mu_i, sigma_i^2)
    B_obs = zeros(numStates, numObservations);
    for i = 1:numStates
        B_obs(i, :) = normpdf(observations, mu_est(i), sigma_est(i));
        % To avoid zero probabilities, a very small epsilon can be added.
        B_obs(i, B_obs(i, :) < eps) = eps; 
    end

    % initialize alpha_0
    alpha(:, 1) = pi_est' .* B_obs(:, 1);
    c(1) = sum(alpha(:, 1));
    alpha(:, 1) = alpha(:, 1) / c(1); % scaling
    
    % Backward recursion for alpha calculation
    for t = 2:numObservations
        alpha(:, t) = (A_est' * alpha(:, t-1)) .* B_obs(:, t);
        c(t) = sum(alpha(:, t));
        alpha(:, t) = alpha(:, t) / c(t);
    end

    % calculate beta
    beta = zeros(numStates, numObservations);
    beta(:, end) = 1 / c(end);

    % Backward recursion for beta calculation
    for t = numObservations-1:-1:1
        beta(:, t) = A_est * (B_obs(:, t+1) .* beta(:, t+1));
        beta(:, t) = beta(:, t) / c(t);
    end
    
    % calculate gamma (the posterior probability of each state at each time step)
    % gamma(i, t) = P(S_t=s_i | O, lambda)
    gamma = alpha .* beta;
    gamma = gamma ./ sum(gamma, 1); % normalization

    % calculate xi (the posterior probability of transitioning from state i at time t to state j at time t+1)
    % xi(i, j, t) = P(S_t=s_i, S_{t+1}=s_j | O, lambda)
    xi = zeros(numStates, numStates, numObservations-1);
    for t = 1:numObservations-1
        denom = 0;
        for i = 1:numStates
            for j = 1:numStates
                xi(i, j, t) = alpha(i, t) * A_est(i, j) * B_obs(j, t+1) * beta(j, t+1);
                denom = denom + xi(i, j, t);
            end
        end
        xi(:, :, t) = xi(:, :, t) / denom;
    end

    % M step：update parameters
    pi_est = gamma(:, 1)';

    % A_est update
    for i = 1:numStates
        for j = 1:numStates
            A_est(i, j) = sum(xi(i, j, :), 3) / sum(gamma(i, 1:end-1));
        end
    end
    A_est = A_est ./ sum(A_est, 2); % 行归一化

    % update mean and variance
    for j = 1:numStates
        gamma_sum = sum(gamma(j, :));
        mu_est(j) = sum(gamma(j, :) .* observations) / gamma_sum;
        var_j = sum(gamma(j, :) .* (observations - mu_est(j)).^2) / gamma_sum;
        sigma_est(j) = sqrt(max(var_j, 1e-3)); 
    end

    % Convergence check
    loglik_current = sum(log(c));
    fprintf('Iteration %d: Log-likelihood = %.6f\n', iter, loglik_current);

    if abs(loglik_current - loglik_prev) < tol
        fprintf('Converged at iteration %d\n', iter);
        break;
    end
    loglik_prev = loglik_current;
end

disp('Optimized state transition matrix A:');
disp(A_est);
disp('Optimized initial state probabilities pi:');
disp(pi_est);
disp('Optimized Gaussian means mu:');
disp(mu_est);
disp('Optimized Gaussian standard deviations sigma:');
disp(sigma_est);

% Infer the hidden states of the observation sequence (using Viterbi algorithm)
B_obs_final = zeros(numStates, numObservations);
for i = 1:numStates
    B_obs_final(i, :) = normpdf(observations, mu_est(i), sigma_est(i));
    B_obs_final(i, B_obs_final(i, :) < eps) = eps;
end

delta = zeros(numStates, numObservations);
psi = zeros(numStates, numObservations);

delta(:, 1) = log(pi_est') + log(B_obs_final(:, 1));

for t = 2:numObservations
    for j = 1:numStates
        [max_val, max_idx] = max(delta(:, t-1) + log(A_est(:, j)));
        delta(j, t) = log(B_obs_final(j, t)) + max_val;
        psi(j, t) = max_idx;
    end
end

inferred_hiddenStates = zeros(1, numObservations);
[~, inferred_hiddenStates(numObservations)] = max(delta(:, numObservations));

for t = numObservations-1:-1:1
     inferred_hiddenStates(t) = psi(inferred_hiddenStates(t+1), t+1);
end
%%
inferred_hiddenStates = reshape(inferred_hiddenStates, numberOfPages, []).'-1;
%%
min_intensity_threshold = 1800; % 4*bkgstd
for track_iter = 1:size(merged_tracks, 1)

    merged_tracks{track_iter, 6} = inferred_hiddenStates(track_iter, :)';

    % on and off time, mean,std,max of intensity
    rna_state = merged_tracks{track_iter, 6};
    changes = [0, find(diff(rna_state))', length(rna_state)];
    lengths = diff(changes);
    values = rna_state(changes(2:end));
    merged_tracks{track_iter, 7} = [values, lengths'];

    start_num = 0;
    for i = 1:size(merged_tracks{track_iter, 7}, 1)
        merged_tracks{track_iter, 7}(i, 3) = mean(rna_intensity_heatmap(track_iter, (start_num+1):(start_num+merged_tracks{track_iter, 7}(i, 2)))); % mean
        merged_tracks{track_iter, 7}(i, 4) = max(rna_intensity_heatmap(track_iter, (start_num+1):(start_num+merged_tracks{track_iter, 7}(i, 2)))); % max
        start_num = start_num+merged_tracks{track_iter, 7}(i, 2);
    end

    % remove max_intensity<4*std burst
    merged_tracks{track_iter, 7}(merged_tracks{track_iter, 7}(:, 4)<min_intensity_threshold, 1) = 0;
    new_state = [];
    for i = 1:size(merged_tracks{track_iter, 7}, 1)
        new_state = [new_state; repmat(merged_tracks{track_iter, 7}(i, 1), [merged_tracks{track_iter, 7}(i, 2), 1])];
    end
    merged_tracks{track_iter, 6} = new_state;

    % on and off time, mean,std,max of intensity
    rna_state = merged_tracks{track_iter, 6};
    changes = [0, find(diff(rna_state))', length(rna_state)];
    lengths = diff(changes);
    values = rna_state(changes(2:end));
    merged_tracks{track_iter, 7} = [values, lengths'];

    start_num = 0;
    for i = 1:size(merged_tracks{track_iter, 7}, 1)
        merged_tracks{track_iter, 7}(i, 3) = mean(rna_intensity_heatmap(track_iter, (start_num+1):(start_num+merged_tracks{track_iter, 7}(i, 2)))); % mean
        merged_tracks{track_iter, 7}(i, 4) = max(rna_intensity_heatmap(track_iter, (start_num+1):(start_num+merged_tracks{track_iter, 7}(i, 2)))); % max
        start_num = start_num+merged_tracks{track_iter, 7}(i, 2);
    end

end

%% check burst intensity
% check specified roi
example_idx = [80, 326];%[4, 9, 23, 25]+218;
example_num = length(example_idx);
figure;
for j = 1:example_num
    track_iter = example_idx(j);
    vcount = merged_tracks{track_iter, 7};
    vcount(:, 3) = vcount(:, 1).*vcount(:, 3);
    output = [];
    for i = 1:size(vcount, 1)
        output = [output; repmat(vcount(i, 3), vcount(i, 2), 1)];
    end
    subplot(example_num, 1, j)
    plot(rna_intensity_heatmap(track_iter, :));
    hold on
    plot(output);
    hold off
    % xlim([1, numberOfPages])
    xlim([1, 120])
    ylim([-1e3, 9e3]);
%     xlabel('Frame(30s)');
%     ylabel('RNA intensity (A.U.)');
end

% check random roi
example_num = 36;
rng(123);
example_idx = randsample(length(merged_tracks), example_num);
figure;
for j = 1:example_num
    track_iter = example_idx(j);
    vcount = merged_tracks{track_iter, 7};
    vcount(:, 3) = vcount(:, 1).*vcount(:, 3);
    output = [];
    for i = 1:size(vcount, 1)
        output = [output; repmat(vcount(i, 3), vcount(i, 2), 1)];
    end
    subplot(12, 3, j)
    plot(rna_intensity_heatmap(track_iter, :));
    hold on
    plot(output);
    hold off
    xlim([1, numberOfPages])
    ylim([-1e3, 4e3]);
%     xlabel('Frame(30s)');
%     ylabel('RNA intensity (A.U.)');
end

%% transcription burst 2-component exponential fitting

on_off_time = cell(0);
temp_on_off_time = [];
for track_iter = [1:218]%[219:414]%[1:218] %1:size(merged_tracks, 1)
    temp_on_off_time = [temp_on_off_time; merged_tracks{track_iter, 7}];
end
on_time = temp_on_off_time(temp_on_off_time(:, 1)==1, 2)*30;
off_time = temp_on_off_time(temp_on_off_time(:, 1)==0, 2)*30;
on_off_time{1, 1} = on_time;
on_off_time{1, 2} = off_time;
mean(on_time)
mean(off_time)
temp_on_off_time = [];
for track_iter = [219:414]%[219:414]%[1:218] %1:size(merged_tracks, 1)
    temp_on_off_time = [temp_on_off_time; merged_tracks{track_iter, 7}];
end
on_time = temp_on_off_time(temp_on_off_time(:, 1)==1, 2)*30;
off_time = temp_on_off_time(temp_on_off_time(:, 1)==0, 2)*30;
on_off_time{2, 1} = on_time;
on_off_time{2, 2} = off_time;

mean(on_time)
mean(off_time)
c1 = cdfplot(on_time/60);
on_x = c1.XData; on_y = 1-c1.YData; close;
c2 = cdfplot(off_time/60);
off_x = c2.XData; off_y = 1-c2.YData; close;
figure;
plot(on_x, on_y);
hold on
plot(off_x, off_y);
hold off
legend({'on time', 'off time'});
ylabel('1-CDF');
xlabel('Time(min)');
xlim([0, 30]);

%% ON and OFF time exponential fitting
group_idx = 1; % 1:ctrl; 2:off;
on_off_idx = 1; % 1:on; 2:off;
on_off_labels = {'On Time(s)', 'OFF Time(s)'};
[y,x] = ecdf(on_off_time{group_idx, on_off_idx}); y = 1-y; x(1)=0;
param_fitresult = zeros(9, 3);

figure;
exp1 = fittype('exp(a*x)');
[f,gof,output] = fit(x,y,exp1, ...
    'StartPoint', [-0.01], ...
    'Lower', [-1], ...
    'Upper', [0]);
subplot(2, 3, 1);
plot(f,x,y)
xlabel(on_off_labels{on_off_idx});
ylabel('1-CDF');
title('1-exp');
param_fitresult(1:3, 1) = [1, -f.a, -1/f.a];

subplot(2, 3, 4);
plot(x,output.residuals);
xlabel(on_off_labels{on_off_idx});
ylabel('Residuals');
title(['R-square: ', num2str(gof.rsquare)]);
ylim([-0.05, 0.05]);

exp2 = fittype('a*exp(b*x) + (1-a)*exp(d*x)');
[f,gof,output] = fit(x,y,exp2, ...
    'StartPoint', [0.5, -0.005, -0.01], ...
    'Lower', [0, -1, -1], ...
    'Upper', [1, 0, 0]);
subplot(2, 3, 2);
plot(f,x,y)
xlabel(on_off_labels{on_off_idx});
ylabel('1-CDF');
title('2-exp');
param_fitresult(1:6, 2) = [f.a, -f.b, -1/f.b, 1-f.a, -f.d, -1/f.d];

subplot(2, 3, 5);
plot(x,output.residuals);
xlabel(on_off_labels{on_off_idx});
ylabel('Residuals');
title(['R-square: ', num2str(gof.rsquare)]);
ylim([-0.05, 0.05]);


exp3 = fittype('a*exp(b*x) + c*exp(d*x) + (1-a-c)*exp(f*x)');
[f,gof,output] = fit(x,y,exp3, ...
    'StartPoint', [1, -0.005, 1, -0.01, -0.05], ...
    'Lower', [0, -1, 0, -1, -1], ...
    'Upper', [1, 0, 1, 0, 0]);
subplot(2, 3, 3);
plot(f,x,y)
xlabel(on_off_labels{on_off_idx});
ylabel('1-CDF');
title('3-exp');
param_fitresult(1:9, 3) = [f.a, -f.b, -1/f.b, f.c, -f.d, -1/f.d, 1-f.a-f.c, -f.f, -1/f.f];

subplot(2, 3, 6);
plot(x,output.residuals);
xlabel(on_off_labels{on_off_idx});
ylabel('Residuals');
title(['R-square: ', num2str(gof.rsquare)]);
ylim([-0.05, 0.05]);

mean(on_off_time{group_idx, on_off_idx})

%%
% on-state fitting
on_off_idx = 2; % 1:on; 2:off;
on_off_labels = {'On Time(s)', 'OFF Time(s)'};
[y,x] = ecdf(on_off_time{on_off_idx}); y = 1-y; x(1)=0;

figure;
subplot(1, 2, 1);
exp1 = fittype('exp(a*x)');
[f,gof,output] = fit(x,y,exp1, ...
    'StartPoint', [-0.01], ...
    'Lower', [-1], ...
    'Upper', [0]);
plot(f,x,y)
xlabel(on_off_labels{on_off_idx});
ylabel('1-CDF');
title(['1-exp; R-square: ', num2str(gof.rsquare)]);
subplot(1, 2, 2);
exp2 = fittype('a*exp(b*x) + (1-a)*exp(d*x)');
[f,gof,output] = fit(x,y,exp2, ...
    'StartPoint', [1, -0.005, -0.01], ...
    'Lower', [0, -1, -1], ...
    'Upper', [1, 0, 0]);
plot(f,x,y)
xlabel(on_off_labels{on_off_idx});
ylabel('1-CDF');
title(['2-exp; R-square: ', num2str(gof.rsquare)]);
ylim([0, 1]);