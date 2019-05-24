function [f, c] = rs_pcm(i_subject)

% Compute the pairwise phase consistency metric (PCM) on the time-course of
% the difference in power between target and non-target stimuli following
% Landau et al (2015)

% Load the data
p = rs_powerdiff(i_subject, 0.1, 'target', false);

% Fourier transform each trial
n_timepoints = floor(length(p.time) / 2); % Only pre-target samples
win = hanning(n_timepoints);

% x = p.powdiff(:,1:n_timepoints); %%%% Old version
x = [p.powdiff_hit(:,1:n_timepoints); p.powdiff_miss(:,1:n_timepoints)];

nfft = 2 ^ ceil(log2(n_timepoints));
Fs = 1 / mean(diff(p.time));
f = Fs * (0:(nfft / 2)) / nfft;
y = fft(x .* win', nfft, 2);
y = y(:,1:nfft/2+1);

% Sort trials into hit and miss
% hit_inx = find(p.trialinfo(:,1) == 1); %%% Old version
% miss_inx = find(p.trialinfo(:,1) == 0);
hit_inx = 1:size(p.powdiff_hit, 1);
miss_inx = (1:size(p.powdiff_miss, 1)) + length(hit_inx);
combos = combvec(hit_inx', miss_inx');

% For each possible pairing of hit and miss trials, calculate the cosine of 
% the  phase  difference, and then average over all possible pairs of hit
% and miss trials.
y_hit = y(combos(1,:),:);
y_miss = y(combos(2,:),:);
phase_diff = angle(y_hit) - angle(y_miss);
cos_dist = cos(phase_diff);
c = nanmean(cos_dist, 1);