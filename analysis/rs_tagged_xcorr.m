function rs_tagged_xcorr(i_subject)

% Compute cross-correlation between power at the two tagged frequencies

% rs_setup
% fname = subject_info.meg{i_subject};
% 
% % Load the data
% d = load([exp_dir 'tfr/trial/' fname '/high']);
% d = d.high_freq_data;
% behav = rs_behavior(i_subject); % For RT

% SIMULATED DATA
fname = 'SIMULATED';
d = load('/Users/geoff/Documents/science/2018_perceptual_sampling/sample_data/simulated_high.mat');
d = d.high_freq_data;
behav = rs_behavior(1);

fsample = 1 / mean(diff(d.time));
max_lag_sec = 1; 
max_lag = round(max_lag_sec * fsample); % XCorr lag in samples
t_lags = -max_lag_sec:(1/fsample):max_lag_sec;

% Toss segments that overlap with or occur after the response
% Or that include the transient response at the beginning of the trial
t = d.time;
padding = 0.2; % Give this much leeway around excluded segments
n_trials = size(d.powspctrm, 1);
xc = nan(n_trials, length(t_lags));
inx_63Hz = find(round(d.freq) == 63);
inx_78Hz = find(round(d.freq) == 78);
for i_trial = 1:n_trials
    disp(i_trial)
    n_trial = d.trialinfo(i_trial, 2);
    % Ignore transient response
    after_beginning = t > (0.5 + padding);
    % Exclude times after the stim turns off
    before_end = t < (4.0 - padding);
    % Exclude times after the response
    resp_time = behav.rt(behav.TrialNumber == n_trial);
    before_resp = t < (resp_time - padding);
    % Which samples to keep
    keep_samps = after_beginning & before_end & before_resp;
    if sum(keep_samps) <= 1
        continue
    end
    x = squeeze(d.powspctrm(i_trial,:,:,keep_samps));
    % Get the data from each frequency
    if d.trialinfo(i_trial, 3) == 63
        a = x(1,inx_63Hz,:);
        b = x(2,inx_78Hz,:);
    elseif d.trialinfo(i_trial, 3) == 78
        a = x(1,inx_78Hz,:);
        b = x(2,inx_63Hz,:);
    else
        error('Did not find the frequency of this trial')
    end
    % Compute the cross-correlation
    xc(i_trial,:) = xcorr(squeeze(a), squeeze(b), max_lag);
end


% % TODO
% % Toss the beginning and end of the trials
% % Toss segments after the response
% % Compute cross-correlation separately for each trial (not append)
% 
% keyboard
% 
% xc = nan([2, length(t)]); % 63-Hz side (L/R) x time
% for i_freq = 1:2
%     left_freq = exp_params.tagged_freqs(i_freq);
%     if i_freq == 1
%         right_freq = exp_params.tagged_freqs(2);
%     elseif i_freq == 2
%         right_freq = exp_params.tagged_freqs(1);
%     end
%     % Select trials with this freq mapping
%     cfg = [];
%     cfg.trials = d.trialinfo(:,3) == left_freq;
%     d_sub = ft_selectdata(cfg, d);
%     % Append all the trials in time: Chan x Freq x Time
%     % This is a bad idea -- it introduces artifacts at trial edges
%     [n_trials, n_chans, n_freqs, n_time] = size(d_sub.powspctrm);
%     pwr = nan(n_chans, n_freqs, n_trials * n_time);
%     start_inx = 1; % Where to place this trial in the new matrix
%     for i_trial = 1:n_trials
%         x = squeeze(d_sub.powspctrm(i_trial,:,:,:)); 
%         pwr(:, :, start_inx:(start_inx + n_time - 1)) = x;
%         start_inx = start_inx + n_time;
%     end
%     clear x start_inx
%     % Run the cross-correlation
%     a = squeeze(pwr(1, round(d_sub.freq) == left_freq, :));
%     b = squeeze(pwr(1, round(d_sub.freq) == right_freq, :));
%     xc(i_freq,:) = nanxcorr(a, b, max_lag);
% end

% Save the data
save([exp_dir 'tfr/trial/' fname '/xcorr'], 'xc', 't_lags')
end


function out = nanxcorr(a, b, maxlag)
% Cross-correlation ignoring NaNs
% This is much slower than the built-in xcorr
lags = -maxlag:maxlag;
out = nan([length(lags) 1]);
for i = 1:length(lags)
    t = lags(i);
    out(i) = corr(a, shift(b, t), 'rows', 'pairwise');
end
end


function out = shift(in, n)
% Shift a vector to the left or right
in = reshape(in, [length(in), 1]);
if n < 0
    out = [in((1 + abs(n)):end); nan(abs(n), 1)];
elseif n == 0
    out = in;
else
    out = [nan(n, 1); in(1:(end - n))];
end
end
