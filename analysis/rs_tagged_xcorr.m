function rs_tagged_xcorr(i_subject)

% Compute cross-correlation between power at the two tagged frequencies

rs_setup
fname = subject_info.meg{i_subject};

% Load the data
% behav = rs_behavior(i_subject);
d = load([exp_dir 'tfr/trial/' fname '/high']);
d = d.high_freq_data;

fsample = 1 / mean(diff(d.time));
max_lag_sec = 1; 
max_lag = round(max_lag_sec * fsample); % XCorr lag in samples
t = -max_lag_sec:(1/fsample):max_lag_sec;

xc = nan([2, length(t)]);
for i_freq = 1:2
    left_freq = exp_params.tagged_freqs(i_freq);
    if i_freq == 1
        right_freq = exp_params.tagged_freqs(2);
    elseif i_freq == 2
        right_freq = exp_params.tagged_freqs(1);
    end
    % Select trials with this freq mapping
    cfg = [];
    cfg.trials = d.trialinfo(:,3) == left_freq;
    d_sub = ft_selectdata(cfg, d);
    % Append all the trials in time: Chan x Freq x Time
    [n_trials, n_chans, n_freqs, n_time] = size(d_sub.powspctrm);
    pwr = nan(n_chans, n_freqs, n_trials * n_time);
    start_inx = 1; % Where to place this trial in the new matrix
    for i_trial = 1:n_trials
        x = squeeze(d_sub.powspctrm(i_trial,:,:,:)); 
        pwr(:, :, start_inx:(start_inx + n_time - 1)) = x;
        start_inx = start_inx + n_time;
    end
    clear x start_inx
    % Run the cross-correlation
    a = squeeze(pwr(1, round(d_sub.freq) == left_freq, :));
    b = squeeze(pwr(1, round(d_sub.freq) == right_freq, :));
    xc(i_freq,:) = nanxcorr(a, b, max_lag);
end

% Save the data
save([exp_dir 'tfr/trial/' fname '/xcorr'], 'xc', 't')
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
