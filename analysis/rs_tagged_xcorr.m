function rs_tagged_xcorr(i_subject)

% Compute cross-correlation between power at the two tagged frequencies

rs_setup
fname = subject_info.meg{i_subject};

% Load the data
behav = rs_behavior(i_subject);
d = load([exp_dir 'tfr/trial/' fname '/high']);
d = d.high_freq_data;

% % Append all the trials in time: Chan x Freq x Time
% [n_trials, n_chans, n_freqs, n_time] = size(d.powspctrm);
% pwr = nan(n_chans, n_freqs, n_trials * n_time);
% start_inx = 1; % Where to place this trial
% for i_trial = 1:n_trials
%     x = squeeze(d.powspctrm(i_trial,:,:,:)); 
%     pwr(:, :, start_inx:(start_inx + n_time - 1)) = x;
%     start_inx = start_inx + n_time;
% end
% clear x
% 
% % Compute the cross-correlations
% %%% THIS DOESN'T ACTUALLY MAKE SENSE
% % We want xcorr between channel NOT WITHIN FREQ, but between RFT freqs
% fsample = mean(diff(d.time));
% max_lag = round(1 * (1 / fsample)); % Maximum XCorr lag in samples
% xc = nan(n_freqs, 2 * max_lag + 1);
% for i_freq = 1:n_freqs
%     a = squeeze(pwr(1, i_freq, :));
%     b = squeeze(pwr(2, i_freq, :));
%     x = nanxcorr(a, b, max_lag);
%     xc(i_freq, :) = x;
% end
% keyboard

% Do this for trials with 63 Hz on the left
a = squeeze(pwr(1, round(d.freq) == 63, :)); % L-RESS, 63 Hz
b = squeeze(pwr(2, round(d.freq) == 78, :)); % R-RESS, 78 Hz
x = nanxcorr(a, b, max_lag);

% And this for trials with 63 Hz on the right
a = squeeze(pwr(2, round(d.freq) == 63, :)); % L-RESS, 63 Hz
b = squeeze(pwr(1, round(d.freq) == 78, :)); % R-RESS, 78 Hz
x = nanxcorr(a, b, max_lag);




% 
% 
% 
% 
% 
% 
% 
% 
% 
% rs_setup
% 
% step_size = 0.04; % Should be divisible by 1/Fs to preserve time-bins
% toi = 0.5:step_size:(exp_params.max_trial_dur-0.5);
% 
% fname = subject_info.meg{i_subject};
% d = rs_preproc(fname, 'trial');
% 
% save_dir = [exp_dir 'xcorr/'];
% [~,~,~] = mkdir(save_dir, fname);
% 
% % TFR around the tagged frequencies
% % time_window = 0.1; % 10 Hz smoothing
% time_window = 0.2; % 5 Hz smoothing
% % time_window = 0.5; % 2 Hz smoothing
% cfg = [];
% cfg.method = 'mtmconvol';
% cfg.taper = 'hanning';
% cfg.toi = toi;
% cfg.keeptrials = 'yes'; 
% cfg.pad = 'nextpow2';
% cfg.padtype = 'zero';
% cfg.output = 'pow';
% cfg.foi = exp_params.tagged_freqs;
% cfg.t_ftimwin = ones(length(cfg.foi), 1) .* time_window;
% d = ft_freqanalysis(cfg, d);
% 
% % Compute the cross-correlation
% maxlag_sec = 1;
% maxlag_samp = round(maxlag_sec / step_size);
% xcorr_length = (maxlag_samp * 2) + 1;
% x = nan(... % xcorr: Trial * Channel * Lag
%     size(d.powspctrm, 1), ...
%     length(d.label), ...
%     xcorr_length);
% nsamps = nan([1 size(d.powspctrm, 1)]); % Number of samples in each trial
% for i_trial = 1:size(d.powspctrm, 1)
%     for i_channel = 1:length(d.label)
%         % Extract the data
%         a = squeeze(d.powspctrm(i_trial,i_channel,1,:))';
%         b = squeeze(d.powspctrm(i_trial,i_channel,2,:))';
%         % Compute xcorr
%         c = nanxcorr(a, b, maxlag_samp);
%         x(i_trial, i_channel, :) = c;
%     end
%     nsamps(i_trial) = sum(~isnan(d.powspctrm(i_trial,1,1,:)));
% end
% 
% % Save the data
% label = d.label;
% time = -maxlag_sec:step_size:maxlag_sec;
% save([save_dir '/' fname '/x'], 'x', 'label', 'time', 'nsamps')
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
