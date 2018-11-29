function rs_tagged_xcorr(i_subject)

% Compute cross-correlation between power at the two tagged frequencies

rs_setup

step_size = 0.04; % Should be divisible by 1/Fs to preserve time-bins
toi = 1.0:step_size:5.0;

fname = subject_info.meg{i_subject};
d = rs_preproc(fname, 'trial');

save_dir = [exp_dir 'xcorr/'];
[~,~,~] = mkdir(save_dir, fname);

% TFR around the tagged frequencies
time_window = 0.2; % Smaller window -> more temporal smoothing
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.toi = toi;
cfg.keeptrials = 'yes'; 
cfg.pad = 'nextpow2';
cfg.padtype = 'zero';
cfg.output = 'pow';
cfg.foi = exp_params.tagged_freqs;
cfg.t_ftimwin = ones(length(cfg.foi), 1).* time_window;
d = ft_freqanalysis(cfg, d);

% Compute the cross-correlation
maxlag_sec = 1;
maxlag_samp = round(maxlag_sec / step_size);
xcorr_length = (maxlag_samp * 2) + 1;
x = nan(... % Trial * Channel * Lag
    size(d.powspctrm, 1), ...
    length(d.label), ...
    xcorr_length);
for i_trial = 1:size(d.powspctrm, 1)
    for i_channel = 1:length(d.label)
        % Extract the data
        a = squeeze(d.powspctrm(i_trial,i_channel,1,:))';
        b = squeeze(d.powspctrm(i_trial,i_channel,2,:))';
        % Compute xcorr
        c = nanxcorr(a, b, maxlag_samp);
        x(i_trial, i_channel, :) = c;
    end
end

% Save the data
label = d.label;
time = -maxlag_sec:step_size:maxlag_sec;
save([save_dir '/' fname '/x'], 'x', 'label', 'time')
end

function out = nanxcorr(a, b, maxlag)
% Cross-correlation ignoring NaNs
% This is much slower than the built-in xcorr
lags = -maxlag:maxlag;
out = nan(size(lags));
for i = 1:length(lags)
    t = lags(i);
    out(i) = corr(a', shift(b, t)', 'rows', 'pairwise');
end
end

function out = shift(in, n)
% Shift a vector to the left or right
if n < 0
    out = [in((1 + abs(n)):end), nan(1, abs(n))];
elseif n == 0
    out = in;
else
    out = [nan(1, n) in(1:(end - n))];
end
end
