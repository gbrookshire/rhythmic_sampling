function rs_tagged_xcorr(i_subject)

% Compute cross-correlation between power at the two tagged frequencies

rs_setup

%{
% SIMULATED DATA
warning('Running analysis with simulated data')
fname = 'SIMULATED';
behav = rs_behavior(1);
%}

tfr_dir = [exp_dir 'tfr/ress/win_0.1s/trial/'];

% Load the data
fname = subject_info.meg{i_subject};
disp(fname)
behav = rs_behavior(i_subject); % For RT

d = load([tfr_dir fname '/high']);
d = d.high_freq_data; 

% Get the time vector for the cross-correlations
fsample = 1 / mean(diff(d.time));
max_lag_sec = 0.5; 
max_lag = round(max_lag_sec * fsample); % XCorr lag in samples
t_lags = -max_lag_sec:(1/fsample):max_lag_sec;

% Compute cross-correlations on all data appended together

% First, pad out trials with extra time so we don't get jumpy
% cross-correlations due to discontinuities from artifact rejection.
pad_dur = 1; % In sec
pad_samp = round(pad_dur * fsample); % How much padding in samples
% Make the padded data
powspctrm_pad = padarray(d.powspctrm, [0 0 0 pad_samp], NaN, 'post');
% Make the time variable
t_pad = (0:(1/fsample):1) + d.time(end);
t_pad = [d.time t_pad(2:end)];
% Combine them into the data structure
d.powspctrm = powspctrm_pad;
d.time = t_pad;
clear pad_dur pad_samp t_pad

% Select trials that were 'hits'
cfg = [];
cfg.trials = d.trialinfo(:,1) == 1;
d = ft_selectdata(cfg, d);

% Get rid of samples that are uninformative (onset transient, after resp)
padding = 0;
for i_trial = 1:size(d.powspctrm, 1)
    trialnum = d.trialinfo(i_trial, 2); % Trial number
    % Ignore transient response
    after_beginning = d.time > (0.5 + padding);
    % Exclude times after the stim turns off
    before_end = d.time < (4.0 - padding);
    % Exclude times after the response
    resp_time = behav.rt(behav.TrialNumber == trialnum);
    if isnan(resp_time) % No resp - keep the whole trial
        before_resp = logical(ones(size(d.time)));
    else
        before_resp = d.time < (resp_time - padding);
    end
    % Which samples to keep
    keep_samps = after_beginning & before_end & before_resp;
    % NaN-out the samples that are uninformative
    d.powspctrm(i_trial, :, :, ~keep_samps) = NaN;
end

% Find which side the stimulus appeared on in each trial
trial_num = d.trialinfo(:,2);
target_side = behav.target_side(trial_num);
freq_left = behav.freq_left(trial_num);
clear trial_num

% Combine the data into one long array & compute cross-correlation
% Do this separately for left-hits and right-hits, and for each RFT mapping
% (i.e. each frequency on each side)
sides = {'left' 'right'};
inx_63Hz = find(round(d.freq) == 63); % Where each tagged freq is
inx_78Hz = find(round(d.freq) == 78);
xc = nan(2, 2, length(t_lags));
for i_target_side = 1:2
    for i_rft_mapping = 1:2
        % Select trials
        target_side_sel = strcmp(target_side, sides{i_target_side});
        rft_map_sel = freq_left == exp_params.tagged_freqs(i_rft_mapping);
        cfg = [];
        cfg.trials = target_side_sel & rft_map_sel;
        d_sel = ft_selectdata(cfg, d);
        x = d_sel.powspctrm;
        % Select freqs and channels
        % Keep in mind that RESS filters are named by *stimulus* side --
        % i.e. the 'left' RESS filter selects channels over the right hemi
        if i_rft_mapping == 1 % 63 Hz on the left side / 78 on right
            left_stim_freq_inx = inx_63Hz;
            right_stim_freq_inx = inx_78Hz;
        elseif i_rft_mapping == 2 % 78 on left / 63 on right
            left_stim_freq_inx = inx_78Hz;
            right_stim_freq_inx = inx_63Hz;
        end
        x_right_hemi = x(:, 1, left_stim_freq_inx, :);
        x_left_hemi = x(:, 2, right_stim_freq_inx, :);
        x_right_hemi = squeeze(x_right_hemi);
        x_left_hemi = squeeze(x_left_hemi);
        % Append into a long vector
        longvec = @(a) reshape(a', [1 numel(a)])';
        x_right_hemi = longvec(x_right_hemi);
        x_left_hemi = longvec(x_left_hemi);
        % Compute cross-correlation
        xx = nanxcorr(x_right_hemi, x_left_hemi, max_lag);
        xc(i_target_side, i_rft_mapping, :) = xx;
     end
end
% Average over RFT mappings
xc = squeeze(mean(xc, 2));

% Save the data
save_dir = [exp_dir 'xcorr/'];
[~,~,~] = mkdir(save_dir, fname);
save([save_dir '/' fname '/xc'], 'xc', 't_lags')


end


function out = nanxcorr(a, b, maxlag)
% Cross-correlation that is robust to NaNs (ignores them)
% This is much slower than the built-in xcorr
lags = -maxlag:maxlag;
out = nan([length(lags) 1]);
for i = 1:length(lags)
    t = lags(i);
    b_shift = shift(b, t);
    if ~isequal(size(a), size(b_shift))
        keyboard
    end
    out(i) = corr(a, b_shift, 'rows', 'pairwise');
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
