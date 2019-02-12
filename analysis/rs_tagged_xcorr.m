function rs_tagged_xcorr(i_subject)

% Compute cross-correlation between power at the two tagged frequencies

rs_setup

%{
% SIMULATED DATA
fname = 'SIMULATED';
behav = rs_behavior(1);
%}

% Load the data
fname = subject_info.meg{i_subject};
behav = rs_behavior(i_subject); % For RT
d = load([exp_dir 'tfr/trial/' fname '/high']);
d = d.high_freq_data; 

fsample = 1 / mean(diff(d.time));
max_lag_sec = 0.5; 
max_lag = round(max_lag_sec * fsample); % XCorr lag in samples
t_lags = -max_lag_sec:(1/fsample):max_lag_sec;

% Toss segments that overlap with or occur after the response
% Or that include the transient response at the beginning of the trial
t = d.time;
padding = 0.2; % Give this much leeway around excluded segments
n_trials = size(d.powspctrm, 1);
xc = nan(n_trials, length(t_lags));
trial_length = nan(n_trials, 1);
inx_63Hz = find(round(d.freq) == 63);
inx_78Hz = find(round(d.freq) == 78);
for i_trial = 1:n_trials
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
    trial_length(i_trial) = sum(keep_samps);
    % Avoid crashing if there aren't enough samples to compute cross-corr
    if sum(keep_samps) <= max_lag
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
    a = squeeze(a);
    b = squeeze(b);
    xc(i_trial,:) = nanxcorr(a, b, max_lag);
end 

% Save the data
save([exp_dir 'tfr/trial/' fname '/xcorr'], 'xc', 't_lags', 'trial_length')
end


function out = nanxcorr(a, b, maxlag)
% Cross-correlation ignoring NaNs
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
