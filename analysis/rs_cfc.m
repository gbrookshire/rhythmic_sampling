function rs_cfc(i_subject)

% Run the CFC analysis for one subject

% Load the data
data = rs_preproc_ress(i_subject, 'trial');
behav = rs_behavior(i_subject); % For RT

% Toss segments that overlap with or occur after the response
% Or that include the transient response at the beginning of the trial
% Load behavioral data
for i_trial = 1:length(data.time)
    % Find the times to keep
    t = data.time{i_trial};
    after_beginning = t > 0.5;
    n_trial = d_side.trialinfo(i_rpt, 2);
    rt = behav.rt(behav.TrialNumber == n_trial);
    before_resp = t < rt;
    before_end = t < exp_params.max_trial_dur;
    keep_samps = after_beginning & before_resp & before_end;
    % Trim the trials
    data.time{i_trial} = data.time{i_trial}(keep_samps);
    data.trial{i_trial} = data.trial{i_trial}(:, keep_samps);
end

% Parameters for the CFC analysis
freq = 55:90;
nfft = 2 ^ 8;
width = 6;

[cfc_data, mod_freq] = cfc(data, freq, nfft, width);
