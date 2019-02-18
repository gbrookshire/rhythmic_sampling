function rs_cfc(i_subject)

% Run the CFC analysis for one subject
% TODO 
%   this can't deal with short segments. Maybe don't downsample first?

rs_setup

% Load the data
fname = subject_info.meg{i_subject};
data = rs_preproc_ress(i_subject, 'trial');
behav = rs_behavior(i_subject); % For RT

%{
% SIMULATED DATA
data = rs_simulate_flicker();
behav = rs_behavior(1); % For RT
fname = 'SIMULATED';
%}

% Set up dir for saving data
save_dir = [exp_dir 'cfc/'];
[~,~,~] = mkdir(save_dir, fname);

% Toss segments that overlap with or occur after the response
% Or that include the transient response at the beginning of the trial
% Load behavioral data
for i_trial = 1:length(data.time)
    % Find the times to keep
    t = data.time{i_trial};
    after_beginning = t > 0.5;
    n_trial = data.trialinfo(i_trial, 2);
    rt = behav.rt(behav.TrialNumber == n_trial);
    if isnan(rt)
        before_resp = ones(size(t));
    else
        before_resp = t < rt;
    end
    before_end = t < exp_params.max_trial_dur;
    keep_samps = after_beginning & before_resp & before_end;
    % Trim the trials
    data.time{i_trial} = data.time{i_trial}(keep_samps);
    data.trial{i_trial} = data.trial{i_trial}(:, keep_samps);
end
% Toss empty trials
cfg = [];
cfg.trials = cellfun(@length, data.time) > 0;
data = ft_selectdata(cfg, data);

% Parameters for the CFC analysis
freq = 55:90;
nfft = 2 ^ 10;
width = 7;

% Compute CFC separately on each freq at its correct side
cfc_data = [];
trial_lens = [];
for side_63 = {'left' 'right'}
    if strcmp(side_63, 'left')
        keep_trials = data.trialinfo(:,3) == 63;
    elseif strcmp(side_63, 'right')
        keep_trials = data.trialinfo(:,3) == 78;
    else
        error('Something went wrong')
    end
    cfg = [];
    cfg.trials = keep_trials;
    data_sub = ft_selectdata(cfg, data);
    [cfc_sub, mod_freq] = cfc2(data_sub, freq, nfft, width);
    cfc_data.(side_63{1}) = cfc_sub; 
    trial_lens.(side_63{1}) = cellfun(@length, data_sub.time);
end

save([save_dir fname '/cfc'], 'cfc_data', 'mod_freq', 'trial_lens')
