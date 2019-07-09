function rs_baseline_alpha(i_subject)

% Identify channels that show a difference in alpha response after stimuli
% appear on the screen.

rs_setup

fname = subject_info.meg{i_subject};

step_size = 0.05; % Should be divisible by 1/Fs to preserve time-bins
toi = -1:step_size:2;

% Set up dir for saving data
save_dir = [exp_dir 'tfr/trial/'];
[~,~,~] = mkdir(save_dir, fname);

% Load preprocessed data
d = load([exp_dir 'preproc/trial/' fname '/preproc']);
d = d.data;

% Compute TFRs
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.toi = toi;
cfg.keeptrials = 'no';
cfg.output = 'pow';
cfg.foi = 3:30;
n_cycles = 4;
cfg.t_ftimwin = n_cycles ./ cfg.foi;
cfg.pad = 7;
cfg.padtype = 'mirror';
tfr = ft_freqanalysis(cfg, d);

save([save_dir '/' fname '/baseline_alpha'], 'tfr')
