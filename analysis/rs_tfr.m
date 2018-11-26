function rs_tfr(i_subject, segment_type)

% Compute and save both low- and high-freq TFRs
% segment_type: How to segment the data ('trial' or 'target' or 'repsonse')

% One reason to think about switching to the Hilbert transform for TFRs:
% the FFT basis functions don't land exactly on the rFT tagging freqs.

rs_setup

if strcmp(segment_type, 'trial')
    toi = -0.5:0.05:1.5; %%%% Would this work for -0.5:0.05:6?
else
    toi = -0.5:0.05:0.5;
end

fname = subject_info.meg{i_subject};
d = rs_preproc(fname, segment_type);

save_dir = [exp_dir 'tfr/' segment_type '/'];
[~,~,~] = mkdir(save_dir, fname);

% Set up the basic cfg options for both freq bands
cfg_base = [];
cfg_base.method = 'mtmconvol';
cfg_base.taper = 'hanning';
cfg_base.toi = toi;
cfg_base.keeptrials = 'yes'; 

% TFR around the tagged frequencies
time_window = 0.1; % Smaller window -> more temporal smoothing
cfg = cfg_base;
cfg.output = 'pow';
cfg.foi = 55:100;
cfg.t_ftimwin = ones(length(cfg.foi), 1).* time_window;
high_freq_data = ft_freqanalysis(cfg, d);
save([save_dir '/' fname '/high'], 'high_freq_data')
clear cfg high_freq_data

% TFR at low freqs (theta, alpha)
n_cycles = 3;
cfg = cfg_base;
cfg.output = 'fourier'; % Get phase with `angle(...)`
cfg.foi = 4:13;
cfg.t_ftimwin = n_cycles ./ cfg.foi;
cfg.pad = 7; % Pad trials out to 7 sec
cfg.padtype = 'mirror'; % Is this OK for estimating phase?
low_freq_data = ft_freqanalysis(cfg, d);
save([save_dir '/' fname '/low'], 'low_freq_data')
   