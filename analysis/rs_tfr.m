function rs_tfr(i_subject, segment_type)

% Compute and save both low- and high-freq TFRs
% segment_type: How to segment the data ('trial' or 'target' or 'repsonse')

% One reason to think about switching to the Hilbert transform for TFRs:
% the FFT basis functions don't land exactly on the rFT tagging freqs.

rs_setup

fname = subject_info.meg{i_subject};

step_size = 0.02; % Should be divisible by 1/Fs to preserve time-bins
if strcmp(segment_type, 'trial')
    toi = -0.5:step_size:(exp_params.max_trial_dur + 0.5);
else
    toi = -0.5:step_size:0.5;
end

%{
% SIMULATED DATA
d = rs_simulate_flicker();
fname = 'SIMULATED';
%}

% Set up dir for saving data
save_dir = [exp_dir 'tfr/' segment_type '/'];
[~,~,~] = mkdir(save_dir, fname);

% Set up the basic cfg options for both freq bands
cfg_base = [];
cfg_base.method = 'mtmconvol';
cfg_base.taper = 'hanning';
cfg_base.toi = toi;
cfg_base.keeptrials = 'yes'; 
% cfg_base.pad = 'nextpow2';
% cfg_base.padtype = 'zero';

% TFR around the tagged frequencies
For HF data, use virtual channels from RESS spatial filters
d = rs_preproc_ress(i_subject, segment_type);
time_window = 0.2; % Smaller window -> more freq smoothing
cfg = cfg_base;
cfg.output = 'pow';
cfg.foi = 55:100;
cfg.t_ftimwin = ones(length(cfg.foi), 1).* time_window;
high_freq_data = ft_freqanalysis(cfg, d);
save([save_dir '/' fname '/high'], 'high_freq_data', '-v7.3')
clear d cfg high_freq_data

% TFR at low freqs (theta, alpha)
% For LF data, use all channels
d = load([exp_dir 'preproc/' segment_type '/' fname '/preproc']);
d = d.data;
n_cycles = 2; % As in Fiebelkorn et al (2018, Neuron)
cfg = cfg_base;
cfg.output = 'fourier'; % Get phase with `angle(...)`
cfg.foi = 3:13;
cfg.t_ftimwin = n_cycles ./ cfg.foi;
cfg.pad = 7; % Pad trials out to 7 sec
cfg.padtype = 'mirror'; % Is this OK for estimating phase?
low_freq_data = ft_freqanalysis(cfg, d);
save([save_dir '/' fname '/low'], 'low_freq_data')
