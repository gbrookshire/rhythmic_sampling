function rs_tfr(i_subject, segment_type)

% Compute and save both low- and high-freq TFRs
% segment_type: How to segment the data ('trial' or 'target' or 'repsonse')

% One reason to think about switching to the Hilbert transform for TFRs:
% the FFT basis functions don't land exactly on the rFT tagging freqs.

rs_setup

fname = subject_info.meg{i_subject};

step_size = 0.01; % Should be divisible by 1/Fs to preserve time-bins
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

%{
% TFR around the tagged frequencies
time_window = 0.1; % Smaller window -> more freq smoothing
cfg = cfg_base;
cfg.output = 'pow';
cfg.foi = 55:100;
cfg.t_ftimwin = ones(length(cfg.foi), 1).* time_window;
% % Using virtual channels from RESS spatial filters
% d = rs_preproc_ress(i_subject, segment_type);
% high_freq_data = ft_freqanalysis(cfg, d);
% save([save_dir '/' fname '/high'], 'high_freq_data', '-v7.3')
% clear d cfg high_freq_data
% Using raw channel responses
% Focus on gradiometers that have high SNR
% cfg.channel = {'MEG2042' 'MEG2043' 'MEG2032' 'MEG2033'};
d = rs_preproc(i_subject, segment_type);
high_freq_data = ft_freqanalysis(cfg, d);
save([save_dir '/' fname '/high'], 'high_freq_data', '-v7.3')
clear d cfg high_freq_data
%}

%
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
save([save_dir '/' fname '/low'], 'low_freq_data', '-v7.3')
%}

%{
% TFR at low freqs (theta, alpha) - for standard LF analyses
% Parameters from Popov, Kastner, Jensen, J Neurosci
% For LF data, use all channels
% time_win = 0.5; % From Popov et al
% time_win = 0.3; % For higher frequency smoothing
d = load([exp_dir 'preproc/' segment_type '/' fname '/preproc']);
d = d.data;
cfg = cfg_base;
cfg.output = 'pow';
cfg.foi = 3:30;
% cfg.t_ftimwin = time_win * ones(size(cfg.foi));
n_cycles = 4;
cfg.t_ftimwin = n_cycles ./ cfg.foi;
if strcmp(segment_type, 'target')
    cfg.toi = -1.5:0.05:1;
elseif strcmp(segment_type, 'trial')
    cfg.toi = toi;
else
    error('segment type <%s> not recognized', segment_type)
end
cfg.pad = 7;
cfg.padtype = 'mirror';
low_freq_data = ft_freqanalysis(cfg, d);
save([save_dir '/' fname '/low_standard'], 'low_freq_data', '-v7.3')

%}
