% Rhythmic sampling
% Fall 2018
% Geoff Brookshire
% Each subject has multiple files
%   [1-2].fif is the QUEST thresholding procedure
%   [3-5].fif are the main experiment

clear variables
close all
rs_setup


%% TODO

% Reject segments with blinks
% Look at microsaccades
%     Check eye-tracking traces
%     Check EOG traces


%% Define the trials
rs_apply_over_subjects(@rs_definetrials, false)


%% Save grad structures for each subject
rs_apply_over_subjects(@rs_grad_struct, false)


%% Make forward model
rs_apply_over_subjects(@rs_forward_model, true)


%% Identify artifacts
% Use the code in rs_identify_artifacts

%% Compute spectra of MEG signals for each trial
% Note that this doesn't reject visually identified artifacts because that
% eliminates whole trials
rs_apply_over_subjects(@rs_spectra, true)


%% Compute the SNR at each channel to get an ROI of power at tagged freqs

rs_setup
for i_subject = 8:height(subject_info)

    fname = subject_info.meg{i_subject};
    spec = load([exp_dir 'spectra/' fname '/spectra']);
    spec = spec.freq_data;

    f_tag = exp_params.tagged_freqs; % Tagged frequencies

    snr = cell(size(f_tag));
    for i_freq = 1:length(f_tag)
        snr{i_freq} = rs_snr(spec, f_tag(i_freq));
    end

    save([exp_dir 'spectra/' fname '/snr'], 'snr')
end


%% Compute TFRs

% One reason to think about switching to the Hilbert transform for TFRs:
% the FFT basis functions don't land exactly on the rFT tagging freqs.

clear variables
rs_setup

segment_type = 'trial'; % 'trial' or 'target' or 'repsonse'

if strcmp(segment_type, 'trial')
    toi = -0.5:0.05:1.5;
else
    toi = -0.5:0.05:0.5;
end

for i_subject = 1:height(subject_info)
    
    if subject_info.exclude(i_subject)
        continue
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
    parsave([save_dir '/' fname '/high'], high_freq_data)
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
    parsave([save_dir '/' fname '/low'], low_freq_data)
    clear cfg low_freq_data
    
end


%% Compute power at the tagged freqs using BP-filters and Hilbert transform

clear variables
rs_setup

segment_type = 'trial'; % Or target

save_dir = [exp_dir 'tfr/hilbert/high_freq/' segment_type '/'];
parfor i_subject = 1:height(subject_info)
    
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    d = rs_preproc(fname, segment_type);

    h = cell(1,2);
    for i_freq = exp_params.tagged_freqs
        f = exp_params.tagged_freqs(i_freq);
        % Get power with bandpass filter and Hilbert transform
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = f + [-2 2];
        %cfg.bpfilttype = 'fir'; % Better than 'but' for Hilbert phase est.
        cfg.hilbert = 'abs';
        cfg.padding = 3;
        cfg.padtype = 'mirror';
        h{i_freq} = ft_preprocessing(d, cfg);
    end
    
    [~,~,~] = mkdir(save_dir, fname);
    parsave([save_dir fname '/high_freq/'], h)
end









%% Does power at the tagged freqs depend on theta/alpha phase?



%% Relationship of power at the tagged frequencies

% Cross-corr of power at the two tagged frequencies.
% Check whether the strength of that cross-corr is modulated by theta/alpha
% phase. Any better way to do this than with a median split?

