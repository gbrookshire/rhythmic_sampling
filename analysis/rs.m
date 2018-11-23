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

rs_apply_over_subjects(@rs_spectra_snr, false)



%% Compute TFRs

tfr_fun = @(i_subj) rs_tfr(i_subj, 'trial');
rs_apply_over_subjects(tfr_fun, true);


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

