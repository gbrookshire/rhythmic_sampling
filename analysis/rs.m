% Rhythmic sampling
% Fall 2018
% Geoff Brookshire
% Each subject has multiple files
%   [1-2].fif is the QUEST thresholding procedure
%   [3-5].fif are the main experiment

%% Define the trials
rs_apply_over_subjects(@rs_definetrials, true)


%% Save grad structures for each subject
rs_apply_over_subjects(@rs_grad_struct, false)


%% Make forward model
rs_apply_over_subjects(@rs_forward_model, true)


%% Identify artifacts
% Use the script `rs_identify_artifacts`


%% Compute spectra of MEG signals for each trial
% Note that this doesn't reject visually identified artifacts because that
% eliminates whole trials
rs_apply_over_subjects(@rs_spectra, true)


%% Compute the SNR at each channel to get an ROI of power at tagged freqs
rs_apply_over_subjects(@rs_spectra_snr, false)


%% Save preprocessed data
for segment_event = {'trial' 'target'}
    preproc_fun = @(i_subj) rs_preproc(i_subj, segment_event{1});
    rs_apply_over_subjects(preproc_fun, true)
end

% Submit this to slurm using slurm_party.py

%% Compute RESS/GEDb spatial filters
rs_apply_over_subjects(@rs_derive_ress, true)

%%%%%%% CURRENTLY HERE IN THE ANALYSES %%%%%%%%%


%% Compute TFRs
% High- and low-frequency
for segment_event = {'trial' 'target'}
    tfr_fun = @(i_subj) rs_tfr(i_subj, segment_event{1});
    rs_apply_over_subjects(tfr_fun, false);
end


%% Does power at the tagged frequencies vary rhythmically?
rs_apply_over_subjects(@rs_tagged_spect, true)


%% Compute CFC - coherence between raw signal and power @ tagged freqs
rs_apply_over_subjects(@rs_cfc, true)

%% Compute cross-correlations between power at tagged freqs
rs_apply_over_subjects(@rs_tagged_xcorr, false)


%% Did power at the tagged frequencies differ b/w hits and misses?
rs_apply_over_subjects(@rs_acc_hfpower_regression, true);


%% Does accuracy vary with the phase of LF oscillations?
rs_apply_over_subjects(@rs_acc_lfphase_fiebelkorn, false)
rs_apply_over_subjects(@rs_acc_lfphase_pbi, false)


%% Does power at the tagged freq depend on the phase of LF oscillations?
rs_apply_over_subjects(@FILL_IN, true)








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



