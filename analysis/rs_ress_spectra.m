function rs_ress_spectra(i_subject)

% Compute the spectra of activity at the RESS/GEDb filters

rs_setup

fname = subject_info.meg{i_subject};
data = load([exp_dir 'preproc/trial/' fname '/preproc']);
data_preproc = data.data; clear data;
grad = load([exp_dir 'grad/' fname '/grad'], 'grad');
behav = rs_behavior(i_subject);
ress_maps = load([exp_dir 'ress/' fname '/ress']);
ress_maps = ress_maps.ress_maps;

spectra = [];
for i_freq = 1:2
    for side = {'left' 'right'}
        filter_freq = exp_params.tagged_freqs(i_freq);
        ress_filter = ress_maps.(side{1}).(['f' num2str(filter_freq)]).ress;

        % Select trials with consistent freq/side mapping
        fieldname = ['freq_' side{1}];
        keep_trials = ismember(...
            data_preproc.trialinfo(:,2), ...
            find(behav.(fieldname) == filter_freq));
        cfg = [];
        cfg.trials = keep_trials;
        data_sub = ft_selectdata(cfg, data_preproc);
        data_ress = rs_applyressfilt(data_sub, ress_filter);

        % Compute the spectra
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'pow';
        cfg.taper = 'hanning';
        cfg.pad = 'nextpow2';
        cfg.padtype = 'zero';
        cfg.polyremoval = 1; % Remove linear trends
        spec = ft_freqanalysis(cfg, data_ress);
        spectra.(side{1}).(['f' num2str(filter_freq)]) = spec;
    end
end

save([exp_dir 'ress/' fname '/spectra'], 'spectra')