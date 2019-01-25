function rs_derive_ress(i_subject)

% Derive the RESS/GEDb spatial filters

rs_setup

fname = subject_info.meg{i_subject};
data = load([exp_dir 'preproc/trial/' fname '/preproc']);
data_preproc = data.data; clear data;
behav = rs_behavior(i_subject);

% Compute and plot the maps and spectra
ress_maps = [];
for i_freq = 1:2
    for side = {'left' 'right'}
        filter_freq = exp_params.tagged_freqs(i_freq);

        % Select trials with consistent freq/side mapping
        fieldname = ['freq_' side{1}];
        keep_trials = ismember(...
            data_preproc.trialinfo(:,2), ...
            find(behav.(fieldname) == filter_freq));
        cfg = [];
        cfg.trials = keep_trials;
        data_sub = ft_selectdata(cfg, data_preproc);

        % Compute RESS components
        [data_ress, maps, ress] = rs_ress(data_sub, filter_freq, 0.5);
        o = [];
        o.maps = maps;
        o.ress = ress;
        o.label = data_preproc.label;
        ress_maps.(side{1}).(['f' num2str(filter_freq)]) = o;
    end
end

save_dir = [exp_dir 'ress/'];
[~,~,~] = mkdir(save_dir, fname);
save([save_dir '/' fname '/ress'], 'ress_maps')
