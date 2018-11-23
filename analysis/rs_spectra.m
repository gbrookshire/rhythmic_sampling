function rs_spectra(i_subject)

rs_setup

fname = subject_info.meg{i_subject};
art_ica = load([exp_dir 'artifacts/' fname '/ica']); % Artifact defs

data_by_block = cell(size(block_info.all));
for i_block = block_info.main

    % Read in the trial definition
    fn = [exp_dir 'trialdef/' fname '/' num2str(i_block) '.mat'];
    if ~exist(fn, 'file')
        warning('No trialdef for sub %s, block %d', fname, i_block)
        continue
    end
    trialdef = load(fn);

    % Preprocess the data
    fn = [exp_dir 'raw/' fname '/' num2str(i_block) '.fif'];
    cfg = [];
    cfg.dataset = fn;
    cfg.event = trialdef.trl.event;
    cfg.trl = trialdef.trl.trial;   
    cfg.channel = art_ica.comp.cfg.channel; % Load the good channels
    cfg.polyremoval = 'yes';
    cfg.polyorder = 1;
    d = ft_preprocessing(cfg);

    % Downsample
    cfg = [];
    cfg.resamplefs = 250;
    d = ft_resampledata(cfg, d);

    % Reject artifact ICs
    cfg = [];
    cfg.component = art_ica.reject_comp;
    d = ft_rejectcomponent(cfg, art_ica.comp, d);

    % Save the preproc data in a cell obj
    data_by_block{i_block} = d;
    clear d
end

% Combine all the data
data_by_block = data_by_block(block_info.main);
data = ft_appenddata([], data_by_block{:});

% Only look at the time while the stimulus is flashing
for i_trial = 1:length(data.trial)
    t = data.time{i_trial};
    t_sel = (t > 0) & (t < (max(t) - 1)); % Assuming 1s post-stim in trial
    d = data.trial{i_trial};
    data.trial{i_trial} = d(:,t_sel);
    data.time{i_trial} = t(:,t_sel);
end

% Compute the spectra
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.channel = 'MEG';
cfg.taper = 'hanning';
cfg.pad = 'nextpow2'; % (2 ^ 11) / (250 Hz) = 8.192 s
cfg.padtype = 'zero';
cfg.polyremoval = 1; % Remove linear trends
freq_data = ft_freqanalysis(cfg, data);
warning('TODO: Make sure this correctly pads trials to 4.096 s')

[~,~,~] = mkdir([exp_dir 'spectra/'], fname);
save([exp_dir 'spectra/' fname '/spectra'], 'freq_data')
