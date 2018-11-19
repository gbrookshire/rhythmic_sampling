% Rhythmic sampling
% Fall 2018
% Geoff Brookshire
% Each subject has multiple files
%   [1-2].fif is the QUEST thresholding procedure
%   [3-5].fif are the main experiment

clear variables
close all

if strcmp(computer, 'PCWIN64')
    cd('Z:/gb/')
else
    cd('/rds/projects/2017/jenseno-02/gb/')
    run('~/startup.m')
end
addpath('rhythmic_sampling/analysis/')
addpath('matlab_helpers/')
rs_setup
cd(exp_dir)

addpath([base_dir 'fieldtrip-20181118/'])
ft_defaults


%% TODO

% Reject segments with blinks
% Look at microsaccades
%     Check eye-tracking traces
%     Check EOG traces


%% Define trials -- 2nd version
rs_setup
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    [~,~,~] = mkdir([exp_dir 'trialdef/'], fname);
    trl_dir = [exp_dir 'trialdef/' fname '/'];

    for i_block = block_info.all
        % Common setup for segmenting based on different events
        dataset = [exp_dir 'raw/' fname '/' num2str(i_block) '.fif'];
        if ~exist(dataset, 'file')
            warning('No MEG data file: %s\\%d.fif', fname, i_block)
            input('Press ENTER to continue')
            continue
        end
        
        trl = rs_trialfun2(dataset);
        save([trl_dir num2str(i_block)], 'trl')
    end
end


%% Save grad structures for each subject
rs_setup
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    dataset = [exp_dir 'raw/' fname '/1.fif'];
    hdr = ft_read_header(dataset);
    grad = hdr.grad;
    [~,~,~] = mkdir([exp_dir 'grad/'], fname);
    save([exp_dir 'grad/' fname '/grad'], 'grad')
end


%% Make forward model

%%%% Check the cfg options for all of this!

rs_setup
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    % Read in the MRI scan
    mri = ft_read_mri(subject_info.T1(i_subject));
    % Segmentation
    cfg = [];
    cfg.write = 'no';
    segmentedmri = ft_volumesegment(cfg, mri);
    % Prepare the head model
    cfg = [];
    cfg.method = 'singleshell';
    headmodel = ft_prepare_headmodel(cfg, segmentedmri);
    % Prepare the lead field
    cfg = [];
    cfg.grad = FREQ_DATA.grad; % Fill in with output from ft_freqanalysis
    cfg.headmodel = headmodel;
    cfg.reducerank = 2;
    cfg.channel = FREQ_DATA.label;
    cfg.grid.resolution = 1; % use a 3-D grid with a 1 cm resolution
    cfg.grid.unit = 'cm';
    grid = ft_prepare_leadfield(cfg);
    % Save the results
    save([exp_dir 'forward_model\' fname '\fm'], 'headmodel', 'grid')
end


%% Identify artifacts
%  When you combine different recordings with ft_appenddata, the sample
%  info gets confused because the numbers for each recording all are
%  indexed to the beginning of the recording, so the sample indices
%  overlap. This doesn't affect the data when plotting it directly, but it
%  does cause ft_databrowser to get confused and jumble the data. The
%  takeaway message is: you can't use ft_databrowser after ft_appenddata. I
%  suspect this also prevents you from re-sampling the data later on.
%
% Work-around:
% - for each subject
%   - for each recording file
%       - visually reject artifacts
%       - save individual artifact defs
%   - combine recordings
%       - ICA

clear variables
i_subject = 14;

%% Visually identify artifacts

rs_setup

fname = subject_info.meg{i_subject};
disp(fname)

bad_chans = [];
cfg_art = [];
% for chan_type = {'grad' 'mag'} % Separately examine grads and mags
for chan_type = {'grad'} % Separately examine grads and mags
    chan_type = chan_type{1};
    disp(chan_type)
    for i_block = block_info.main
        disp(['Block ' num2str(i_block)])

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
        cfg.channel = chan.(chan_type).names;
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [48 52];
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 100;
        cfg.padding = 8;
        cfg.padtype = 'data';
        cfg.polyremoval = 'yes';
        cfg.polyorder = 2;
        d = ft_preprocessing(cfg);

        % Plot the summary stats and look for bad channels
        cfg = [];
        cfg.method = 'summary';
        cfg.layout = chan.(chan_type).layout;
        data_sel_rej = ft_rejectvisual(cfg, d);
        rej_chans_summ = setdiff(d.label, data_sel_rej.label);
        
        % Detailed look at the signals
        % Tag segments that are artifacts and look for other bad channels
        cfg = [];
        cfg.channel = data_sel_rej.label;
        cfg.viewmode = 'vertical';
        cfg.continuous = 'no';
        cfg.layout = chan.(chan_type).layout;
        cfg_art.(chan_type){i_block} = ft_databrowser(cfg, d);
        
        % Put together the list of bad channels
        rej_chans_visual = input('Other bad channels: ');
        rej_chans_visual = reshape(rej_chans_visual, ... % Reshape to
            [length(rej_chans_visual), 1]);              % combine
        bad_ones = [rej_chans_summ rej_chans_visual];
        bad_ones = reshape(bad_ones, [1 length(bad_ones)]);
        bad_chans.(chan_type){i_block} = bad_ones;
    end
end

[~,~,~] = mkdir([exp_dir 'artifacts\'], fname);
save([exp_dir 'artifacts\' fname '\visual'], 'bad_chans', 'cfg_art')


%% Use ICA to identify eye and cardiac artifacts

fname = subject_info.meg{i_subject};
a = load([exp_dir 'artifacts\' fname '\visual']);

% Only keep MEG channels and reject bad channels
if isfield(a.bad_chans, 'mag')
    bad_chans = union( ...
        union([a.bad_chans.grad{:}], {}), ...
        union([a.bad_chans.mag{:}], {}));
else
    bad_chans = union([a.bad_chans.grad{:}], {});
end
bad_chans = cellfun(@(s) ['-' s], bad_chans, 'UniformOutput', false);
chan_sel = [chan.all.names bad_chans'];

data_by_block = cell(size(block_info.all));

% Reject visually-identified artifacts identified above
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
    cfg.channel = chan_sel;
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [48 52];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 100;
    cfg.padding = 8;
    cfg.padtype = 'data';
    cfg.polyremoval = 'yes';
    cfg.polyorder = 2;
    d = ft_preprocessing(cfg);

    % Reject the visually-tagged artifacts
    cfg = [];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.visual_grad = a.cfg_art.grad{i_block}.artfctdef.visual;
    %cfg.artfctdef.visual_mag = a.cfg_art.mag{i_block}.artfctdef.visual;
    data_by_block{i_block} = ft_rejectartifact(cfg, d);
end

% Combine data
% Delete the sample info so it doesn't mix up data from diff trials
for i_block = block_info.main
    data_by_block{i_block} = rmfield(data_by_block{i_block}, 'sampleinfo');
end
% Only keep the data from the main blocks
data_by_block = data_by_block(block_info.main);
data_all_blocks = ft_appenddata([], data_by_block{:});
clear data_by_block

% Downsample
cfg = [];
cfg.resamplefs = 250;
data_all_blocks = ft_resampledata(cfg, data_all_blocks);

% Run ICA
cfg = [];  
cfg.method = 'runica';  
cfg.runica.maxsteps = 100;  
comp = ft_componentanalysis(cfg, data_all_blocks);   

% Visualize the components
cfg = [];
cfg.component = 1:20;
cfg.layout = chan.mag.layout;
cfg.comment = 'no';
ft_topoplotIC(cfg, comp)

% Plot the components with their timecourses
cfg = [];
cfg.channel = 1:15; 
cfg.continuous = 'no';
cfg.viewmode = 'component'; 
cfg.layout = chan.mag.layout;
cfg.compscale = 'local';
ft_databrowser(cfg, comp);

% Save the list of components to reject
reject_comp = input('Components to reject: ');
save([exp_dir 'artifacts\' fname '\ica'], 'comp', 'reject_comp')


%% Look for irregularities in the photo-diode
% Identify them by looking for changes in frequency
rs_setup
ignore_before = 0.5; % Ignore anomalies before this time in the trial
ignore_after = 1.5; % Ignore anom after 0.5 s from the end of the trial
exclude_pre = 0.1; % Exclude this much time before the photodiode jump
exclude_post = 0.1; % Exclude this much time after

fname = subject_info.meg{i_subject};

% Flag any jumps in the photodiode above this threshold
threshold = 0.3;

% Read in the raw photodiode recordings
phot = cell(size(block_info.all));
for i_block = block_info.all

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
    cfg.channel = 'MISC004'; % Only the photodiode
    cfg.detrend = 'yes'; % Remove slow drifts 
    d = ft_preprocessing(cfg);
    phot{i_block} = d;
    d
end

anomalies = cell(size(phot));
photo_artfctdef = cell(size(phot));

for i_block = block_info.all
    % Get the instantaneous phase
    phi = cellfun(@(x) angle(hilbert(x)), ...
        phot{i_block}.trial, ...
        'UniformOutput', false);
    block_anomalies = cell(1, length(phi));
    i_trial = 1;
    while i_trial <= length(phi)
        % Get the instantaneous frequency
        dphi_dt = diff(phi{i_trial});
        % Wrap negative values around for easy comparison
        dphi_dt(dphi_dt < 0) = dphi_dt(dphi_dt < 0) + (2 * pi);
        % Look for changes in freq
        dphi2_d2t = diff(dphi_dt);
        % Threshold to find anomalous values - i.e. skipped frames
        anom = abs(dphi2_d2t) > threshold; 
        % Add back in the samples lost to derivative operations (diff)
        anom_pad = [anom false false];
        % Plot the photodiode anomalies
        subplot(2,1,1)
        plot(phot{i_block}.trial{i_trial}) % Raw photodiode rec
        hold on
        t = 1:length(phot{i_block}.trial{i_trial});
        plot(t(anom_pad), ...
            phot{i_block}.trial{i_trial}(anom_pad), ...
            'or')
        hold off
        ylabel('Raw photodiode (V)')
        xlim([1 6000])
        ylim([-0.025 0.025])
        subplot(2,1,2)
        plot(abs(dphi2_d2t)) % Changes in frequency
        hold on
        t = 1:length(dphi2_d2t);
        plot(t(anom), abs(dphi2_d2t(anom)), 'or')
        hold off
        ylabel('abs(d^{2}\phi / dt^{2}) (rad)')
        xlim([1 6000])
        ylim([0 6])
        fprintf('%i: ', i_trial)
        r = input('[N]ext or [p]revious: ', 's');
        switch r
            case {'', 'n'} % Save the result and go on
                block_anomalies{i_trial} = anom_pad;
                i_trial = i_trial + 1;
            case 'p' % Go to the previous trial
                i_trial = i_trial - 1;
        end
    end
    anomalies{i_block} = block_anomalies;
    
    % Make an artifact def object as would be output by ft_databrowser
    warning('Make sure photodiode artifacts are correctly tagged')
    photo_artifact = [];
    fsample = phot{i_block}.fsample;
    for i_trial = 1:length(phot{i_block}.trial)
        anom = anomalies{i_block}{i_trial}; % Photodiode errors
        t = phot{i_block}.time{i_trial}; % Time in the trial
        % Ignore changes at the very beginning of the trial
        anom(t < ignore_before) = 0;
        % Ignore changes at the end of the trial
        anom(t > (phot{i_block}.time{i_trial}(end)) - ignore_after) = 0;
        trial_start_samp = phot{i_block}.sampleinfo(i_trial, 1);
        art_samps = find(anom) + trial_start_samp;
        % Toss tagged samples that are too close to each other
        overlap_thresh = 10;
        art_samps([inf diff(art_samps)] < overlap_thresh) = [];
        % Build a mat to put in cfg.artfctdef.xxx.artifact
        % Nx2: [start_samp end_samp]
        art_def = [art_samps' - (exclude_pre * fsample), ...
            art_samps' + (exclude_post * fsample)];
        photo_artifact = [photo_artifact; art_def];
    end
    photo_artfctdef{i_block} = photo_artifact;
    
end

save([exp_dir 'artifacts\' fname '\photodiode'], ...
    'anomalies', 'photo_artfctdef')


%% Compute spectra of MEG signals for each trial
% Note that this doesn't reject visually identified artifacts!

clear variables
rs_setup

for i_subject = 1:height(subject_info)

    fname = subject_info.meg{i_subject};
    art_ica = load([exp_dir 'artifacts\' fname '\ica']); % Artifact defs

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

    [~,~,~] = mkdir([exp_dir 'spectra\'], fname);
    save([exp_dir 'spectra\' fname '\spectra'], 'freq_data')
end


%% Compute the SNR at each channel to get an ROI of power at tagged freqs

rs_setup
for i_subject = 8:height(subject_info)

    fname = subject_info.meg{i_subject};
    spec = load([exp_dir 'spectra\' fname '\spectra']);
    spec = spec.freq_data;

    f_tag = exp_params.tagged_freqs; % Tagged frequencies

    snr = cell(size(f_tag));
    for i_freq = 1:length(f_tag)
        snr{i_freq} = rs_snr(spec, f_tag(i_freq));
    end

    save([exp_dir 'spectra\' fname '\snr'], 'snr')
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

    save_dir = [exp_dir 'tfr/%s/' segment_type '/'];
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
    parsave([sprintf(save_dir, 'high') fname '/x'], high_freq_data)
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
    parsave([sprintf(save_dir, 'low') fname '/x'], high_freq_data)
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

