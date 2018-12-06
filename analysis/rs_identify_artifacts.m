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
i_subject = 16;

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

[~,~,~] = mkdir([exp_dir 'artifacts/'], fname);
save([exp_dir 'artifacts/' fname '/visual'], 'bad_chans', 'cfg_art')


%% Use ICA to identify eye and cardiac artifacts

rs_setup

fname = subject_info.meg{i_subject};
a = load([exp_dir 'artifacts/' fname '/visual']);

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
save([exp_dir 'artifacts/' fname '/ica'], 'comp', 'reject_comp')


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

save([exp_dir 'artifacts/' fname '/photodiode'], ...
    'anomalies', 'photo_artfctdef')


%% Look for eye-blinks in the Eyetracker recordings

clear variables
rs_setup

amp_threshold = 0.2;
size_threshold = 5;
exclude_pre = 0.1;
exclude_post = 0.1;

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};

    % Read in the raw eyelink recordings
    eyes_artfctdef = cell(size(block_info.all));
    for i_block = block_info.main

        % Read in the trial definition
        fn = [exp_dir 'trialdef/' fname '/' num2str(i_block) '.mat'];
        trialdef = load(fn);

        % Preprocess the data
        fn = [exp_dir 'raw/' fname '/' num2str(i_block) '.fif'];
        cfg = [];
        cfg.dataset = fn;
        cfg.event = trialdef.trl.event;
        cfg.trl = trialdef.trl.trial;
        cfg.channel = {'MISC001' 'MISC002'};
        cfg.polyremoval = 'yes';
        cfg.polyorder = 0;
        d = ft_preprocessing(cfg);

        % Find the eye blink artifacts in this subject
        eye_artifact = nan(0, 2);
        for i_trial = 1:length(d.trial)
            flag_samps = all(abs(d.trial{i_trial}) > amp_threshold, 1);
            flag_regions = bwconncomp(flag_samps);
            flag_region_size = cellfun(@length, flag_regions.PixelIdxList);
            blink_regions = find(flag_region_size >= size_threshold);
            blink_samps = flag_regions.PixelIdxList(blink_regions);
            blink_starts = cellfun(@(x) x(1), blink_samps);
            blink_ends = cellfun(@(x) x(end), blink_samps);
            art_def = [blink_starts' - (exclude_pre * d.fsample), ...
                       blink_ends' + (exclude_post * d.fsample)];
            % Adjust for the sample that the trial starts on
            art_def = art_def + d.sampleinfo(i_trial, 1);
            eye_artifact = [eye_artifact; art_def];
        end
        
        % Merge overlapping artifacts
        eye_art_old = eye_artifact;
        eye_art_merged = eye_art_old(1,:);
        clear eye_artifact
        for i_art = 2:size(eye_art_old, 1)
            if eye_art_old(i_art, 1) < eye_art_merged(end, 2)
                eye_art_merged(end, 2) = eye_art_old(i_art, 2);
            else
                eye_art_merged(end + 1, :) = eye_art_old(i_art, :);
            end
        end
        
        % Add to the overall cell object
        eyes_artfctdef{i_block} = eye_art_merged;

        % % Check the length of the artifacts
        % hist(diff(eye_artifact, 1, 2), 50);
    end

    save([exp_dir 'artifacts/' fname '/eye'], 'eyes_artfctdef')
    
end
    

