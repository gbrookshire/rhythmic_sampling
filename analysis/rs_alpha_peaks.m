% Compute the TFR around alpha peaks

% TODO
% Split by RFT frequency
% Dynamically adjust the colorscale
% Check for FIXMEs

%{
Following Spaak et al (2012)
- Only select segments in which alpha power is > 60th percentile for >= 800 ms
    - This selected ~3% of the data
- Trimmed out 1 s of data centered on these alpha bursts
- BP-filter the alpha signal
    - two-pass FIR-LS filter b/w 7 and 14 Hz (order: 426)
- TFRs of the signal
    - Hanning taper
    - 20-300 Hz in steps of 2 Hz
    - 7 cycles
%}



rs_setup
segment_duration = 1; % s

% Load behavioral data (to get RTs)
behav = rs_behavior(i_subject);

% % Load data
% fname = subject_info.meg{i_subject};
% data_preproc = load([exp_dir 'preproc/trial/' fname '/preproc']);
% data_preproc = data_preproc.data;

data_preproc = rs_preproc_ress(i_subject, 'trial');
fsample = 1 / mean(diff(data_preproc.time{1}));
segment_width = round((segment_duration/2) * fsample);

% Get the band-passed alpha oscillations
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [7 14];
cfg.bpfilttype = 'but'; %FIRLS get a warning: not recommended for neural signals
data_alpha = ft_preprocessing(cfg, data_preproc);

% Get alpha power in order to trim out high-power segments
cfg = [];
cfg.hilbert = 'abs';
data_alpha_pow = ft_preprocessing(cfg, data_alpha);

% Convenience function to concatenate all trials
conc = @(x) cat(2, x{:});

% Get the alpha power cutoffs at each channel 
thrsh_alpha = 60; % Percentile
cutoffs = prctile(conc(data_alpha_pow.trial), thrsh_alpha, 2);

% Look for contiguous segments of data with alpha power above that threshold
thrsh_time = 0.8; % seconds
thrsh_samp = thrsh_time * fsample; % minimum number of samples to keep
% high_power = cellfun(@(x) x > cutoffs, ...
%     data_alpha_pow.trial, ...
%     'UniformOutput', false);
clear segment_info
seg_counter = 0;
for i_trial = 1:length(data_preproc.trial)
    rt = behav(behav.TrialNumber == i_trial,:).rt;
    for i_chan = 1:length(data_preproc.label)
        x = data_alpha_pow.trial{i_trial}(i_chan,:);
        x_thresh = x > cutoffs(i_chan);
        x_change = [0 diff(x_thresh)];
        onsets = (x_thresh == 1) & (x_change == 1);
        offsets = (x_thresh == 0) & (x_change == -1);
        if x_thresh(1) == 1
            onsets(1) = 1;
        end
        if x_thresh(end) == 1
            offsets(end) = 1;
        end
        segments = [find(onsets)' find(offsets)'];
        seg_dur = diff(segments, 1, 2);
        for i_seg = 1:size(segments, 1)
            % Check whether the segment is too early in the trial
            seg_time = data_preproc.time{i_trial};
            seg_time = seg_time(segments(i_seg,:));
            too_early = ~any(seg_time < 0.5);
            % Check whether the segment is long enough
            long_enough = seg_dur(i_seg) > thrsh_samp;
            % Check whether the segment happens before the response
            too_late = any(seg_time > rt);
            if long_enough && ~too_early && ~too_late
                s = [];
                s.n_trial = i_trial;
                s.n_chan = i_chan;
                s.seg_dur = seg_dur(i_seg);
                s.seg_midpoint = round(mean(segments(i_seg,:)));
                s.beg = s.seg_midpoint - segment_width;
                s.end = s.seg_midpoint + segment_width;
                if s.beg > 0 % If this starts after the ft.trial
                    seg_counter = seg_counter + 1;
                    segment_info(seg_counter) = s;
                end
            end
        end
    end
end
clear s seg_time seg_dur segments seg_counter offsets onsets

% Align segments to the nearest alpha peak
max_shift = 0; % Keep track of the biggest shift
for i_seg = 1:length(segment_info)
    s = segment_info(i_seg);
    x = data_alpha.trial{s.n_trial}(s.n_chan, s.beg:s.end);
    % Find the closest alpha peak to the middle of this segment
    [~, locs] = findpeaks(x);
    dist_from_center = abs(locs - segment_width);
    [~, center_inx] = min(dist_from_center);
    % Realign to this peak
    shift = segment_width - locs(center_inx);
    s.seg_midpoint = s.seg_midpoint + shift;
    s.beg = s.beg - shift;
    s.end = s.end - shift;
    segment_info_shifted(i_seg) = s;
    max_shift = max(shift, max_shift); 
end

% Cut out those segments from the alpha data
t = -(segment_duration/2):(1/fsample):(segment_duration/2);
data_seg = [];
data_seg.label = data_preproc.label;
data_seg.trial = {};
data_seg.time = {};
data_seg.trialinfo = nan([0 2]); % Keep track of trial num & targeted chan
for i_seg = 1:length(segment_info_shifted)
    s = segment_info_shifted(i_seg);
    try
        x = data_preproc.trial{s.n_trial}(:, s.beg:s.end);
    catch me
        disp(me.identifier)
        continue
    end
    data_seg.trial{end+1} = x;
    data_seg.time{end+1} = t;
    data_seg.trialinfo(end+1,:) = [s.n_trial, s.n_chan];
end

%{
% Weird -- the alpha phase is pretty different between the two RESS channels
for i_seg = 1:length(data_seg.trial)
    n_chan = segment_info_shifted(i_seg).n_chan;
    subplot(2,1,1) % Plot the selected channel
    plot(t, data_seg.trial{i_seg}(n_chan,:))
    hold on
    subplot(2,1,2) % Plot the non-selected channel
    n_chan = mod(n_chan, 2) + 1;
    plot(t, data_seg.trial{i_seg}(n_chan,:))
    hold on
end
subplot(2,1,1),hold off
subplot(2,1,2),hold off
%}

% Average over segments
% Do this separately for segments at each channel
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [7 14];
cfg.bpfilttype = 'but'; %FIRLS get a warning: not recommended for neural signals
data_seg_alpha = ft_preprocessing(cfg, data_seg);
cfg = [];
cfg.trials = data_seg.trialinfo(:,2) == 1;
data_seg_avg_1 = ft_timelockanalysis(cfg, data_seg_alpha); 
cfg = [];
cfg.trials = data_seg.trialinfo(:,2) == 2;
data_seg_avg_2 = ft_timelockanalysis(cfg, data_seg_alpha);

% Compute the TFR for these segments
% Do this separately for each RESS channel, and each RFT freq
cfg = [];
cfg.method = 'mtmconvol';
cfg.foi = 20:2:100;
cfg.taper = 'hanning';
cfg.t_ftimwin = 7 ./ cfg.foi;
cfg.toi = 'all';

cfg.trials = data_seg.trialinfo(:,2) == 1;
data_seg_tfr_1 = ft_freqanalysis(cfg, data_seg);

cfg.trials = data_seg.trialinfo(:,2) == 2;
data_seg_tfr_2 = ft_freqanalysis(cfg, data_seg);


% Plot it
x_lim = 0.3;
subplot(2,2,1)
cfg = [];
cfg.channel = 1;
cfg.xlim = [-1 1] * x_lim;
cfg.ylim = [40 90];
cfg.colorbar = 'no';
cfg.zlim = [0 4e-28];
ft_singleplotTFR(cfg, data_seg_tfr_1)

subplot(2,2,2)
cfg.channel = 2;
cfg.zlim = [0 1e-28];
ft_singleplotTFR(cfg, data_seg_tfr_2)

subplot(2,2,3) 
plot(t, data_seg_avg_1.avg(1,:))
xlim( [-1 1] * x_lim)

subplot(2,2,4) 
plot(t, data_seg_avg_2.avg(2,:))
xlim( [-1 1] * x_lim)





