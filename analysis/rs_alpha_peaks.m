function rs_alpha_peaks(i_subject)

% Compute the TFR around alpha peaks

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
segment_duration = 0.8; % s
thrsh_alpha = 40; % Keep segments of data with alpha pow above this percentile
thrsh_time = 0.8; % Keep segments above the alpha threshold for this dur (s)

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
cfg = []; % Look at slower fluctuations in power (widen the window)
cfg.lpfilter = 'yes';
cfg.lpfreq = 3;
data_alpha_pow = ft_preprocessing(cfg, data_alpha_pow);

% Convenience function to concatenate all trials
conc = @(x) cat(2, x{:});

% Get the alpha power cutoffs at each channel 
cutoffs = prctile(conc(data_alpha_pow.trial), thrsh_alpha, 2);

% Look for contiguous segments of data with alpha power above that threshold
thrsh_samp = thrsh_time * fsample; % minimum number of samples to keep
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

% Detrend each segment
cfg = [];
cfg.detrend = 'yes';
data_seg = ft_preprocessing(cfg, data_seg);

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
% Keep track of resuls in a cell: Channel * TaggedFreq
avg_sels = cell([length(data_seg.label) length(exp_params.tagged_freqs)]);
avg_counts = nan(size(avg_sels));

for i_chan = 1:length(data_seg.label)
    for i_tagfreq = 1:2 % Which tagged freq is on the left side
        % Select trials where the segment occured in this channel
        chan_trial_sel = data_seg.trialinfo(:,2) == i_chan;
        % Select trials where the segment occured at this tagged freq
        frq = exp_params.tagged_freqs(i_tagfreq);
        trial_num_by_seg = data_seg.trialinfo(:,1);
        rft_trial_sel = behav{trial_num_by_seg, 'freq_left'} == frq;
        % Put together the trial selections
        trial_sel = chan_trial_sel & rft_trial_sel;
        % Select the channel
        chan_sel = zeros(size(data_seg.label));
        chan_sel(i_chan) = 1;
        chan_sel = logical(chan_sel);
        % Save the relevant information
        sels = [];
        sels.trial = trial_sel;
        sels.channel = chan_sel;
        avg_sels{i_chan, i_tagfreq} = sels;
        avg_counts(i_chan, i_tagfreq) = sum(trial_sel);
    end
end

fname = subject_info.meg{i_subject};
fn = [exp_dir 'alpha_peaks/' strrep(fname, '/', '_')];
save(fn, 'data_seg', 'avg_sels', 'avg_counts')
