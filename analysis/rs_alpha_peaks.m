function rs_alpha_peaks(i_subject)

% Compute the TFR around alpha peaks

%{
v2
- Look at all segments of data
- Timelock to each alpha peak
- Then start to throw away segments with lower power
- Compute power on all data
- Then re-segment
%}

rs_setup
segment_duration = 0.8; % s

% Load behavioral data (to get RTs)
behav = rs_behavior(i_subject);

% Load preprocessed data
data_preproc = rs_preproc_ress(i_subject, 'trial');

fsample = 1 / mean(diff(data_preproc.time{1}));
segment_width = round((segment_duration/2) * fsample);

% Compute TFRs over all data
cfg.method = 'mtmconvol';
cfg.foi = 55:1:90;
cfg.taper = 'hanning';
cfg.t_ftimwin = 6 ./ cfg.foi;
cfg.toi = 'all';
cfg.keeptrials = 'yes';
data_tfr = ft_freqanalysis(cfg, data_preproc);

% Get the band-passed alpha oscillations
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [7 14];
cfg.bpfilttype = 'but'; %FIRLS get a warning: not recommended for neural signals
data_alpha = ft_preprocessing(cfg, data_preproc);

% Find overlapping segments centered around alpha peaks
% Go through each trial, and look for alpha peaks
% Identify center, beginning, and ending of each segment
% We need to be careful to keep different kinds of 'trials' distinc
% - exp_trial: trial numbers as counted by the experiment script
% - preproc_trial: FT-trial numbers as counted by rs_preproc
% - segments: Snippets of data (coded as 'trials' in a FT object)
segments = nan([0 7]); % Preproc Trial, chan, center, begin, end, t_begin, t_end
for i_preproc_trial = 1:length(data_alpha.trial)
    t = data_alpha.time{i_preproc_trial};
    exp_trial_num = data_alpha.trialinfo(i_preproc_trial,2); % trial in orig exp
    for i_chan = 1:length(data_alpha.label)
        [~, centers] = findpeaks(data_alpha.trial{i_preproc_trial}(i_chan,:), ...
            'MinPeakHeight', 0);
        beginnings = centers - segment_width;
        endings = centers + segment_width;
        % Chuck segments that would start before there's data, or that
        % would end after the data ends
        too_early = beginnings <= 0;
        too_late = endings > length(t);
        chuck_segments = too_early | too_late;
        centers(chuck_segments) = [];
        beginnings(chuck_segments) = [];
        endings(chuck_segments) = [];
        % Check whether the segments overlap with resps, stim onset, etc
        t_beg = t(beginnings);
        t_end = t(endings);
        after_stim_onset = t_beg > 0.5; % Avoid transients from stim onset
        rt = behav(behav.TrialNumber == exp_trial_num, :).rt;
        before_rt = t_end < rt;
        before_trial_end = t_end < 4.0;
        chuck_segments = ~(after_stim_onset & before_rt & before_trial_end);
        centers(chuck_segments) = [];
        beginnings(chuck_segments) = [];
        endings(chuck_segments) = [];
        t_beg(chuck_segments) = [];
        t_end(chuck_segments) = [];
        % Round to the nearest ms and go slightly before
        % To make sure this catches the right samples
        t_beg = round(t_beg, 3);
        t_end = round(t_end, 3);
        % Save the remaining segments
        o = ones(length(centers), 1);
        s = [i_preproc_trial*o i_chan*o centers' beginnings' endings' t_beg' t_end'];
        segments = [segments; s];
    end
end

% Snip the MEG data out

% Set up the time variables for the segment and whole trials
t_seg = -(segment_duration/2):(1/fsample):(segment_duration/2);
t_whole = round(data_tfr.time, 3);

% Set up Fieldtrip-style data struct for alpha
data_seg_alpha = [];
data_seg_alpha.label = data_preproc.label;
data_seg_alpha.time = {};
data_seg_alpha.trial = {};

% Set up FT-style data struct for TFR
data_seg_tfr = data_tfr;
data_seg_tfr.time = t_seg;
data_seg_tfr.powspctrm = nan(size(segments, 1),... % rpt_chan_freq_time
    length(data_tfr.label), ...
    length(data_tfr.freq), ...
    length(t_seg));

% Stick segments into the data structs
for i_segment = 1:size(segments, 1)
    % Info about each seg
    seg_info = segments(i_segment,:);
    preproc_trial_num = seg_info(1); % Number of the trial in data_preproc
    chan = seg_info(2);
    centersamp = seg_info(3);
    beginsamp = seg_info(4);
    endsamp = seg_info(5);
    t_beg = seg_info(6);
    t_end = seg_info(7);
    % Get the alpha-filtered data
    x = data_alpha.trial{preproc_trial_num}(:,beginsamp:endsamp);
    data_seg_alpha.trial{i_segment} = x;
    data_seg_alpha.time{i_segment} = t_seg;
    % Get the TFR data
    t_inx = (t_beg <= t_whole) & (t_whole <= t_end);
    x = data_tfr.powspctrm(preproc_trial_num,chan,:,t_inx);
    data_seg_tfr.powspctrm(i_segment, chan, :, :) = x;
end

% Average within each condition, b/c the TFR objects are huge (>2.2 GB)
cond_counts = nan(2,2);
cond_alpha = cell(2,2);
cond_tfr = cell(2,2);
% Get the RFT freq on the left side for each segment
preproc_trial_nums = segments(:,1);
preproc_left_rft_freqs = data_preproc.trialinfo(:,3);
seg_left_rft_freqs = preproc_left_rft_freqs(preproc_trial_nums);
rft_freqs = exp_params.tagged_freqs;
% Make the averages for each condition
for i_chan = 1:length(data_preproc.label)
    for i_rft_freq = 1:2
        % Select segments with the alpha peak in this (virtual) channel
        chan_trial_sel = segments(:,2) == i_chan;
        % Select segments with each RFT frequency at this chan
        % First get trials with this RFT frequency on the LEFT
        rft_trial_sel = seg_left_rft_freqs == rft_freqs(i_rft_freq);
        % If we're actually looking at the right virtual channel, invert it
        if i_chan == 2
            rft_trial_sel = ~rft_trial_sel
        end
        % Combine trial selections
        trial_sel = chan_trial_sel & rft_trial_sel;
        cond_counts(i_chan, i_rft_freq) = sum(trial_sel);
        % Average over alpha activity
        cfg = [];
        cfg.trials = trial_sel;
        cfg.channel = data_seg_alpha.label(i_chan);
        avg_alpha = ft_timelockanalysis(cfg, data_seg_alpha);
        cond_alpha{i_chan, i_rft_freq} = avg_alpha;
        % Average over TFRs
        avg_tfr = data_seg_tfr;
        avg_tfr.label = {'chan'};
        x = avg_tfr.powspctrm;
        x = x(trial_sel, i_chan, :, :);
        x = mean(x, 1);
        avg_tfr.powspctrm = x;
        cond_tfr{i_chan, i_rft_freq} = avg_tfr;
    end
end


fname = subject_info.meg{i_subject};
fn = [exp_dir 'alpha_peaks/' strrep(fname, '/', '_')];
save(fn, 'cond_counts', 'cond_alpha', 'cond_tfr', 'segments')


%% Plot the results
side_labels = {'left' 'right'};
x_lim = 0.3;
close all
i_cond = 1;
for i_chan = 1:2
    for i_rft_freq = 1:2
        % Plot the TFR
        width = 0.19;
        spacing = 0.05;
        lpos = spacing + (width + spacing) * (i_cond-1);
        subplot('position', [lpos, 0.35, width, 0.55])
        d = cond_tfr{i_chan, i_rft_freq};
        x = squeeze(d.powspctrm);
        % Divide by mean in each freq
        for i_freq = 1:length(d.freq)
            x(i_freq,:) = x(i_freq,:) / nanmean(x(i_freq,:));
        end
        imagesc(d.time, d.freq, x)
        set(gca, 'YDir', 'normal')
        xlim([-1 1] * x_lim)
        if i_cond == 1
            ylabel('Frequency (Hz)')
        end
        colorbar('SouthOutside')
        
        % Add a title
        title(sprintf('%s, %i Hz, n=%i', ...
            side_labels{i_chan}, ...
            rft_freqs(i_rft_freq), ...
            cond_counts(i_chan, i_rft_freq)))

        % Plot the alpha power
        subplot('position', [lpos, 0.13, width, 0.17])
        d = cond_alpha{i_chan, i_rft_freq};
        plot(d.time, d.avg, '-b')
        hold on
        plot(0, 0, '+k')
        hold off
        xlim([-1 1] * x_lim)
        
        if i_cond == 1
            ylabel('Amplitude (T)')
            xlabel('Time (s)')
        end
        
        i_cond = i_cond + 1;
    end
end

width = 25;
height = 10;
set(gcf,'units','centimeters')
set(gcf,'paperunits','centimeters')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize', [width height])
set(gcf,'paperposition',[0,0,width,height])
set(gcf, 'renderer', 'painters');

fname = subject_info.meg{i_subject};
fn = [exp_dir 'plots/alpha_peaks/' strrep(fname, '/', '_')];
print('-dpng', fn)

