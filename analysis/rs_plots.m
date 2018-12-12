%% Plot the target threshold over the course of the experiment

rs_setup
close all

clrs = [0 0.6 0; 0.6 0 0.6];

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    % Load the behavioral data
    fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
    behav = rs_behavior(fn);
    subplot(4,4,i_subject)
    hold on
    for s = {'left' 'right'}
        for f = [63 78]
            if f == 63
                clr = clrs(1,:);
            else
                clr = clrs(2,:);
            end
            if strcmp(s, 'left')
                clr = clr + (1 - max(clr)); % Left - lighter
            else
                clr = clr * 0.8; % Right - darker
            end
            inx = strcmp(behav.target_side,s) & (behav.target_side_freq==f);
            plot(behav.target_opacity(inx), '-', 'color', clr)
            ylim([0 0.4])
        end
    end
end

subplot(4,4,1)
ylabel('Target Opacity')
xlabel('Trial')

subplot(4,4,5)
text(1,1,'Left', 'color', [1 1 1] * 0.7)
text(2,1,'Right', 'color', [1 1 1] * 0.3)
xlim([1 3])
ylim([0 2])
axis off

subplot(4,4,5)
text(1,1, 'Left', 'color', [1 1 1] * 0.7)
text(2,1, 'Right', 'color', [1 1 1] * 0.3)
text(1,0, '63 Hz', 'color', clrs(1,:))
text(2,0, '78 Hz', 'color', clrs(2,:))
xlim([1 3])
ylim([0 2])
axis off

print('-dpng', [exp_dir 'plots/behav/opacity'])


%% Behavioral performance over course of expt

rs_setup

% Subject * Side (left/right) * Freq (63/78)
accuracy = nan(height(subject_info), 2, 2);

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
    behav = rs_behavior(fn);

    trials_per_block = size(behav,1) / 5;
    block_num = repelem((1:10)', trials_per_block / 2, 1);
    behav.block_num = block_num;

    % Only look at the the 'real' trials
    behav = behav(behav.target_t > 1.0, :);
    
    % Get the numbers of hits and FAs for each block
    hits = accumarray(behav.block_num, ...
        behav.hit, ...
        [], @mean);
    fas = accumarray(behav.block_num, ...
        behav.false_alarm, ...
        [], @mean);

    
    cols = [0 0.6 0; 0.6 0 0.6];

    %Plot the results
    close all
    figure('position', [500, 500, 500, 200])
    subplot(1,2,1)
    plot(1:10, hits, '-b')
    hold on
    plot([1 10], [0.5 0.5], '--b') % Expected performance
    plot(1:10, fas, '-r')

    % Make bars for different parts of the experiment
    rectangle('Position', [1 0.9 3.5 0.1], ...
        'EdgeColor', [1 1 1], 'FaceColor', 0.5 * [1 1 1])
    text(1.5, 0.95, 'Thresh', 'color', [1 1 1])

    rectangle('Position', [4.75 0.9 6 0.1], ...
        'EdgeColor', [1 1 1], 'FaceColor', 0.5 * [1 1 1])
    text(6.5, 0.95, 'Test', 'color', [1 1 1])

    hold off
    text(8, 0.8, 'Hit', 'color', 'b')
    text(8, 0.7, 'FA', 'color', 'r')
    xticks(1:10)
    xlim([1 10])
    ylim([0 1])
    xlabel('Block')
    ylabel('Proportion')

    % Plot accuracy in each condition in the main test phase
    subplot(1,2,2)
    means = nan(2,2);
    sides = {'left' 'right'};
    for i_side = 1:2
        for i_freq = 1:2
            mask = behav.block_num >= 5;
            mask = mask & behav.target_t > 1.5;
            mask = mask & strcmp(behav.target_side, sides{i_side});
            mask = mask & ...
                behav.target_side_freq == exp_params.tagged_freqs(i_freq);
            x = mean(behav.hit(mask));
            means(i_side, i_freq) = x;
            accuracy(i_subject, i_side, i_freq) = x;
        end
    end
    b = bar(means, 'EdgeColor', 'none');
    %hold on
    %plot([0.5 2.5], [0.5 0.5], '--k')
    %hold off
    ylim([0 1])
    xlim([0.5 2.5])
    ylabel('Proportion hits')
    xticklabels(sides)
    xlabel('Side')
    text(2, 0.9, '63 Hz', 'color', b(1).FaceColor)
    text(2, 0.8, '78 Hz', 'color', b(2).FaceColor)

    fname = subject_info.behav{i_subject};
    print('-dpng', [exp_dir 'plots/behav/' strrep(fname,'/','_')])

end


%% Plot the overall accuracy
figure('position', [500, 500, 200, 200])
means = squeeze(nanmean(accuracy, 1));
sterrs = squeeze(nanstd(accuracy, 1) ./ sqrt(sum(~isnan(accuracy), 1)));
b = bar(means, 'EdgeColor', 'none');
hold on
% Add error bars
for i_side = 1:2
    for i_freq = 1:2
        y = means(i_side, i_freq) + (sterrs(i_side, i_freq) * [-1 1]);
        x = i_side + ((i_freq - 1.5) * ([1 1] * 0.3));
        plot(x, y, '-k')
    end
end
hold off
ylim([0 1])
xlim([0.5 2.5])
ylabel('Proportion hits')
xticklabels(sides)
xlabel('Side')
text(2, 0.9, '63 Hz', 'color', b(1).FaceColor)
text(2, 0.8, '78 Hz', 'color', b(2).FaceColor)

print('-dpng', [exp_dir 'plots/behav/overall_accuracy'])


%% Look at artifact counts for each subject
% Plot over the course of the trial:
% - Trial start
% - Trial end
% - Density of artifacts at each timepoint
% This doesn't seem to reflect trialcounts after rejecting artifacts

clear variables
close all
rs_setup


%%% This doesn't reflect outcome of rs_preproc -- don't use it
%{
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    art_dir = [exp_dir 'artifacts/' fname];
    eye_art = load([art_dir '/eye']);
    photo_art = load([art_dir '/photodiode']);

    % Make objects to keep track of events
    fsample = 1000;
    t = -1:(1/fsample):6; % Time in seconds
    t_samp = t * fsample; % Time in samples
    trial_on = zeros(size(t));
    blink_on = zeros(size(t));
    photo_on = zeros(size(t));
    
    % Go through each recording
    for i_rec = block_info.main
        trialdef = load([exp_dir 'trialdef/' fname '/' num2str(i_rec)]);
        for i_trial = 1:size(trialdef.trl.(evt), 1) % Go through the trials
            start_samp = trialdef.trl.(evt)(i_trial, 1);
            if isnan(start_samp)
                continue
            end
            end_samp = trialdef.trl.(evt)(i_trial, 2);
            trial_dur = end_samp - start_samp;
            % Log how much of time this trial covers
            slice = ((start_samp + 1):end_samp) - start_samp;
            trial_on(slice) = trial_on(slice) + 1;
            % Select the artifacts in the current trial
            eye_art_rel = eye_art.eyes_artfctdef{i_rec} - start_samp;
            eye_curr = (eye_art_rel(:,1)>0) & (eye_art_rel(:,1)<trial_dur);
            eye_curr = find(eye_curr);
            photo_art_rel = photo_art.photo_artfctdef{i_rec} - start_samp;
            photo_curr = (photo_art_rel(:,1)>0) & (photo_art_rel(:,1)<trial_dur);
            photo_curr = find(photo_curr);
            % Log how much time those artifacts cover
            if ~isempty(eye_curr)
                for i = eye_curr'
                    this_art = eye_art_rel(i,:);
                    slice = this_art(1):this_art(2);
                    slice(slice > max(t_samp)) = [];
                    blink_on(slice) = blink_on(slice) + 1;
                end
            end
            if ~isempty(photo_curr)
                for i = photo_curr'
                    this_art = photo_art_rel(i,:);
                    slice = this_art(1):this_art(2);
                    slice(slice > max(t_samp)) = [];
                    photo_on(slice) = photo_on(slice) + 1;
                end
            end
        end
    end
    
    subplot(4, 4, i_subject)
    plot(t, trial_on, '-k')
    hold on
    plot(t, blink_on, '-g')
    plot(t, photo_on, '-b')
    plot(t, trial_on - blink_on - photo_on, '-r')
    hold off
    box off
    ylim([0 400])
    yticks([0 336])
    if strcmp(evt, 'trial')
        xlim([-1 5])
        xticks([-1 0 5])
    else
        xlim([-1 1])
        xticks([-1 0 1])
    end
end

subplot(4,4,1)
xlabel('Time (s)')
ylabel('Number of trials')

subplot(4,4,5)
text(1, 4, 'Trials overall', 'color', 'k')
text(1, 3, 'Photodiode artifact', 'color', 'b')
text(1, 2, 'Blink artifact', 'color', 'g')
text(1, 1, 'Trials after exclusion', 'color', 'r')
axis off
xlim([1 4])
ylim([0 5])

print('-dpng', [exp_dir 'plots/artifacts/trial_counts_' evt])
%}

close all
for segment_type = {'target' 'trial'}
    segment_type = segment_type{1};
    for i_subject = 1:height(subject_info)
        if subject_info.exclude(i_subject)
            continue
        end

        %fname = subject_info.meg{i_subject};
        %data = rs_preproc(fname, 'trial');

        % Load preprocessed data
        fname = subject_info.meg{i_subject};
        disp(fname)
        data = load([exp_dir 'preproc/' segment_type '/' fname '/preproc']);
        data = data.data;

        trial_len = cellfun(@(x) size(x,2), data.trial);
        nan_counts = zeros([1 max(trial_len)]);
        trial_counts = nan_counts;
        for i_trial = 1:length(data.trial)
            x = data.trial{i_trial}(1,:); % Look for NaNs in one channel
            nan_counts(1:length(x)) = nan_counts(1:length(x)) + isnan(x);
            trial_counts(1:length(x)) = trial_counts(1:length(x)) + ~isnan(x);
        end
        subplot(4,4,i_subject)
        plot(nan_counts, '-r')
        hold on
        plot(trial_counts, '-k')
        hold off
        xlim([0 max(trial_len)])
        ylim([0 336])
        xticks([])
        yticks([])
    end

    subplot(4,4,13)
    xlabel('Time (sample)')
    ylabel('Count')
    xticks([0 max(trial_len)])
    yticks([0 336])

    subplot(4,4,5)
    text(1, 2, 'Data after exclusion', 'color', 'k')
    text(1, 1, 'NaN', 'color', 'r')
    axis off
    xlim([1 4])
    ylim([0 5])

    print('-dpng', [exp_dir 'plots/artifacts/trial_counts_' segment_type])
end

%% Scalp topo of the SNR for the tagged frequencys
% SNR: Compare power in a narrow frequency to power at nearby frequencies

clear variables
rs_setup

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    close all

    spec = load([exp_dir 'spectra\' fname '\spectra']);
    spec = spec.freq_data;
    snr = load([exp_dir 'spectra\' fname '\snr']);
    snr = snr.snr;
    grad = load([exp_dir 'grad\' fname '\grad'], 'grad'); % To combine grads

    f_tag = exp_params.tagged_freqs; % Tagged frequencies

    i_plot = 1;
    for sensor_type = {'grad_cmb' 'mag'}
        sensor_type = sensor_type{1}; %#ok<FXSET>
        for i_freq = 1:length(f_tag)
            x = snr{i_freq};
            if strcmp(sensor_type, 'grad_cmb')
                x.grad = grad.grad;
                cfg = [];
                cfg.method = 'sum';
                x = ft_combineplanar(cfg, x);
            end
            subplot(3, length(f_tag), i_plot)
            title(sprintf('%s, %d Hz', ...
                sensor_type, f_tag(i_freq)), ...
                'Interpreter', 'none')
            if strcmp(sensor_type, 'mag')
                cfg = [];
                cfg.channel = chan.mag.names;
                x = ft_selectdata(cfg, x);
            end
            cfg = [];
            cfg.layout = chan.(sensor_type).layout;
            cfg.zlim = [1 max(x.powspctrm)];
            cfg.colorbar = 'no';
            cfg.style = 'straight';
            cfg.comment = 'no';
            cfg.shading = 'interp';
            cfg.markersymbol = '.';
            cfg.gridscale = 200;
            ft_topoplotTFR(cfg, x)
            c = colorbar();
            c.Label.String = 'SNR';
            i_plot = i_plot + 1;
            clear x
        end
    end

    % Plot the spectrum separately for each channel
    subplot(3,1,3)
    plot(spec.freq, 20 * log10(spec.powspctrm))
    hold on
    plot(f_tag, [1 1] * -595, '^r', 'linewidth', 3)
    hold off
    xlim([0 100])
    ylim([-600 -450])
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')

    print('-dpng', '-r300', [exp_dir 'plots/snr_topo/' strrep(fname,'/','_')])
end


%% Plot power at the tagged freqs time-locked to stimulus onset

% Instead of plotting power at tagged freq, plot instantaneous SNR at the
%  tagged frequency? Not at first.

clear variables
rs_setup

approx_eq = @(x,y) abs(x - y) < 0.1;

% Hold onto the data for all subject
% Subject * Freq * Time
overall_data = nan(height(subject_info), 46, 51);

close all
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    disp(fname)
    hf_fname = [exp_dir 'tfr/trial/' fname '/high.mat'];
    %finfo = dir(hf_fname);
    %disp(finfo.date)
    hf = load(hf_fname);
    freq_data = hf.high_freq_data;
    clear hf;

    % Select channels in this roi
    chan_sel = ismember(freq_data.label, snr_roi);
    
    % Exclude the trials that are all NaNs
    % Those appeared while computing TFRs of trials with NaNs
    nan_trial = all(all(all(isnan(freq_data.powspctrm), 2), 3), 4);
    fprintf('%i NaN; %i numeric \n', ...
        sum(nan_trial), sum(~nan_trial))

    cfg = [];
    cfg.trials = find(~nan_trial);
    cfg.avgoverrpt = 'yes';
    cfg.nanmean = 'yes';
    freq_data = ft_selectdata(cfg, freq_data);
    x = mean(freq_data.powspctrm(chan_sel,:,:), 1);
    overall_data(i_subject,:,:) = x;
    clear x
    
    % Plot of frequency tagging response from trial onset
    subplot(2,1,1)
    cfg = [];
    cfg.channel = freq_data.label(chan_sel);
    cfg.baseline = [-0.5 -0.1];
    cfg.baselinetype = 'relative';
    cfg.title = ' ';
    cfg.colorbar = 'no';
    ft_singleplotTFR(cfg, freq_data);
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    hold on
    scatter([0 0], exp_params.tagged_freqs, 'w>', 'filled')
    title(sprintf('%i trials', sum(~nan_trial)))
    hold off
    % zlabel('Power')

    % Plot of power at the tagged frequencies over time
    subplot(2,1,2)
    cfg = [];
    cfg.baseline = [-0.5 -0.1];
    cfg.baselinetype = 'relative';
    freq_data = ft_freqbaseline(cfg, freq_data);

    clr = [0.9 0 0.9; 0 0.5 1];
    for i_freq = 1:2
        freq_inx = approx_eq(freq_data.freq, ...
            exp_params.tagged_freqs(i_freq));
        x = freq_data.powspctrm(chan_sel, freq_inx, :);
        x = squeeze(nanmean(x, 1)); % Average over channels
        plot(freq_data.time, x, 'LineWidth', 2, 'color', clr(i_freq,:));
        ypos = max(x);
        text(-0.4, ypos, ...
            [num2str(exp_params.tagged_freqs(i_freq)) ' Hz'], ...
            'color', clr(i_freq,:))
        hold on
    end
    xlabel('Time (s)')
    ylabel('Power (relative change)')
    plot([-0.5 1.5], [1 1], '--k')
    hold off

    print('-dpng', '-r300', ...
        [exp_dir 'plots/stim_onset/' strrep(fname, '/', '_')])
end

% Plot the averages over all subjects
freq_data.powspctrm = nanmean(overall_data, 1);
freq_data.label = {'averaged'};
subplot(2,1,1)
cfg = [];
% cfg.channel = union(roi{:});
cfg.baseline = [-0.5 -0.1];
cfg.baselinetype = 'relative';
cfg.title = ' ';
cfg.colorbar = 'no';
ft_singleplotTFR(cfg, freq_data);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
hold on
scatter([0 0], exp_params.tagged_freqs, 'w>', 'filled')
hold off

% % Plot of power at the tagged frequencies over time
subplot(2,1,2)
cfg = [];
cfg.baseline = [-0.5 -0.1];
cfg.baselinetype = 'relative';
freq_data = ft_freqbaseline(cfg, freq_data);

clr = [0.9 0 0.9; 0 0.5 1];
for i_freq = 1:2
    freq_inx = approx_eq(freq_data.freq, ...
        exp_params.tagged_freqs(i_freq));
    x = freq_data.powspctrm(:, freq_inx, :);
    x = squeeze(nanmean(x, 1)); % Average over channels
    plot(freq_data.time, x, 'LineWidth', 2, 'color', clr(i_freq,:));
    ypos = max(x);
    text(-0.4, ypos, ...
        [num2str(exp_params.tagged_freqs(i_freq)) ' Hz'], ...
        'color', clr(i_freq,:))
    hold on
end
xlabel('Time (s)')
ylabel('Power (relative change)')
plot([-0.5 1.5], [1 1], '--k')
hold off

print('-dpng', '-r300', ...
    [exp_dir 'plots/stim_onset/overall'])


%% Plot difference in HF power between hits and misses

clear variables
rs_setup

roi_type = 'anatomical'; % 'anatomical' or 'functional'

approx_eq = @(x,y) abs(x - y) < 0.1;

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    disp(fname)
    hf_fname = [exp_dir 'tfr/target/' fname '/high.mat'];
    hf = load(hf_fname);
    freq_data = hf.high_freq_data;
    clear hf;

    if strcmp(roi_type, 'anatomical')
        % Anatomical ROI
        roi = {occip_roi occip_roi}; % Duplicate for 2 frequencies
    elseif strcmp(roi_type, 'functional')
        roi = rs_roi(fname, 1);     
    end
    
    % Load the behavioral data
    fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
    behav = rs_behavior(fn);
    behav = behav(225:end, :); % main blocks of trials

    % Load the trialdefs to get hits/misses
    hit = [];
    missing = [];
    for i_rec = block_info.main
        fn = [exp_dir 'trialdef/' fname '/' num2str(i_rec) '.mat'];
        trialdef = load(fn);
        missing = [missing isnan(trialdef.trl.target(:,1))']; % Trials w/o a target
        hit = [hit trialdef.trl.trial(:, 4)'];
    end
    hit = hit(~missing);
    
    % Select hits/misses for each target frequency
    % Along the way, exclude the trials that are all NaNs
    % (How did that happen? Check those trials)
    nan_trial = all(all(all(isnan(freq_data.powspctrm), 2), 3), 4)';
    
    % Plot the overall TFR
    clr = [0.9 0 0.9; 0 0.5 1];
    close all
    for i_freq = 1:2
        targ_freq = exp_params.tagged_freqs(i_freq);
        freq_inx = (behav.target_side_freq(~missing) == targ_freq)';
        cfg = [];
        cfg.avgoverrpt = 'yes';
        % Hits
        cfg.trials = (hit == 1) & ~nan_trial & freq_inx;
        hit_data = ft_selectdata(cfg, freq_data);
        % Misses
        cfg.trials = (hit == 0) & ~nan_trial & freq_inx;
        miss_data = ft_selectdata(cfg, freq_data);
        % Diff
        diff_data = hit_data;
        diff_data.powspctrm = hit_data.powspctrm - miss_data.powspctrm;
    
        % Plot of frequency tagging response from trial onset
        cfg = [];
        cfg.channel = union(roi{:});
        cfg.title = ' ';
        cfg.colorbar = 'no';

        figure(1) % Overall TFR
        subplot(2, 3, 1 + (3 * (i_freq - 1)))
        ft_singleplotTFR(cfg, hit_data);
        subplot(2, 3, 2 + (3 * (i_freq - 1)))
        ft_singleplotTFR(cfg, miss_data);
        subplot(2, 3, 3 + (3 * (i_freq - 1)))
        ft_singleplotTFR(cfg, diff_data);
        
        figure(2) % Time-course of power at tagged frequencies
        hold on
        cfg = [];
        cfg.frequency = [-0.5 0.5] + targ_freq;
        cfg.avgoverfreq = 'yes';
        cfg.channel = union(roi{:});
        cfg.avgoverchan = 'yes';
        cfg.nanmean = 'yes';
        hit_data = ft_selectdata(cfg, hit_data);
        miss_data = ft_selectdata(cfg, miss_data);
        diff_data = hit_data;
        diff_data.powspctrm = hit_data.powspctrm - miss_data.powspctrm;
        plot(diff_data.time, squeeze(diff_data.powspctrm), ...
            'color', clr(i_freq,:))
        hold off
    end
    
    figure(1)
    subplot(2,3,1)
    title('Hit')
    subplot(2,3,2)
    title('Miss')
    subplot(2,3,3)
    title('Difference')
    
    subplot(2,3,1)
    text(-1, 70, '63 Hz', 'Rotation', 90)
    ylabel('Frequency (Hz)')
    subplot(2,3,4)
    text(-1, 70, '78 Hz', 'Rotation', 90)
    ylabel('Frequency (Hz)')
    
    subplot(2, 3, 4)
    xlabel('Time (s)')

    print('-dpng', '-r300', ...
        [exp_dir 'plots/hit_miss/tfr_' strrep(fname, '/', '_')])

    figure(2)
    xlabel('Time (s)')
    ylabel('Power (Hit - Miss)')
    hold on
    plot([-0.5 0.5], [0 0], '--k')
    hold off
    
    print('-dpng', '-r300', ...
        [exp_dir 'plots/hit_miss/timecourse_' strrep(fname, '/', '_')])
       
end


%% Results of the HF power regression

clear variables
rs_setup
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    x = load([exp_dir 'tfr/target/' fname '/high_power_acc_stats']);
    
    % Channels in the ROI
    keep_chans = ismember(x.label, occip_roi);
    
    % Extract the regression coefficients for hit/miss predicting HF power
    coef = cellfun(@(c) c{'Hit', 'Estimate'}, x.stats);
    avg_coef = mean(coef(:,:,keep_chans), 3);
    % Extract p-values
    pval = cellfun(@(c) c{'Hit', 'pValue'}, x.stats);
    avg_pval = geomean(pval(:,:,keep_chans), 3); %%% <-------------- FIX ME
    
    % Plot it
    close all
    
    figure(1)
    imagesc(x.time, x.freq, avg_coef');
    set(gca, 'YDir', 'normal')
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    c = colorbar;
    ylabel(c, 'Beta')
    colormap('parula')
    
    figure(2)
    imagesc(x.time, x.freq, log10(avg_pval)');
    set(gca, 'YDir', 'normal')
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    caxis([min(min(log10(avg_pval))) 0])
    c = colorbar;
    ylabel(c, 'log10 p-value')
    colormap(flipud(bone))
    
    print(FILL_IN)
end
    
% Plot the averages over subjects
    
    
    
%% Does hit-rate differ by LF phase?

clear variables
rs_setup

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    x = load([exp_dir 'tfr/target/' fname '/low_power_acc_stats.mat']);
    
    % Anatomical ROI
    keep_chans = ismember(x.label, occip_roi);
    
    hit_rate = x.hit_rate(:, keep_chans, :);
    hit_rate = squeeze(mean(hit_rate, 2));
    
    bin_start = linspace(-pi, pi, x.n_bins);
    imagesc(bin_start, x.freq, hit_rate')
    
    xticks([-pi 0 pi])
    xticklabels({'-\pi', '0', '\pi'})
    xlabel('Phase')
    ylabel('Frequency (Hz)')
    set(gca, 'YDir', 'normal')
    
    cbh = colorbar('v');
    set(cbh,'YTick', [min(min(hit_rate)) max(max(hit_rate))])
    ylabel(cbh, 'Prop. hits')
  
    print('-dpng', '-r300', ...
        [exp_dir 'plots/hit_miss/lf_' strrep(fname, '/', '_')])
end


%% Cross-correlation of power at the tagged frequencies

% Load the data
clear variables
rs_setup
x_overall = nan(height(subject_info), 51); % Subject * Time
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    disp(fname)
    x = load([exp_dir 'xcorr/' fname '/x']);
    keepchans = ismember(x.label, occip_roi);
    x_subj = squeeze(mean(nanmean(x.x(:,keepchans,:), 1), 2));
    x_overall(i_subject,:) = x_subj;
end

close all
figure('position', [500, 500, 500, 200])

% Plot the cross-correlations
subplot(1,2,1)
plot(x.time, x_overall, '-'); %, 'color', [0.7 0.7 1]);
hold on
plot(x.time, nanmean(x_overall, 1), '-b', 'LineWidth', 2.5)
plot([-1 1], [0 0], '--k')
hold off
xlabel('Lag (s)')
ylabel('Correlation')

% Compute and plot the FFT of the cross-correlation
% uk.mathworks.com/help/matlab/examples/fft-for-spectral-analysis.html
subplot(1,2,2)
nfft = 2^6;
sample_per = mean(diff(x.time));
Fs = 1 / sample_per;
f = (1/sample_per) * (0:(nfft / 2)) / nfft;
y = fft(x_overall, nfft, 2);
Pyy = 1 / (nfft * Fs) * abs(y(:,1:nfft/2+1)) .^ 2; % Power spectrum

plot(f, db(Pyy, 'power'), '-'); %, 'color', [0.7 0.7 1]);
hold on
plot(f, db(nanmean(Pyy, 1), 'power'), '-b', 'LineWidth', 2.5)
hold off
xlabel('Frequency (Hz)')
ylabel('Power (dB/Hz)')
xlim([0 12])

print('-dpng', '-r300', ...
    [exp_dir 'plots/xcorr'])


%% Plot the autocorrelation of RTF power and its spectrum











%%%%%%
%%%%%% Beyond this point are scripts from the old expt.




%% plot theta phase at target time for hits/misses

clear variables
close all
rs_setup

i_subject = 2;
fname = subject_info.meg{i_subject};
d = load([exp_dir 'theta\' fname '\phase_targets']);
froi = load([exp_dir 'froi\' fname '\froi']);
fn = [exp_dir 'behav\' subject_info.behav{i_subject} '_main.csv'];
[behav_trials, behav_targets] = rs_behavior(fn);

phi = ft_appenddata([], d.d_blocks{:});
clear d

cfg = [];
cfg.channel = froi.froi.grad; % Choose the grad chans with best SNR
cfg.avgoverchan = 'no';
cfg.latency = [-0.0001 0.0001]; % Time of the target onset
cfg.avgovertime = 'yes';
phi = ft_selectdata(cfg, phi);
phi = cell2mat(phi.trial);

% Plot - Histogram of theta phase, separately for hits and misses
% This might not be very informative, since we don't have a ton of data
%mathworks.com/help/matlab/creating_plots/polar-axes-grid-lines-and-labels.html
n_bins = 10;
for i_chan = 1:size(phi, 1)
    subplot(3, 3, i_chan)
    bins = linspace(-pi, pi, n_bins);
    h = polarhistogram(phi(i_chan, ~behav_targets.hit), bins, ...
        'FaceColor', 'red', 'EdgeColor', 'red', 'LineWidth', 2, ...
        'Normalization', 'probability');
    hold on
    polarhistogram(phi(i_chan, behav_targets.hit), bins, ...
        'FaceColor', 'blue', 'EdgeColor', 'blue',  'LineWidth', 2, ...
        'Normalization', 'probability')
    hold off
    thetaticks([])
    rticks([])
end

% Plot - Sliding window of Prob(hit) as a function of theta phase
% As in Fiebelkorn et al (2018, Neuron)