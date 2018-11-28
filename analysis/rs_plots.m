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


%% Behavioral performance over course of expt

rs_setup
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
        end
    end
    b = bar(means, 'EdgeColor', 'none');
    hold on
    plot([0.5 2.5], [0.5 0.5], '--k')
    hold off
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


%% Plot power at the tagged freqs time-locked to stimulus onset

% Instead of plotting power at tagged freq, plot instantaneous SNR at the
%  tagged frequency? Not at first.

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
    hf_fname = [exp_dir 'tfr/trial/' fname '/high.mat'];
    finfo = dir(hf_fname);
    disp(finfo.date)
    hf = load(hf_fname);
    freq_data = hf.high_freq_data;
    clear hf;

    if strcmp(roi_type, 'anatomical')
        % Anatomical ROI
        roi = [2012 2013 ... % occipital gradiometers
            2022 2023 ...
            2032 2033 ...
            2042 4043 ...
            2112 2113];
        roi = cellfun(@(n) ['MEG' num2str(n)], num2cell(roi), ...
            'UniformOutput', false);
        roi = {roi roi}; % Duplicate for 2 frequencies
    elseif strcmp(roi_type, 'functional')
        roi = rs_roi(fname, 1);     
    end
    
    % Exclude the trials that are all NaNs
    % (How did that happen? Check those trials)
    nan_trial = all(all(all(isnan(freq_data.powspctrm), 2), 3), 4);
    cfg = [];
    cfg.trials = find(~nan_trial);
    cfg.avgoverrpt = 'yes';
    freq_data = ft_selectdata(cfg, freq_data);
    
    close all
    % Plot of frequency tagging response from trial onset
    subplot(2,1,1)
    cfg = [];
    cfg.channel = union(roi{:});
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
    % zlabel('Power')

    % Plot of power at the tagged frequencies over time
    subplot(2,1,2)
    cfg = [];
    cfg.baseline = [-0.5 -0.1];
    cfg.baselinetype = 'relative';
    freq_data = ft_freqbaseline(cfg, freq_data);

    clr = [0.9 0 0.9; 0 0.5 1];
    for i_freq = 1:2
        chan_sel = ismember(freq_data.label, roi{i_freq});
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


%% Plot difference between hits and misses


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
        roi = [2012 2013 ... % occipital gradiometers
            2022 2023 ...
            2032 2033 ...
            2042 4043 ...
            2112 2113];
        roi = cellfun(@(n) ['MEG' num2str(n)], num2cell(roi), ...
            'UniformOutput', false);
        roi = {roi roi}; % Duplicate for 2 frequencies
    elseif strcmp(roi_type, 'functional')
        roi = rs_roi(fname, 1);     
    end
    
    % Load the behavioral data
    fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
    behav = rs_behavior(fn);
    behav = behav(225:end, :); % main blocks of trials
    assert(all(behav.hit == freq_data.trialinfo))
    continue

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
        cfg.trials = hit & ~nan_trial & freq_inx;
        hit_data = ft_selectdata(cfg, freq_data);
        % Misses
        cfg.trials = ~hit & ~nan_trial & freq_inx;
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





%% Check eye-tracking
% Don't save plots, just go through and plot a lot of examples

clear variables
close all
rs_setup
i_subject = 4;

fname = subject_info.meg{i_subject};

data_by_block = cell(size(block_info.all));

for i_block = 1:length(block_info.main)
    n_block = block_info.main(i_block);
    
    % Read in the trial definition
    fn = [exp_dir 'trialdef\' fname '\trials_' num2str(n_block) '.mat'];
    if ~exist(fn, 'file')
        warning('No trialdef for sub %s, block %d', fname, n_block)
        continue
    end

    % Preprocess the data
    cfg = load(fn);
    cfg = cfg.cfg;
    cfg.channel = {'BIO001' 'BIO002' 'MISC001' 'MISC002'};
    d = ft_preprocessing(cfg);

    data_by_block{i_block} = d;
end
d = data_by_block;
clear data_by_block

i_block = 1;

cfg = [];
cfg.channel = 'BIO*';
eog = ft_selectdata(cfg, d{i_block});

cfg.channel = 'MISC*';
eyelink = ft_selectdata(cfg, d{i_block});


for i_trial = 1:length(eog.trial)
    subplot(2,1,1)
    plot(eog.time{i_trial}, eog.trial{i_trial})
    subplot(2,1,2)
    plot(eyelink.time{i_trial}, eyelink.trial{i_trial})
    input('')
end









%%%%%%
%%%%%% Beyond this point are scripts from the old expt.




%% Split trials into hits/misses with the target at each freq
data_by_cond = cell(2,2); % Freq x (Hit/Miss)
hit_vals = [true false];
for i_hit = 1:length(hit_vals)
    hit_sel = behav_targets.hit == hit_vals(i_hit);
    for i_freq = 1:length(power)
        freq_sel = behav_targets.freq == d.tagged_freqs(i_freq);
        cfg = [];
        cfg.trials = hit_sel & freq_sel;
        cfg.avgoverrpt = 'no';
        cfg.latency = [-0.0001 0.0001]; % Time of the target onset
        cfg.avgovertime = 'yes';
        cfg.channel = froi.froi.grad; % Select the best gradiometers
        cfg.avgoverchan = 'yes';
        data_by_cond{i_freq,i_hit} = ft_selectdata(cfg, power{i_freq});
    end
end

% Make boxplots
d = []; % Data vector
g = []; % Grouping vector
counter = 1;
for i_freq = 1:length(tagged_freqs)
    for i_hit = 1:length(hit_vals)
        x = cell2mat(data_by_cond{i_freq,i_hit}.trial);
        d = [d x];
        g = [g counter*ones(size(x))];
        counter = counter + 1;
    end
end

h = boxplot(d, g, ...
    'Positions', [1 2 4 5], 'Colors', 'rbrb');
set(h, ...
    'LineStyle', '-', ...
    'LineWidth', 2)

xticks([1.5 4.5])
xticklabels(tagged_freqs)
xlabel('Tagged frequency')
ylabel('Amplitude')

text(4.5, max(d) * 1.15, 'Hit', 'Color', 'r')
text(4.5, max(d) * 1.05, 'Miss', 'Color', 'b')
ylim([min(d) * 0.8, max(d) * 1.2])


%% Plot timecourse of amplitude at tagged freq

data_by_cond = cell(2,2); % Freq x (Hit/Miss)
hit_vals = [true false];
for i_hit = 1:length(hit_vals)
    hit_sel = behav_targets.hit == hit_vals(i_hit);
    for i_freq = 1:length(power)
        freq_sel = behav_targets.freq == tagged_freqs(i_freq);
        cfg = [];
        cfg.trials = hit_sel & freq_sel;
        cfg.avgoverrpt = 'no';
        cfg.channel = froi.froi.grad; % Select the best gradiometers
        cfg.avgoverchan = 'no';
        data_by_cond{i_freq,i_hit} = ft_selectdata(cfg, power{i_freq});
    end
end

matrify = @(i_freq, i_hit) ...
    squeeze(nanmean(cat(3, data_by_cond{i_freq,i_hit}.trial{:}), 1));
sterr = @(x) std(x,1) / sqrt(size(x,1));

for i_freq = 1:length(tagged_freqs)
    subplot(length(tagged_freqs), 1, i_freq)
    t = data_by_cond{1,1}.time{1};
    m = matrify(i_freq, 1);
    fill([t fliplr(t)], ...
        [mean(m,2) + sterr(m); fliplr(mean(m,2) + sterr(m))]', ...
        'r')
    hold on
    plot(data_by_cond{1,1}.time{1}, ...
        nanmean(m, 2), '-r')
    plot(data_by_cond{1,1}.time{1}, ...
        nanmean(matrify(i_freq,2), 2), '-b')
    hold off
end

%%% Why is this so jagged?

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