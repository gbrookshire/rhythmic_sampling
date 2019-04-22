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
% Look at the number of NaNs in each sample over the trials

clear variables
close all
rs_setup

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


%% Plot the spectra of the raw signals

clear variables
rs_setup

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject}; 
    spec = load([exp_dir 'spectra/' fname '/spectra']);
    spec = spec.freq_data;

    % Average over occipital channels
    chans = [2042 2032 2112];
    chans = [chans (chans + 1)];
    chans = cellfun(@(n) ['MEG' num2str(n)],...
        num2cell(chans),...
        'UniformOutput', false);
    cfg = [];
    cfg.channel = chans;
    cfg.avgoverchan = 'yes';
    spec = ft_selectdata(cfg, spec);

    p = spec.powspctrm;
    %disp(size(p))
    %p = p / sum(p); % Normalize
    s(i_subject,:) = p;
end


close all
figure('position', [50, 50, 200, 200])
plot(spec.freq, db(s, 'power'), '-', 'color', [1 1 1] * 0.6)
hold on
plot(spec.freq, db(nanmean(s, 1), 'power'), '-b', 'LineWidth', 1.5)

xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0 100])
box('off')
%plot(exp_params.tagged_freqs, [-38 -38], '^r')
hold off
print('-depsc', [exp_dir 'plots/spectra'])


%% Scalp topo of the SNR for the tagged frequencys
% SNR: Compare power in a narrow frequency to power at nearby frequencies

clear variables
rs_setup

% Object to hold topographies across subjects
topos = nan(height(subject_info), 204, 2); % Subject x Chan x Freq

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    close all

    spec = load([exp_dir 'spectra/' fname '/spectra']);
    spec = spec.freq_data;
    snr = load([exp_dir 'spectra/' fname '/snr']);
    snr = snr.snr;
    grad = load([exp_dir 'grad/' fname '/grad'], 'grad'); % To combine grads

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
                topos(i_subject,:,i_freq) = x.powspctrm; % Save the powspect
                topo_struct = x;
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

    %print('-dpng', '-r300', [exp_dir 'plots/snr_topo/' strrep(fname,'/','_')])
end

% Plot the average topography across subjects & frequencies
close all
figure('position', [200, 200, 200, 200])
topo_struct.powspctrm = squeeze(mean(nanmean(x, 1), 2));
cfg = [];
cfg.layout = chan.(sensor_type).layout;
cfg.zlim = [1 max(topo_struct.powspctrm)];
cfg.colorbar = 'no';
cfg.style = 'straight';
cfg.comment = 'no';
cfg.shading = 'interp';
cfg.markersymbol = '.';
cfg.gridscale = 200;
ft_topoplotTFR(cfg, topo_struct)
c = colorbar();
c.Label.String = 'SNR';
print('-dpng', '-r300', [exp_dir 'plots/snr_topo/avg'])

%% Plot the RESS spatial filters 

clear variables
rs_setup

ress_maps = cell([1 height(subject_info)]);
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    grad = load([exp_dir 'grad/' fname '/grad'], 'grad');
    data_ress = load([exp_dir 'ress/' fname '/ress']);
    ress_maps{i_subject} = data_ress.ress_maps;
end

% Load list of channel labels
labels = load([exp_dir 'spectra/' subject_info.meg{1} '/spectra']);
labels = labels.freq_data.label;

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    subj_map = ress_maps{i_subject};
    i_plot = 1;
    for freq = exp_params.tagged_freqs
        for side = {'left' 'right'}
            maps = subj_map.(side{1}).(['f' num2str(freq)]).maps(:,1);
            % Make data structure to show maps
            d_maps = [];
            d_maps.label = labels;
            d_maps.time = 1;
            d_maps.avg = real(maps);
            d_maps.dimord = 'chan_time';
            d_maps.grad = grad.grad;
            % Combine planar gradiometers
            cfg = [];
            cfg.method = 'sum';
            d_maps = ft_combineplanar(cfg, d_maps);

            % Plot it
            subplot(2,2,i_plot)
            cfg = [];
            cfg.marker = 'on';
            cfg.markersymbol = 'o';
            cfg.markersize = 1;
            cfg.comment = 'no';
            cfg.style = 'straight';
            cfg.layout = chan.grad_cmb.layout;
            cfg.gridscale = 200;
            ft_topoplotER(cfg, d_maps)
            title(sprintf('%s, %i Hz', side{1}, freq))
            i_plot = i_plot + 1;
        end
    end
    print('-dpng', '-r300', ...
        [exp_dir 'plots/ress_maps/' ...
        strrep(fname, '/', '_')])
    
end

% Plot filters averaged over subjects
ress_maps(cellfun(@isempty, ress_maps)) = []; % Delete empty cells
i_plot = 1;
for side = {'left' 'right'}
    
    m = cellfun(@(c) c.(side{1}).f63.maps(:,1), ...
        ress_maps, 'UniformOutput', false);
    m = [m{:}]; % Concatenate into a matrix

    % Make data structure to show maps
    d_maps = [];
    d_maps.label = labels;
    d_maps.time = 1:size(m,2); % Actually not time, but subject
    d_maps.avg = real(m);
    d_maps.dimord = 'chan_time';
    d_maps.grad = grad.grad;
    % Combine planar gradiometers
    cfg = [];
    cfg.method = 'sum';
    d_maps = ft_combineplanar(cfg, d_maps);

    % Plot it
    subplot(1,2,i_plot)
    cfg = [];
    cfg.marker = 'on';
    cfg.markersymbol = 'o';
    cfg.markersize = 1;
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.layout = chan.grad_cmb.layout;
    cfg.gridscale = 200;
    ft_topoplotER(cfg, d_maps)
    i_plot = i_plot + 1;
    title([side{1} ' stim'])
end
print('-dpng', '-r300', [exp_dir 'plots/ress_maps/avg'])


%% Plot the spectra after RESS spatial filters using 63 Hz
% Requires a lot of space & moving data around - do on the cluster
clear variables
rs_setup

spectra = cell([1 height(subject_info)]); % spectra for all subjects

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    close all
    fname = subject_info.meg{i_subject};
    data_ress = rs_preproc_ress(i_subject, 'trial');

    % Compute the spectrum
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.taper = 'hanning';
    cfg.pad = 'nextpow2';
    cfg.padtype = 'zero';
    cfg.polyremoval = 1; % Remove linear trends
    cfg.keeptrials = 'yes';
    spec = ft_freqanalysis(cfg, data_ress);
    spectra{i_subject} = spec;

    i_plot = 1;
    for freq = exp_params.tagged_freqs
        for side = {'left' 'right'} % Which stim side
            % Average over trials
            cfg = [];
            if strcmp(side{1}, 'left')
                cfg.trials = spec.trialinfo(:,3) == freq;
            else
                cfg.trials = spec.trialinfo(:,3) ~= freq;
            end
            cfg.avgoverrpt = 'yes';
            cfg.channel = side;
            data_sub = ft_selectdata(cfg, spec);
            % Plot it
            subplot(2,2,i_plot)
            plot(data_sub.freq, db(data_sub.powspctrm))
            xlim([0 100])
            title(sprintf('%i Hz, %s', freq, side{1}))
            xlabel('Frequency (Hz)')
            ylabel('Power (dB)')

            i_plot = i_plot + 1;
        end
    end
    print('-dpng', '-r300', ...
        [exp_dir 'plots/ress_maps/spec_63hzMap_' strrep(fname,'/','_')])
end

% Avg over subjects

close all
i_plot = 1;
for freq = exp_params.tagged_freqs
    for side = {'left' 'right'} % Which stim side
        % Combine data for all subjects
        spec_cond = nan([height(subject_info) length(spec.freq)]);
        for i_subject = 1:height(subject_info)
            if subject_info.exclude(i_subject) 
                continue
            end
            % Average over trials
            cfg = [];
            if strcmp(side{1}, 'left')
                cfg.trials = spectra{i_subject}.trialinfo(:,3) == freq;
            else
                cfg.trials = spectra{i_subject}.trialinfo(:,3) ~= freq;
            end
            cfg.avgoverrpt = 'yes';
            cfg.channel = side;
            data_sub = ft_selectdata(cfg, spectra{i_subject});
            % Normalize by area under the curve
            pwr = data_sub.powspctrm;
            pwr = pwr ./ sum(pwr);
            spec_cond(i_subject,:) = pwr;
        end
        spec_mean = nanmean(spec_cond, 1);
    
        [~, ~, spec_ci] = ttest(spec_cond);
        % Plot it
        subplot(2,2,i_plot)
        fill([spec.freq fliplr(spec.freq)], ...
            db([spec_ci(1,:), fliplr(spec_ci(2,:))]), ...
            'b', 'FaceALpha', 0.4, 'EdgeColor', 'none')
        hold on
        plot(spec.freq, db(spec_mean), '-b')
        hold off
        xlim([0 100])
        title(sprintf('%i Hz, %s', freq, side{1}))
        xlabel('Frequency (Hz)')
        ylabel('Power (dB)')
        
        i_plot = i_plot + 1;
        
    end
end

print('-dpng', '-r300', ...
    [exp_dir 'plots/ress_maps/spec_63hzMap_overall'])


%% Plot power at the tagged freqs time-locked to stimulus onset

clear variables
rs_setup

win_size = 0.1;
win_str = sprintf('win_%.1fs', win_size);


approx_eq = @(x,y) abs(x - y) < 0.1;

% Hold onto the data for all subject
% Subject * Freq * Time
overall_data = nan(height(subject_info), 46, 501);

close all
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    disp(fname)
    hf_fname = [exp_dir 'tfr/' win_str '/trial/' fname '/high.mat'];
    %finfo = dir(hf_fname);
    %disp(finfo.date)
    hf = load(hf_fname);
    freq_data = hf.high_freq_data;
    clear hf;
    
    % Exclude the trials that are all NaNs
    % Those appeared while computing TFRs of trials with NaNs
    nan_trial = all(all(all(isnan(freq_data.powspctrm), 2), 3), 4);
    fprintf('%i NaN; %i numeric \n', ...
        sum(nan_trial), sum(~nan_trial))

    % Average over subjects and channels
    cfg = [];
    cfg.trials = find(~nan_trial);
    cfg.avgoverrpt = 'yes';
    cfg.nanmean = 'yes';
    cfg.avgoverchan = 'yes';
    freq_data = ft_selectdata(cfg, freq_data);
    overall_data(i_subject,:,:) = freq_data.powspctrm;
    
    % Plot of frequency tagging response from trial onset
    subplot(2,1,1)
    cfg = [];
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
    plot([-0.5 4.5], [1 1], '--k')
    xlim([-0.5 4.5])
    xticks(0:4)
    hold off

    print('-dpng', '-r300', ...
        [exp_dir 'plots/stim_onset/' win_str '/' strrep(fname, '/', '_')])
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
plot([-0.5 4.5], [1 1], '--k')
xlim([-0.5 4.5])
xticks(0:4)
hold off

print('-dpng', '-r300', ...
    [exp_dir 'plots/stim_onset/' win_str '/overall'])


%% Simulated data - example trial

d = rs_simulate_flicker();

subplot(2, 1, 1)
plot(d.time{1}, d.trial{1})
xlabel('Time (s)')
ylabel('Amplitude')

subplot(2, 1, 2)
plot(d.time{1}, d.trial{1})
xlabel('Time (s)')
ylabel('Amplitude')
xlim([2 3])

print('-dpng', '-r300', ...
    [exp_dir 'plots/simulated/example'])


%% Plot spectra the envelope of tagged frequencies

clear variables
rs_setup

% Spectra for all subjects for each cond. Subj x CarrierFreq x ModFreq
all_spectra = [];

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    data_ress = load([exp_dir 'tfr/trial/' fname '/spect']);
    data_ress = data_ress.data_ress;

    close all
    i_cond = 1;
    for freq = exp_params.tagged_freqs
        for side = {'left' 'right'}
            
            % Get subset of the data
            if strcmp(side, 'left')
                trial_inx = data_ress.left.trialinfo(:,3) == freq;
            else
                trial_inx = data_ress.left.trialinfo(:,3) ~= freq;
            end
            cfg = [];
            cfg.trials = trial_inx;
            d = ft_selectdata(cfg, data_ress.(side{1}));
            
            % Plot
            subplot(2, 2, i_cond)
            car_freq = str2double(d.label);
            mod_freq = d.freq;
            spectra = d.powspctrm;
            spectra = squeeze(mean(spectra, 1)); % Avg over trials
            % Normalize across frequencies
            for i_row = 1:size(spectra, 1)
                % Dividing by the total area under the FFT curve
                spectra(i_row,:) = spectra(i_row,:) / sum(spectra(i_row,:));
            end
            all_spectra.(side{1}).(['f' num2str(freq)])(i_subject,:,:) = spectra;
            %spectra = db(spectra, 'power'); % Convert to dB
            imagesc(mod_freq, car_freq, spectra)
            set(gca,'YDir','normal')
            title(sprintf('%s, %i Hz', side{1}, freq))
            i_cond = i_cond + 1;
        end
    end
    subplot(2, 2, 3)
    ylabel('Carrier Frequency (Hz)')
    xlabel('Modulation Frequency (Hz)')
    
    print('-dpng', '-r300', ...
        [exp_dir 'plots/tagged_spect/' strrep(fname, '/', '_')])
end

% Plot the average over subjects
i_cond = 1;
for freq = exp_params.tagged_freqs
    for side = {'left' 'right'}
        subplot(2, 2, i_cond);
        imagesc(mod_freq, car_freq, ...
            squeeze(nanmean(all_spectra.(side{1}).(['f' num2str(freq)]), 1)))
        set(gca,'YDir','normal')
        title(sprintf('%s, %i Hz', side{1}, freq))
        i_cond = i_cond + 1;
    end
end
subplot(2, 2, 3)
ylabel('Carrier Frequency (Hz)')
xlabel('Modulation Frequency (Hz)')

print('-dpng', '-r300', [exp_dir 'plots/tagged_spect/avg'])


%% Plot CFC - RESS

clear variables
rs_setup

freqs = 30:90;
sides = {'left' 'right'};

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    d = load([exp_dir 'cfc/ress/' fname '/cfc']);
    mod_freq_sel = d.mod_freq < 30;
    mf = d.mod_freq(mod_freq_sel);
    
    for i_freq = 1:2
        for i_side = 1:2
            subplot(2, 2, i_side + (2 * (i_freq - 1)))
            if i_side == 2 % Convert from side of space to RESS
                freq_inx = 1 + mod(i_freq, 2);
            elseif i_side == 1
                freq_inx = i_freq;
            end
            x = d.cfc_data.(sides{i_side})(:,mod_freq_sel,freq_inx);
            overall_cfc(i_subject,freq_inx,i_side,:,:) = x;
            imagesc(mf, freqs, x)
            set(gca, 'YDir', 'normal')
            %colorbar;
            title(sprintf('%i Hz, %s', ...
                exp_params.tagged_freqs(i_freq), ...
                sides{i_side}))
        end
    end    
    ylabel('Carrier frequency (Hz)')
    xlabel('Modulation frequency (Hz)')

%     print('-dpng', '-r300', ...
%         [exp_dir 'plots/cfc/' strrep(fname, '/', '_')])
end

% Avg over subjects
tagged = nan(2, 2, length(mf)); % Side * TagFreq * ModFreq
for i_freq = 1:2
    for i_side = 1:2
        subplot(2, 2, i_side + (2 * (i_freq - 1)))
        if i_side == 2 % Convert from side of space to RESS
            freq_inx = 1 + mod(i_freq, 2);
        elseif i_side == 1
            freq_inx = i_freq;
        end
        x = squeeze(mean(overall_cfc(:,freq_inx,i_side,:,:), 1));
        imagesc(mf, freqs, x)
        set(gca, 'YDir', 'normal')
        %colorbar;
        title(sprintf('%i Hz, %s', ...
            exp_params.tagged_freqs(i_freq), ...
            sides{i_side}))
        
        % Take out modulation at the tagged frequency
        tagged_freq_inx = find(freqs == exp_params.tagged_freqs(i_freq));
        tagged(i_side,i_freq,:) = x(tagged_freq_inx,:);
    end
end    
ylabel('Carrier frequency (Hz)')
xlabel('Modulation frequency (Hz)')
% print('-dpng', '-r300', [exp_dir 'plots/cfc/avg'])

% Plot a line-plot of power fluctuations at the tagged frequency
close all
plot(mf, db(squeeze(mean(mean(tagged, 1), 2)), 'power'), 'linewidth', 2)
xlabel('Modulation frequency (Hz)')
ylabel('Power (dB)')
xlim([0 15])
box off

width = 5;
height = 4;
set(gcf,'units','centimeters')
set(gcf,'paperunits','centimeters')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize', [width height])
set(gcf,'paperposition',[0,0,width,height])
set(gcf, 'renderer', 'painters');

print('-depsc', [exp_dir 'plots/cfc/avg_line_tagged'])


%% Plot CFC - computed over raw channels

clear variables
rs_setup

freqs = 30:90;
sides = {'left' 'right'};

% Read in 1 preprocessed datafile to get channel names
d = load('preproc/target/181009_b46d/181009/preproc.mat');
channel_names = d.data.label;
clear d

% Make an occipital ROI based on SNR topographies
snr_roi = [2032 2033 ...
    2042 2043 ...
    2112 2113];
snr_roi = cellfun(@(n) ['MEG' num2str(n)], ...
    num2cell(snr_roi), ...
    'UniformOutput', false);
snr_roi_inx = ismember(channel_names, snr_roi);

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    d = load([exp_dir 'cfc/all_chans/' fname '/cfc']);
    mf = d.mod_freq;
    mf_sel = mf < 50;
    mf = mf(mf_sel);
    
    % Average over occipital ROI
    x = mean(d.cfc_data(:,mf_sel,snr_roi_inx), 3);
    overall_cfc(i_subject,:,:) = x;
    
    imagesc(mf, freqs, x)
    set(gca, 'YDir', 'normal')
    %colorbar;
    ylabel('Carrier frequency (Hz)')
    xlabel('Modulation frequency (Hz)')
    xlim([0 50])

    print('-dpng', '-r300', ...
        [exp_dir 'plots/cfc/all_chans/' strrep(fname, '/', '_')])
end


% Plot average over subjects
imagesc(mf, freqs, squeeze(nanmean(overall_cfc, 1)))
set(gca, 'YDir', 'normal')
%colorbar;
ylabel('Carrier frequency (Hz)')
xlabel('Modulation frequency (Hz)')

print('-dpng', '-r300', [exp_dir 'plots/cfc/all_chans/avg'])

%% Plot cross-correlation of power at the two tagged frequencies

clear variables
rs_setup

close all

xc = nan([height(subject_info), 101]); % Subj x Lag
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    x = load([exp_dir 'tfr/trial/' fname '/xcorr']);
    xc(i_subject,:) = nanmean(x.xc, 1);
end

% Plot avg over subjects
subplot(2,1,1)
plot(x.t_lags, xc, '-', 'color', [0.7 0.7 1])
hold on
plot(x.t_lags, nanmean(xc, 1), '-b', 'LineWidth', 2.5)
plot([-1 1], [0 0], '-k')
plot([0 0], [-1 1], '-k')
hold off
ylim([-0.05 0.3])
xlim([-0.5 0.5])
xlabel('Lag (s)')
ylabel('Correlation')

% % Compute and plot the FFT of the cross-correlation
% % uk.mathworks.com/help/matlab/examples/fft-for-spectral-analysis.html
subplot(2,1,2)
nfft = 2^7;
sample_per = mean(diff(x.t_lags));
Fs = 1 / sample_per;
f = (1/sample_per) * (0:(nfft / 2)) / nfft;
y = fft(xc, nfft, 2);
Pyy = 1 / (nfft * Fs) * abs(y(:,1:nfft/2+1)) .^ 2; % Power spectrum

plot(f, db((Pyy), 'power'), '-', 'color', [0.7 0.7 1]);
hold on
plot(f, db(nanmean(Pyy, 1), 'power'), '-b', 'LineWidth', 2.5)
hold off
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0 12])

print('-dpng', '-r300', [exp_dir 'plots/xcorr'])


%% Compare HF power in hits vs misses
clear variables
close all
rs_setup

% Which TFR window to use?
win_size = 0.2;
win_str = sprintf('win_%.1fs', win_size);
tfr_dir = [exp_dir 'tfr/' win_str '/'];

% Array for all data: Subj x Side x TaggedFreq x Hit X TFRfreq x Time
agg_data = nan([height(subject_info), 2, 2, 2, 46, 101]);

for i_subject = 1:height(subject_info)
        if subject_info.exclude(i_subject)
        continue
    end

    % Read in the data segmented around targets
    behav = rs_behavior(i_subject);
    fname = subject_info.meg{i_subject};
    fn = [tfr_dir 'target/' fname '/high'];
    d = load(fn);
    d = d.high_freq_data;

    % Information about each trial
    trial_numbers = d.trialinfo(:,2);
    targ_side_per_trial = behav.target_side(trial_numbers);
    targ_freq_per_trial = behav.target_side_freq(trial_numbers);

    % Put the data together
    targ_side_labels = {'left' 'right'};
    for i_side = 1:2
        targ_side = targ_side_labels{i_side};
        sel_targ_side = strcmp(targ_side_per_trial, targ_side);
        for i_freq = 1:2
            targ_freq = exp_params.tagged_freqs(i_freq);
            sel_targ_freq = targ_freq_per_trial == targ_freq;
            for hit = [0 1]
                sel_hit = d.trialinfo(:,1) == hit;
                trial_sel = sel_targ_side & sel_targ_freq & sel_hit;
                cfg = [];
                cfg.trials = trial_sel;
                cfg.channel = targ_side;
                d_sub = ft_selectdata(cfg, d);
                x = squeeze(nanmean(d_sub.powspctrm, 1)); % Avg over trials
                agg_data(i_subject,i_side,i_freq,hit+1,:,:) = x;
            end
        end
    end
end
agg_data_orig = agg_data;

% Normalize power within each participant, side, freq
% Divide by the area under the curve
for i_side = 1:2
    for i_freq = 1:2
        for i_subject = 1:size(agg_data, 1)
            x = agg_data_orig(i_subject,i_side,i_freq,:,:,:); % Select data
            x_bar = mean(x, 4); % Avg over hits and misses
            auc = sum(sum(x_bar, 5), 6); % Area under the curve 
            agg_data(i_subject,i_side,i_freq,:,:,:) = x / auc;
        end
    end
end 

% plot it

targ_side_labels = {'left' 'right'};
hit_labels = {'miss' 'hit'};

i_plot = 1;
for i_side = 1:2
    targ_side = targ_side_labels{i_side};
    for i_freq = 1:2
        targ_freq = exp_params.tagged_freqs(i_freq);
        

        for hit = [0 1]
            x = agg_data(:,i_side,i_freq,hit+1,:,:);
            x = nanmean(x, 1); % Average over subjects
            x = squeeze(x);
            hit_miss(:,:,hit+1) = x;
            
            subplot(4, 3, i_plot)
            imagesc(d_sub.time, d_sub.freq, x)
            set(gca, 'YDir', 'normal')
            title(sprintf('%s, %i Hz, %s', ...
                targ_side, targ_freq, hit_labels{hit+1}))
            
            i_plot = i_plot + 1;
        end
        
        subplot(4, 3, i_plot)
        xdiff = hit_miss(:,:,2) - hit_miss(:,:,1);
        clims = [-1 1] * max(max(abs(xdiff)));
        imagesc(d_sub.time, d_sub.freq, xdiff, clims)
        set(gca, 'YDir', 'normal')
        title(sprintf('%s, %i Hz, Diff', ...
            targ_side, targ_freq))

        i_plot = i_plot + 1;
    end
end
print('-dpng', '-r300', ...
    [exp_dir 'plots/accuracy/high_freq/hf_acc_by_side-' win_str '.png'])


%% Plot it, collapsing over side
close all
% figure
i_plot = 1;
for i_freq = 1:2
    targ_freq = exp_params.tagged_freqs(i_freq);

    for hit = [0 1]
        x = agg_data(:,:,i_freq,hit+1,:,:);
        x = nanmean(x, 2); % Average over sides
        x = nanmean(x, 1); % Average over subjects
        x = squeeze(x);
        hit_miss(:,:,hit+1) = x;

        subplot(2, 3, i_plot)
        imagesc(d_sub.time, d_sub.freq, x)
        set(gca, 'YDir', 'normal')
%         title(sprintf('%i Hz, %s', ...
%             targ_freq, hit_labels{hit+1}))
        title(hit_labels{hit+1})

        i_plot = i_plot + 1;
    
    end

    subplot(2, 3, i_plot)
    xdiff = hit_miss(:,:,2) - hit_miss(:,:,1);
    clims = [-1 1] * max(max(abs(xdiff)));
    imagesc(d_sub.time, d_sub.freq, xdiff, clims)
    set(gca, 'YDir', 'normal')
%     title(sprintf('%i Hz, Diff', ...
%         targ_freq))
    title('diff')

    i_plot = i_plot + 1;
    
end

subplot(2,3,4)
ylabel('Frequency (Hz)')
xlabel('Time (s)')

width = 12;
height = 6;
% set(gcf,'units','centimeters')
set(gcf,'paperunits','centimeters')
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf,'papersize', [width height])
set(gcf,'paperposition',[0,0,width,height])
% set(gcf, 'renderer', 'painters');

print('-dpng', '-r300', ...
    [exp_dir 'plots/accuracy/high_freq/hf_acc-' win_str '.png'])



%% Plot the difference as the ratio of power in hits to misses

i_plot = 1;
for i_freq = 1:2
    targ_freq = exp_params.tagged_freqs(i_freq);

    for hit = [0 1]
        x = agg_data_orig(:,:,i_freq,hit+1,:,:);
        x = nanmean(x, 2); % Average over sides
        x = nanmean(x, 1); % Average over subjects
        x = squeeze(x);
        hit_miss(:,:,hit+1) = x;

        subplot(2, 3, i_plot)
        imagesc(d_sub.time, d_sub.freq, x)
        set(gca, 'YDir', 'normal')
        title(hit_labels{hit+1})

        i_plot = i_plot + 1;
    
    end

    subplot(2, 3, i_plot)
    h = hit_miss(:,:,2);
    m = hit_miss(:,:,1);
    %xcomp = h ./ m;
    xcomp = (h-m) ./ (h+m);
    %clims = [-1 1] * max(max(abs(xdiff)));
    imagesc(d_sub.time, d_sub.freq, xcomp) %, clims)
    set(gca, 'YDir', 'normal')
    title('diff')

    i_plot = i_plot + 1;
    
end

subplot(2,3,4)
ylabel('Frequency (Hz)')
xlabel('Time (s)')

width = 12;
height = 6;
set(gcf,'paperunits','centimeters')
set(gcf,'paperposition',[0,0,width,height])

print('-dpng', '-r300', ...
    [exp_dir 'plots/accuracy/high_freq/hf_acc-' win_str '_ratio.png'])

% Run some very simple statistics
% Get a map of t-statistics 
% agg_data: Subj x Side x TaggedFreq x Hit X TFRfreq x Time
x = nanmean(agg_data_orig, 2); % Average over sides
h = x(:,:,:,2,:,:);
m = x(:,:,:,1,:,:);
%xcomp = h ./ m;
xcomp = (h-m) ./ (h+m);
xcomp = squeeze(xcomp);

[h,p,ci,stats] = ttest(xcomp, 0, 'dim', 1, 'alpha', 0.01);

for i_freq = 1:2

    f = exp_params.tagged_freqs(i_freq);

    % Plot t-values
    subplot(2,2,(i_freq-1)*2 + 1)
    t = squeeze(stats.tstat(:,i_freq,:,:));
    %clims = [-1 1] * max(max(abs(t)));
    clims = [-1 1] * 2.5;
    imagesc(d_sub.time, d_sub.freq, t, clims)
    set(gca, 'YDir', 'normal')
    title(sprintf('%i Hz, t', f))
    %colorbar;
    hold on
    plot(0.5 * [-1 1], [f f], '--w')
    hold off

    % Plot points where this sample is signif
    subplot(2,2,(i_freq-1)*2 + 2)
    imagesc(d_sub.time, d_sub.freq, squeeze(h(:,i_freq,:,:)))
    set(gca, 'YDir', 'normal')
    title('Signif')
    hold on
    plot(0.5 * [-1 1], [f f], '--w')
    hold off
end

subplot(2,2,3)
ylabel('Frequency (Hz)')
xlabel('Time (s)')

width = 8;
height = 6;
set(gcf,'paperunits','centimeters')
set(gcf,'paperposition',[0,0,width,height])

print('-dpng', '-r300', ...
    [exp_dir 'plots/accuracy/high_freq/hf_acc-' win_str '_ratio_stats.png'])


%% Compare power at the target/non-target stimulus between hits and misses

clear variables
close all
rs_setup

win_size = 0.1; % Size of the TFR window used

% Read in the data
powdiff = nan(height(subject_info), 2, 101); % Subject x Accuracy x Time
powdiff_all = cell([1 height(subject_info)]);
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    p = rs_powerdiff(i_subject, win_size, 'target');
    powdiff_all{i_subject} = p;
    for hit = 0:1
        x = nanmean(p.powdiff(p.trialinfo(:, 1) == hit, :), 1);
        powdiff(i_subject, hit+1, :) = x;
    end
end
    
% Plot it

% Individual subjects
for i_subject = 1:height(subject_info)
    
    % Misses
    x = squeeze(powdiff(i_subject,1,:));
    plot(p.time, x, '-b', 'LineWidth', 2)
    hold on
    
    % Hits
    x = squeeze(powdiff(i_subject,2,:));
    plot(p.time, x, '-r', 'LineWidth', 2)

    % Labels
    xlabel('Time (s)')
    ylabel(sprintf('Targ - Non-targ power\n(Z)'))
    text(-0.4, 0.15, 'Hit', 'color', 'r')
    text(-0.4, 0.11, 'Miss', 'color', 'b')
    hold off
    
    fn = sprintf('targ-non-diff_win%.1fs_%s', ...
        win_size, strrep(subject_info.meg{i_subject}, '/', '_'));
    print('-dpng', '-r300', [exp_dir 'plots/accuracy/high_freq/' fn '.png'])
end
close all

%% All subjects

sterr = @(x) nanstd(x, 1) / sqrt(size(x, 1));
subplot(2, 1, 1)

% Misses
x = squeeze(powdiff(:,1,:));
x_mean = nanmean(x, 1);
x_se = sterr(x);
fill([p.time fliplr(p.time)], [x_mean + x_se, fliplr(x_mean - x_se)], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on
plot(p.time, x_mean, '-b', 'LineWidth', 2)

% Hits
x = squeeze(powdiff(:,2,:));
x_mean = nanmean(x, 1);
x_se = sterr(x);
fill([p.time fliplr(p.time)], [x_mean + x_se, fliplr(x_mean - x_se)], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
plot(p.time, x_mean, '-r', 'LineWidth', 2)
plot(p.time(logical(h)), zeros(sum(h)) - 0.09, 'k*')

% Labels
xlabel('Time (s)')
ylabel(sprintf('Targ - Non-targ power\n(Z)'))
text(-0.4, 0.15, 'Hit', 'color', 'r')
text(-0.4, 0.11, 'Miss', 'color', 'b')

plot([min(p.time) max(p.time)], [0 0], '-k')
hold off

% Difference: Hit - Miss
subplot(2, 1, 2)
[h,pval,ci] = ttest2(squeeze(powdiff(:,2,:)), squeeze(powdiff(:,1,:)));
x = squeeze(powdiff(:,2,:) - powdiff(:,1,:));
x_mean = nanmean(x, 1);
x_se = sterr(x) * 1.96; % 95% CI
fill([p.time fliplr(p.time)], ... %[x_mean + x_se, fliplr(x_mean - x_se)], ...
    [ci(2,:) fliplr(ci(1,:))], ...
    [0.5 0 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on
plot(p.time, x_mean, '-', 'LineWidth', 2, 'color', [0.5 0 0.5])
plot(p.time(logical(h)), zeros(sum(h)) - 0.19, 'k*')

text(-0.4, 0.3, 'Difference (Hit - Miss)', 'color', [0.5 0 0.5])
xlabel('Time (s)')
ylabel(sprintf('Targ - Non-targ power\n(Z)'))
plot([min(p.time) max(p.time)], [0 0], '-k')
hold off

% % Plot significance for each sample
% hold on
% plot(p.time(h == 1), ...
%     ones([1 sum(h)]) * min(min(ci)), ...
%     'ok', 'LineWidth', 2)
% hold off

fn = sprintf('targ-non-diff_win%.1fs', win_size);
print('-dpng', '-r300', [exp_dir 'plots/accuracy/high_freq/' fn '.png'])


% How many significant points did we get?
sum(h)
% How many would we expect to get by chance?
length(h) * 0.05

%% Spectrum of avg target/non-target power-diffs

% Test whether there are stropnger rhythms during hits than during misses.
% Does this analysis actually make sense? I think what we'd actually expect
% is not that the strength of oscillations differ between hits and misses,
% but that the *phase* differs.

n_timepoints = floor(size(powdiff, 3) / 2);
nfft = 2 ^ ceil(log2(n_timepoints));
sample_per = mean(diff(p.time));
Fs = 1 / sample_per;
f = (1/sample_per) * (0:(nfft / 2)) / nfft;
y = fft(powdiff(:,:,1:n_timepoints), nfft, 3);
y = y(:,:,1:nfft/2+1);
Pyy = 1 / (nfft * Fs) * abs(y(:,:,1:nfft/2+1)) .^ 2; % Power spectrum

% Half-assed stats
ttest2(squeeze(Pyy(:,1,:)), squeeze(Pyy(:,2,:)))

% Plot the power spectrum averaged over participants
line_colors = {'b' 'r'};
tag = {'Miss' 'Hit'};
subplot(1, 2, 1)
for i_hit = 1:2
    plot(f, log10(squeeze(nanmean(Pyy(:, i_hit, :), 1))), ...
        'color', line_colors{i_hit})
    hold on
end
hold off
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Power (log10)')
text(15, -3.4, 'Hit', 'color', line_colors{2})
text(15, -3.5, 'Miss', 'color', line_colors{1})
title('Power of Target-Nontarget env.')

% Difference between the spectra by subject
subplot(1, 2, 2)
Pyy_diff = squeeze(log10(Pyy(:,2,:)) - log10(Pyy(:,1,:)));
plot(f, Pyy_diff, ...
    'color', [0.8 0.5 0.8])
hold on
plot(f, nanmean(Pyy_diff, 1), ...
    'linewidth', 2, ...
    'color', [0.5 0 0.5])
plot([0 50], [0 0], '--k')
hold off
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Power (log10)')
title('Diff in Targ-Nontarg spectra (Hit-Miss)')


fn = sprintf('targ-non-diff_win%.1fs_spectra', win_size);
print('-dpng', '-r300', [exp_dir 'plots/accuracy/high_freq/' fn '.png'])


%% Look at differences in phase

% I think we want to read in the trialwise data first
% For each trial, do an FFT
% Get diff in complex FFT coefficient for each trial
% Average across trials (within each frequency)
% Compute the angle and radius of the averaged coefficient

% Compare phase using the Hilbert transform
h = nan(size(powdiff));
for i_subject = 1:height(subject_info)
    for i_hit = 1:2
        x = squeeze(powdiff(i_subject, i_hit, :));
        h(i_subject, i_hit, :) = hilbert(x);
    end
end
phi = angle(h);
phi_diff = squeeze(phi(:,2,:) - phi(:,1,:));
phase_cons = mean(exp(1j * phi_diff), 2);

% Plot it
plot(real(phase_cons), imag(phase_cons), 'bo')
hold on
a = nanmean(phase_cons); % Average phase difference
plot([0 real(a)], [0 imag(a)], '-r', 'linewidth', 1.5)
plot([0 0], [-20 20], '-k')
plot([-20 20], [0 0], '-k')
xlim([-1 1] * 0.3)
ylim([-1 1] * 0.3)
xlabel('Re')
ylabel('Im')
hold off

fn = sprintf('targ-non-diff_win%.1fs_phase', win_size);
print('-dpng', '-r300', [exp_dir 'plots/accuracy/high_freq/' fn '.png'])


%% Hit rate as a function of LF phase (analysis like Fiebelkorn et al 2018)

clear variables
close all
rs_setup

% Read in the data
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    fn = [exp_dir 'tfr/win_0.1s/target/' fname '/lfphase_acc_fieb'];
    d = load(fn);
    grad = load([exp_dir 'grad/' fname '/grad'], 'grad'); % To combine grads
    d.grad = grad.grad;
    d_all(i_subject) = d;
end

% Plot it

% subject * chan * freq
z_planar = nan([height(subject_info) 204 size(d.z, 3)]);

for i_subject = 1:length(d_all)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    
    d = d_all(i_subject);

    % dependence of hit-rate on phase across channels
    close all
    z = squeeze(mean(d.z, 1));
    imagesc(d.freq, 1:length(d.label), z)
    set(gca, 'YDir', 'normal')
    xlabel('Frequency (Hz)')
    ylabel('Channel')
    colorbar;
    print('-dpng', '-r300', ...
        [exp_dir 'plots/accuracy/low_freq/hr_amp_' strrep(fname, '/', '_')])
    
    % topo of rhythmicity at each frequency
    d_topo = []; % Make the data structure
    d_topo.label = d.label;
    d_topo.freq = d.freq;
    d_topo.grad = d.grad;
    d_topo.time = 0;
    d_topo.powspctrm = z;
    d_topo.dimord = 'chan_freq_time';
    cfg = [];
    cfg.method = 'sum';
    d_topo = ft_combineplanar(cfg, d_topo);
    z_planar(i_subject,:,:) = d_topo.powspctrm;
    close all
    figure('position', [200, 200, 2000, 400])
    for i_freq = 1:length(d.freq)
        subplot(1, length(d.freq), i_freq)
        cfg = [];
        cfg.layout = chan.grad_cmb.layout;
        cfg.ylim = d.freq(i_freq) + [-0.1 0.1];
        cfg.zlim = [0 max(max(d_topo.powspctrm))];
        cfg.colorbar = 'no';
        cfg.style = 'straight';
        cfg.comment = 'no';
        cfg.shading = 'interp';
        cfg.markersymbol = '.';
        cfg.gridscale = 200;
        ft_topoplotTFR(cfg, d_topo)
        title(sprintf('%i Hz', d.freq(i_freq)))
    end
    print('-dpng', '-r300', ...
        [exp_dir 'plots/accuracy/low_freq/hr_amp_topo_' strrep(fname, '/', '_')])
end

% Average over subjects

close all
imagesc(d.freq, 1:length(d.label), mean(cat(3, d_all(:).z), 3))
set(gca, 'YDir', 'normal')
xlabel('Frequency (Hz)')
ylabel('Channel')
colorbar;
print('-dpng', '-r300', ...
    [exp_dir 'plots/accuracy/low_freq/hr_amp_avg'])

d_topo.powspctrm = squeeze(nanmean(z_planar, 1));
close all
figure('position', [200, 200, 2000, 400])
for i_freq = 1:length(d.freq)
    subplot(1, length(d.freq), i_freq)
    % Plot it
    cfg = [];
    cfg.layout = chan.grad_cmb.layout;
    cfg.ylim = d.freq(i_freq) + [-0.1 0.1];
    cfg.zlim = [0 max(max(d_topo.powspctrm))];
    cfg.colorbar = 'no';
    cfg.style = 'straight';
    cfg.comment = 'no';
    cfg.shading = 'interp';
    cfg.markersymbol = '.';
    cfg.gridscale = 200;
    ft_topoplotTFR(cfg, d_topo)
    title(sprintf('%i Hz', d.freq(i_freq)))
end
print('-dpng', '-r300', ...
    [exp_dir 'plots/accuracy/low_freq/hr_amp_topo_avg'])


%% Example plot of hit rate by phase for one subject & freq

i_subject = 1;
i_freq = 3;

fname = subject_info.meg{i_subject};
fn = [exp_dir 'tfr/target/' fname '/lfphase_acc_fieb'];
d = load(fn);

% Plot the avg hit rate for one freq in each phase bin
imagesc(1:d.n_bins, 1:length(d.label), d.hit_rate(:,:,i_freq)')
ylabel('Channel')
xlabel('Phase bin')
c = colorbar;
c.Label.String = 'Accuracy';
title(sprintf('Subject %i, %i Hz', i_subject, d.freq(i_freq)))

print('-dpng', '-r300', ...
    [exp_dir 'plots/accuracy/low_freq/example_hr_by_phase'])


%% PBI

clear variables
close all
rs_setup

% Read in the data
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    fn = [exp_dir 'tfr/target/' fname '/lfphase_acc_pbi'];
    d = load(fn);
    grad = load([exp_dir 'grad/' fname '/grad'], 'grad'); % To combine grads
    d.grad = grad.grad;
    d_all(i_subject) = d;
    
    close all
    for i_freq = 1:length(d.freq)
        subplot(3, 4, i_freq)
        x = squeeze(d.pbi(:,i_freq,:));
        clims = [-1 1] * max(max(abs(x)));
        imagesc(d.time, 1:length(d.label), x, clims)
        xlabel('Time (s)')
        ylabel('Channel')
        title(sprintf('%i Hz', d.freq(i_freq)))
    end
    print('-dpng', '-r300', ...
        [exp_dir 'plots/accuracy/low_freq/pbi_' strrep(fname, '/', '_')])
end

% Average over subjects
pbi_all = cat(4, d_all(:).pbi);
close all
for i_freq = 1:length(d.freq)
    subplot(3, 4, i_freq)
    x = squeeze(mean(pbi_all(:,i_freq,:,:), 4));
    clims = [-1 1] * max(max(abs(x)));
    imagesc(d.time, 1:length(d.label), x, clims)
    xlabel('Time (s)')
    ylabel('Channel')
    title(sprintf('%i Hz', d.freq(i_freq)))
end
print('-dpng', '-r300', [exp_dir 'plots/accuracy/low_freq/pbi_avg'])
    
%     % topo of rhythmicity at each frequency
%     d_topo = []; % Make the data structure
%     d_topo.label = d.label;
%     d_topo.freq = d.freq;
%     d_topo.grad = d.grad;
%     d_topo.dimord = d.dimord;
%     d_topo.time = 0;
%     d_topo.powspctrm = d.pbi;
%     cfg = [];
%     cfg.method = 'sum';
%     d_topo = ft_combineplanar(cfg, d_topo);
%     pbi_planar(i_subject,:,:,:) = d_topo.powspctrm;
%     close all
%     figure('position', [200, 200, 2000, 400])
%     for i_freq = 1:length(d.freq)
%         subplot(1, length(d.freq), i_freq)
%         cfg = [];
%         cfg.layout = chan.grad_cmb.layout;
%         cfg.ylim = d.freq(i_freq) + [-0.1 0.1];
%         cfg.zlim = [0 max(max(d_topo.powspctrm))];
%         cfg.colorbar = 'no';
%         cfg.style = 'straight';
%         cfg.comment = 'no';
%         cfg.shading = 'interp';
%         cfg.markersymbol = '.';
%         cfg.gridscale = 200;
%         ft_topoplotTFR(cfg, d_topo)
%         title(sprintf('%i Hz', d.freq(i_freq)))
%     end
%     print('-dpng', '-r300', ...
%         [exp_dir 'plots/accuracy/low_freq/pbi_topo_' strrep(fname, '/', '_')])
% end


%% Alpha peaks

clear all_tfr
clear all_alpha

for i_subject = 1:size(subject_info, 1)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    
    fn = [exp_dir 'alpha_peaks/' strrep(fname, '/', '_')];
    load(fn)
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [7 14];
    cfg.bpfilttype = 'but'; %FIRLS warning: not recommended for neural signals
    data_seg_alpha = ft_preprocessing(cfg, data_seg);

    close all
    for i_chan = 1:length(data_seg.label)
        for i_tagfreq = 1:2 % Which tagged freq is on the left side
            % Avg alpha activity over selected channels & trials
            cfg = [];
            cfg.trials = avg_sels{i_chan,i_tagfreq}.trial;
            cfg.channel = data_seg.label(avg_sels{i_chan,i_tagfreq}.channel);
            d_alpha = ft_timelockanalysis(cfg, data_seg_alpha);
            all_alpha(i_subject,i_chan,i_tagfreq,:) = d_alpha.avg;
            % Compute TFRs
            cfg.method = 'mtmconvol';
            cfg.foi = 55:2:90;
            cfg.taper = 'hanning';
            cfg.t_ftimwin = 7 ./ cfg.foi;
            cfg.toi = 'all';
            d_tfr = ft_freqanalysis(cfg, data_seg);
            x = squeeze(d_tfr.powspctrm);
            all_tfr(i_subject,i_chan,i_tagfreq,:,:) = x;

            i_cond = ((i_chan - 1) * 2) + i_tagfreq;
            x_lim = 0.3;

            % Plot TFR
            subplot(2, 4, i_cond)
            imagesc(d_tfr.time, d_tfr.freq, squeeze(d_tfr.powspctrm))
            set(gca, 'YDir', 'normal')
            xlim([-1 1] * x_lim)
            title(sprintf('%s, %i Hz, n=%i', ...
                data_preproc.label{i_chan}, ...
                exp_params.tagged_freqs(i_tagfreq), ...
                avg_counts(i_chan, i_tagfreq)))

            % Plot alpha
            subplot(2, 4, i_cond + 4)
            plot(t, d_alpha.avg)
            hold on
            plot(0, 0, 'k+')
            hold off
            xlim([-1 1] * x_lim)

        end
    
    % Adjust color scale within each channel
    x = all_tfr(i_subject,i_chan,:,:,:);
    cmax = max(reshape(x, [1 numel(x)]));
    for i_plot = (1:2) + ((i_chan - 1) * 2)
        subplot(2, 4, i_plot)
        caxis([0 cmax])
    end
        
    end

    subplot(2, 4, 1)
    ylabel('Frequency (Hz)')

    subplot(2, 4, 5)
    ylabel('Amplitude (T)')
    xlabel('Time (S)')

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

end


%% Get the avg over subjects

close all
for i_chan = 1:length(data_seg.label)
    for i_tagfreq = 1:2 % Which tagged freq is on the left side
        i_cond = ((i_chan - 1) * 2) + i_tagfreq;
        x_lim = 0.3;

        % Plot TFR
        subplot(2, 4, i_cond)
        x = squeeze(nanmean(all_tfr(:,i_chan,i_tagfreq,:,:), 1));
        imagesc(d_tfr.time, d_tfr.freq, x)
        set(gca, 'YDir', 'normal')
        xlim([-1 1] * x_lim)
        title(sprintf('%s, %i Hz', ...
            data_preproc.label{i_chan}, ...
            exp_params.tagged_freqs(i_tagfreq)))

        % Plot alpha
        subplot(2, 4, i_cond + 4)
        x = squeeze(nanmean(all_alpha(:,i_chan,i_tagfreq,:), 1));
        plot(t, x)
        hold on
        plot(0, 0, 'k+')
        hold off
        xlim([-1 1] * x_lim)

    end
        
    % Adjust color scale
    x = nanmean(all_tfr, 1);
    cmax = max(reshape(x, [1 numel(x)]));
    for i_plot = 1:4
        subplot(2, 4, i_plot)
        caxis([0 cmax])
    end

end

subplot(2, 4, 1)
ylabel('Frequency (Hz)')

subplot(2, 4, 5)
ylabel('Amplitude (T)')
xlabel('Time (S)')

width = 25;
height = 10;
set(gcf,'units','centimeters')
set(gcf,'paperunits','centimeters')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize', [width height])
set(gcf,'paperposition',[0,0,width,height])
set(gcf, 'renderer', 'painters');

fname = subject_info.meg{i_subject};
fn = [exp_dir 'plots/alpha_peaks/avg'];
print('-dpng', fn)














%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% OLD SCRIPTS PAST HERE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
