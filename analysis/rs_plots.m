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
close all

figure('position', [50, 50, 300, 200])
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
    p = p / sum(p); % Normalize
    hold on
    plot(spec.freq, db(p, 'power'))

end
xlabel('Frequency (Hz)')
ylabel('Power')
xlim([0 100])
plot(exp_params.tagged_freqs, [-38 -38], '^r')
hold off
print('-dpng', [exp_dir 'plots/spectra'])


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

approx_eq = @(x,y) abs(x - y) < 0.1;

% Hold onto the data for all subject
% Subject * Freq * Time
overall_data = nan(height(subject_info), 46, 251);

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
plot([-0.5 4.5], [1 1], '--k')
xlim([-0.5 4.5])
xticks(0:4)
hold off

print('-dpng', '-r300', ...
    [exp_dir 'plots/stim_onset/overall'])


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


%% Plot CFC

clear variables
rs_setup

freqs = 30:90;
sides = {'left' 'right'};

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    d = load([exp_dir 'cfc/' fname '/cfc']);
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

    print('-dpng', '-r300', ...
        [exp_dir 'plots/cfc/' strrep(fname, '/', '_')])
end

% Avg over subjects
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
    end
end    
ylabel('Carrier frequency (Hz)')
xlabel('Modulation frequency (Hz)')
print('-dpng', '-r300', [exp_dir 'plots/cfc/avg'])



%% Plot cross-correlation of power at the two tagged frequencies

clear variables
rs_setup

close all

xc = nan([height(subject_info), 51]); % Subj x Lag
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

plot(f, Pyy, '-', 'color', [0.7 0.7 1]);
hold on
plot(f, nanmean(Pyy, 1), '-b', 'LineWidth', 2.5)
hold off
xlabel('Frequency (Hz)')
ylabel('Power')
xlim([0 12])

print('-dpng', '-r300', [exp_dir 'plots/xcorr'])


%% Compare HF power in hits vs misses
clear variables
close all
rs_setup

% Array for all data: Subj x Side x TaggedFreq x Hit X TFRfreq x Time
agg_data = nan([height(subject_info), 2, 2, 2, 46, 51]);

for i_subject = 1:height(subject_info)
        if subject_info.exclude(i_subject)
        continue
    end

    % Read in the data segmented around targets
    behav = rs_behavior(i_subject);
    fname = subject_info.meg{i_subject};
    fn = [exp_dir 'tfr/target/' fname '/high'];
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
                x = squeeze(nanmean(d_sub.powspctrm, 1));
                agg_data(i_subject,i_side,i_freq,hit+1,:,:) = x;
            end
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
print('-dpng', '-r300', [exp_dir 'plots/hf_acc_by_side'])


% Plot it, collapsing over side
figure
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
        title(sprintf('%i Hz, %s', ...
            targ_freq, hit_labels{hit+1}))

        i_plot = i_plot + 1;
    end

    subplot(2, 3, i_plot)
    xdiff = hit_miss(:,:,2) - hit_miss(:,:,1);
    clims = [-1 1] * max(max(abs(xdiff)));
    imagesc(d_sub.time, d_sub.freq, xdiff, clims)
    set(gca, 'YDir', 'normal')
    title(sprintf('%i Hz, Diff', ...
        targ_freq))

    i_plot = i_plot + 1;
end

subplot(2,3,4)
ylabel('Frequency (Hz)')
xlabel('Time (s)')

print('-dpng', '-r300', [exp_dir 'plots/hf_acc'])


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
    fn = [exp_dir 'tfr/target/' fname '/lfphase_acc_fieb'];
    d = load(fn);
    grad = load([exp_dir 'grad/' fname '/grad'], 'grad'); % To combine grads
    d.grad = grad.grad;
    d_all(i_subject) = d;
end

% Plot it

z_planar = nan([height(subject_info) 204 size(d.z, 2)]);

for i_subject = 1:length(d_all)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    
    d = d_all(i_subject);

    % dependence of hit-rate on phase across channels
    close all
    imagesc(d.freq, 1:length(d.label), d.z)
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
    d_topo.powspctrm = d.z;
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
    
    for i_freq = 1:length(d.freq)
        subplot(3, 4, i_freq)
        x = squeeze(d.pbi(:,i_freq,:));
        clims = [-1 1] * max(max(abs(x)));
        imagesc(d.time, d.freq, x, clims)
        xlabel('Time (s)')
        ylabel('Channel')
        title(sprintf('%i Hz', d.freq(i_freq)))
    end
    print('-dpng', '-r300', ...
        [exp_dir 'plots/accuracy/low_freq/pbi_' strrep(fname, '/', '_')])

    
    % topo of rhythmicity at each frequency
    d_topo = []; % Make the data structure
    d_topo.label = d.label;
    d_topo.freq = d.freq;
    d_topo.grad = d.grad;
    d_topo.dimord = d.dimord;
    d_topo.time = 0;
    d_topo.powspctrm = d.pbi;
    cfg = [];
    cfg.method = 'sum';
    d_topo = ft_combineplanar(cfg, d_topo);
    pbi_planar(i_subject,:,:) = d_topo.powspctrm;
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
        [exp_dir 'plots/accuracy/low_freq/pbi_topo_' strrep(fname, '/', '_')])
end



















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
