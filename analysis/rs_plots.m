%% Plot the target threshold over the course of the experiment

rs_setup
close all

clrs = [0 0.6 0; 0.6 0 0.6];

thrsh = nan(height(subject_info), 2, 2); % Subj * Side * RFT_Freq

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    % Load the behavioral data
    fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
    behav = rs_behavior(i_subject);
    subplot(4,4,i_subject)
    hold on
    sides = {'left' 'right'};
    for i_side = 1:2
        s = sides{i_side};
        freqs = [63 78];
        for i_freq = 1:2
            f = freqs(i_freq);
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
            x = behav.target_opacity(inx);
            plot(x, '-', 'color', clr)
            ylim([0 0.4])
            % Store the value
            thrsh(i_subject, i_side, i_freq) = x(end);
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

%% Plot the thresholds for each subject
for i_side = 1:2
    for i_freq = 1:2
        if i_freq == 1
            clr = clrs(1,:);
        else
            clr = clrs(2,:);
        end
        if i_side == 1
            clr = clr + (1 - max(clr)); % Left - lighter
        else
            clr = clr * 0.8; % Right - darker
        end
        plot(1:height(subject_info), thrsh(:,i_side,i_freq), ...
            'o', 'color', clr, 'LineWidth', 2)
        hold on

    end
end
hold off
xlabel('Subject')
ylabel('Threshold')
xlim([0 height(subject_info) + 1])
print('-dpng', [exp_dir 'plots/behav/opacity-threshold'])


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


%% RT histograms

% clear variables
% close all
% rs_setup

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    subplot(4, 4, i_subject)
    behav = rs_behavior(i_subject);
    rt = behav.rt - behav.target_t;
    histogram(rt, 0:0.1:2, ...
        'Normalization', 'count', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 1.5)
%         'EdgeAlpha', 0.4)
%     hold on
end
% hold off
subplot(4,4,1)
ylabel('Count')
xlabel('RT (s)')

print('-dpng', [exp_dir 'plots/behav/rt_hist'])


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
            % Take the first column - first RESS component is the only one
            % that's meaningful.
            maps = subj_map.(side{1}).(['f' num2str(freq)]).maps(:,1);
            % Make data structure to show maps
            d_maps = [];
            d_maps.label = labels;
            d_maps.time = 1;
            d_maps.avg = real(maps);
            d_maps.dimord = 'chan_time';
            d_maps.grad = grad.grad;
            
            % Plot the raw filter coefficients
            subplot(2,2,i_plot)
            mags = endsWith(d_maps.label, '1');
            inx = 1:length(d_maps.label);
            plot(inx(mags), d_maps.avg(mags), 'or')
            hold on
            plot(inx(~mags), d_maps.avg(~mags), 'ob')
            hold off
            
%             % Combine planar gradiometers
%             cfg = [];
%             cfg.method = 'sum';
%             d_maps_gradcmb = ft_combineplanar(cfg, d_maps);
% 
%             % Plot gradiometers
%             figure(1)
%             subplot(2,2,i_plot)
%             cfg = [];
%             cfg.marker = 'on';
%             cfg.markersymbol = 'o';
%             cfg.markersize = 1;
%             cfg.comment = 'no';
%             cfg.style = 'straight';
%             cfg.layout = chan.grad_cmb.layout;
%             cfg.gridscale = 200;
%             ft_topoplotER(cfg, d_maps_gradcmb)
%             title(sprintf('%s, %i Hz', side{1}, freq))
%             
%             % Plot magnetometers
%             figure(2)
%             subplot(2,2,i_plot)
%             cfg.layout = chan.mag.layout;
%             ft_topoplotER(cfg, d_maps)
%             title(sprintf('%s, %i Hz', side{1}, freq))
            
            i_plot = i_plot + 1;
        end
    end
    
    print('-dpng', '-r300', ...
        [exp_dir 'plots/ress_maps/raw-coefs-' ...
        strrep(fname, '/', '_')])

%     figure(1)
%     print('-dpng', '-r300', ...
%         [exp_dir 'plots/ress_maps/grad-' ...
%         strrep(fname, '/', '_')])
% 
%     figure(2)
%     print('-dpng', '-r300', ...
%         [exp_dir 'plots/ress_maps/mag-' ...
%         strrep(fname, '/', '_')])
    
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


%% Simulated data - theta-band modulation

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


%% Simulated data - alpha phase coding

d = rs_simulate_alpha_peaks();

i_trial = 1;

t = d.time{i_trial};
plot(t, d.trial{i_trial}(1,:), '-r', ...
    t, d.trial{i_trial}(2,:) + 1, '-', ...
    t, d.trial{i_trial}(3,:) + 2, '-b')
xlim([1 2])
xlabel('Time (s)')
ylabel('Amplitude')
yticks((1:3) - 0.5)
yticklabels(d.label)
box('off')

print('-dpng', [exp_dir 'plots/alpha_peaks/simulation']);

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
        colorbar('EastOutside');
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
xlim([1 15])
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

% Only keep gradiometers - magnetometers have weird modulated noise
snr_roi = {'MEG2032' 'MEG2033' 'MEG2042' 'MEG2043' 'MEG2112' 'MEG2113'};
snr_roi_inx = ismember(channel_names, snr_roi);

vers = 'all_chans/nodetrend/';

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    d = load([exp_dir 'cfc/' vers fname '/cfc']);
    mf = d.mod_freq;
    mf_sel = (1 < mf) & (mf < 30);
    mf = mf(mf_sel);
    
    % Average over ROI
    x = mean(d.cfc_data(:,mf_sel,snr_roi_inx), 3);
    overall_cfc(i_subject,:,:) = x;
    
    imagesc(mf, freqs, x)
    set(gca, 'YDir', 'normal')
    colorbar;
    ylabel('Carrier frequency (Hz)')
    xlabel('Modulation frequency (Hz)')

    print('-dpng', '-r300', ...
        [exp_dir 'plots/cfc/' vers strrep(fname, '/', '_')])
end


% Plot average over subjects
imagesc(mf, freqs, squeeze(nanmean(overall_cfc, 1)))
set(gca, 'YDir', 'normal')
colorbar;
ylabel('Carrier frequency (Hz)')
xlabel('Modulation frequency (Hz)')

print('-dpng', '-r300', [exp_dir 'plots/cfc/' vers 'avg'])

%% Plot cross-correlation of power at the two tagged frequencies -- OLD

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

%% Cross-correlations -- NEW

clear variables
rs_setup

close all

clrs = {'b' 'r'};

xc = nan([height(subject_info), 2, 101]); % Subj * Hit side * Lag
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    x = load([exp_dir 'xcorr/' fname '/xc']);
    xc(i_subject,:,:) = x.xc;

    plot(x.t_lags, x.xc)
    hold on
    plot([-0.5 0.5], [0 0], '-k')
    hold off
    xlabel('Lag (s)')
    ylabel('Corr')
    legend('Left hits', 'Right hits')
    
    fn = [exp_dir 'plots/xcorr/' strrep(fname, '/', '_')];
    print('-dpng', '-r300', fn);
end

% Plot average over subjects
n = sum(subject_info.exclude == 0);
for i_hit_side = 1:2
    xcx = xc(:,i_hit_side,:);
    m = squeeze(nanmean(xcx, 1));
    sem = squeeze(nanstd(xcx, 1)) ./ sqrt(n);
    plot(x.t_lags, m, clrs{i_hit_side})
    hold on
    plot([-0.5 0.5], [0 0], '-k')
    fill([x.t_lags fliplr(x.t_lags)], ...
        [m + sem; flipud(m - sem)], ...
        clrs{i_hit_side}, ...
        'edgecolor', 'none')
    alpha(0.2)
end    
hold off
xlabel('Lag (s)')
ylabel('Corr')
legend('Left hits', 'Right hits')

fn = [exp_dir 'plots/xcorr/avg'];
print('-dpng', '-r300', fn);



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
%

win_size = 0.1; % Size of the TFR window used
bootstrap = false;

if bootstrap
    k = 1e4;
    warning('Bootstrapping - permuting data')
else
    k = 0;
end

% Read in the data
powdiff = nan(height(subject_info), 2, 101); % Subject x Accuracy x Time
powdiff_all = cell([1 height(subject_info)]);
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    p = rs_powerdiff(i_subject, ...
        win_size, ...
        'target', ...
        k); % Bootstrap samples
    powdiff_all{i_subject} = p;

    
    powdiff(i_subject, 1, :) = nanmean(p.powdiff_miss, 1);
    powdiff(i_subject, 2, :) = nanmean(p.powdiff_hit, 1);
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
    text(-0.4, 0.11, 'Miss', 'color', 'b')
    text(-0.4, 0.15, 'Hit', 'color', 'r')
    hold off
    
    if ~bootstrap
        fn = sprintf('targ-non-diff_win%.1fs_%s', ...
            win_size, strrep(subject_info.meg{i_subject}, '/', '_'));
        print('-dpng', '-r300', ...
            [exp_dir 'plots/accuracy/high_freq/' fn '.png'])
    end
end
% close all

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
% plot(p.time(logical(h)), zeros(sum(h)) - 0.09, 'k*')

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
plot(p.time(logical(h)), zeros(sum(h)) + min(min(min(powdiff))), 'k*')

text(-0.4, 0.3, 'Difference (Hit - Miss)', 'color', [0.5 0 0.5])
xlabel('Time (s)')
ylabel(sprintf('Targ - Non-targ power\n(Z)'))
plot([min(p.time) max(p.time)], [0 0], '-k')
hold off


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

%% Pairwise phase consistency metric (PCM) following Landau et al 2015

% From Landau et al (2015)
% This phase consistency metric (PCM) assumes a value of -1 if target 
% events leading to hits versus misses are preceded by perfectly opposing 
% LGA phases.  It  assumes  a  value  of  1  if  target  events  leading  
% to hits versus misses are preceded by perfectly aligned LGA phases.

clear variables
close all
rs_setup

% Read in the data
pcm = nan(height(subject_info), 33); % Subject * Freq
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    [f, c] = rs_pcm(i_subject);
    pcm(i_subject,:) = c;   
end

% Simple stats
[h, p, ci] = ttest(pcm);

% Plot it
% subplot(1,2,1)
fill([f, fliplr(f)], [ci(1,:) fliplr(ci(2,:))], ...
    'b', 'edgecolor', 'none', 'facealpha', 0.3)
hold on
% plot(f, pcm, '-', 'LineWidth', 0.5, 'color', [1 1 1] * 0.6)
hold on
plot(f, nanmean(pcm, 1), '-b', 'LineWidth', 2)
plot([min(f) max(f)], [0 0], '-k')
hold off
xlim([1 20])
xlabel('Frequency (Hz)')
ylabel('PCM')

% % Plot correlation of PCM and RFT amplitude
% rri = rs_rft_responsivness();
% subplot(1,2,2)
% plot(pcm(:,f==4.6875), rri, 'ok')
% xlabel('PCM')
% ylabel('RRI')

print('-dpng', '-r300', [exp_dir 'plots/accuracy/high_freq/pcm.png'])



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

clear variables
close all

vers = 'prct100'; % Which version of the analysis to run
%vers = 'teta_test';
%vers = 'beta_test';

rs_setup
rft_freqs = exp_params.tagged_freqs;

all_data = cell([1 height(subject_info)]);

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    else
        fname = subject_info.meg{i_subject};
        fn = [exp_dir 'alpha_peaks/' vers '/' strrep(fname, '/', '_')];
        data = load(fn, 'cond_counts', 'cond_alpha', 'cond_tfr', 'segments');
        all_data{i_subject} = data;
    end

    side_labels = {'left' 'right'};
    x_lim = 0.3;
    %close all
    i_cond = 1;
    for i_chan = 1:2
        for i_rft_freq = 1:2
            % Plot the TFR
            width = 0.19;
            spacing = 0.05;
            lpos = spacing + (width + spacing) * (i_cond-1);
            subplot('position', [lpos, 0.35, width, 0.55])
            d = data.cond_tfr{i_chan, i_rft_freq};
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
                data.cond_counts(i_chan, i_rft_freq)))

            % Plot the alpha power
            subplot('position', [lpos, 0.13, width, 0.17])
            d = data.cond_alpha{i_chan, i_rft_freq};
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

    fig_width = 25;
    fig_height = 10;
    set(gcf,'units','centimeters')
    set(gcf,'paperunits','centimeters')
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'papersize', [fig_width fig_height])
    set(gcf,'paperposition',[0,0,fig_width,fig_height])
    set(gcf, 'renderer', 'painters');

    fn = [exp_dir 'plots/alpha_peaks/' strrep(fname, '/', '_')];
    print('-dpng', fn)
end

% Look for rhythms in HF power, when locked to alpha peaks
% Is this double-dipping?

clear all_spectra
close all

for i_subject = [1:height(subject_info) 0]
    if i_subject == 0
        % Average across subjects
        fname = 'avg';
    elseif subject_info.exclude(i_subject)
        continue
    else
        fname = subject_info.meg{i_subject};
        data = all_data{i_subject};        
    end

    side_labels = {'left' 'right'};
    x_lim = 0.3;
    %close all
    i_cond = 1;
    for i_chan = 1:2
        for i_rft_freq = 1:2
            subplot(2,2,i_cond)
            
            data_cond = data.cond_tfr{i_chan, i_rft_freq};
            t_inx = (-x_lim < data_cond.time) & (data_cond.time < x_lim);
            x = squeeze(data_cond.powspctrm);
            x = x(:,t_inx);
            % Divide by mean in each freq
            for i_freq = 1:length(data_cond.freq)
               x(i_freq,:) = x(i_freq,:) / nanmean(x(i_freq,:));
            end
           
            % FFT of power fluctuations at freq
            nfft = 2 ^ nextpow2(size(x, 2));
            nfft = size(x, 2);
            sample_per = mean(diff(data_cond.time));
            Fs = 1 / sample_per;
            f = (1/sample_per) * (0:(nfft / 2)) / nfft;
            f_inx = 1 < f & f < 31;
            y = fft(x, nfft, 2);
            Pyy = 1 / (nfft * Fs) * abs(y(:,1:round(nfft/2)+1)) .^ 2; % Power spectrum
            Pyy = Pyy(:,f_inx);
            f = f(f_inx);
            if i_subject == 0 %% Averages
                Pyy = squeeze(nanmean(all_spectra(:,i_cond,:,:), 1));
            else
                all_spectra(i_subject,i_cond,:,:) = Pyy;
            end
            imagesc(f, data_cond.freq, db((Pyy), 'power'))
            set(gca, 'YDir', 'normal')
            xlabel('FFT Frequency (Hz)')
            ylabel('TFR Frequency (Hz)')
            xlim([1 30])
            colorbar;
            
            % Add a title
            title(sprintf('%s, %i Hz', ...
                side_labels{i_chan}, ...
                rft_freqs(i_rft_freq)))
            
            i_cond = i_cond + 1;
        end
    end
    
    fn = [exp_dir 'plots/alpha_peaks/spect_' strrep(fname, '/', '_')];
    print('-dpng', fn)

end

%{
% simulated data
x = rand(size(x));
x = x .* sin(data_cond.time(t_inx) * 10 * 2 * pi);
nfft = 2 ^ nextpow2(size(x, 2));
sample_per = mean(diff(data_cond.time));
Fs = 1 / sample_per;
f = (1/sample_per) * (0:(nfft / 2)) / nfft;
y = fft(x, nfft, 2);
Pyy = 1 / (nfft * Fs) * abs(y(:,1:nfft/2+1)) .^ 2; % Power spectrum
imagesc(f, data_cond.freq, db((Pyy), 'power'))
set(gca, 'YDir', 'normal')
xlim([1 30])
%}

%% Alpha-peaks - compare the two tagged frequencies

clear variables
close all

vers = 'overlap/indiv_chans/alpha/';
x_lim = 0.3;

rs_setup
rft_freqs = exp_params.tagged_freqs;

clear all_spectra
clear all_x

for i_subject = [1:height(subject_info) 0] %inf
    if i_subject == 0 
        fname = 'avg';
    elseif i_subject == inf
        fname = 'SIMULATED';
        fn = [exp_dir 'alpha_peaks/' vers '/' strrep(fname, '/', '_')];
        data = load(fn);
    elseif subject_info.exclude(i_subject)
        continue
    else
        fname = subject_info.meg{i_subject};
        fn = [exp_dir 'alpha_peaks/' vers '/' strrep(fname, '/', '_')];
        data = load(fn);
    end
    
    if (i_subject == 1) || (i_subject == inf)
        t = data.cond_tfr{1,1}.time;
        t_inx = (-x_lim <= t) & (t <= x_lim);
        f_tfr = data.cond_tfr{1,1}.freq;
        sample_per = mean(diff(t));
        Fs = 1 / sample_per;
    end
        
    side_labels = {'left' 'right'};
    nchans = size(data.cond_alpha, 1);
    for i_chan = 1:nchans
        % We're comparing the stim at 63 Hz to the stim at 78 Hz.
        % Because these were always on opposite sides, this is equivalent
        % to looking within each frequency and comparing stim (L-R)/(L+R)
        d63 = squeeze(data.cond_tfr{i_chan, 1}.powspctrm);
        d78 = squeeze(data.cond_tfr{i_chan, 2}.powspctrm);
        x = (d63 - d78) ./ (d63 + d78);
        x = x(:,t_inx);
        if i_subject == 0 %% Averages
            x = squeeze(nanmean(all_x(:,i_chan,:,:), 1));
        elseif i_subject == inf
            disp('Simulated data')
        else
            all_x(i_subject,i_chan,:,:) = x;
        end
        
        chan_name = data.cond_alpha{i_chan,1}.label{1};

        % Plot details
        n_rows = 3;
        n_cols = 3;
        space = 0.07;
        subplot_width = (1 / n_rows) - space - 0.01;
        subplot_height = (1 / n_cols) - space - 0.01;
        tfr_height = subplot_height * 0.67;
        trace_height = subplot_height * 0.33;
        i_col = ceil(i_chan / n_cols);
        i_row = mod(i_chan, n_rows);
        if i_row == 0
            i_row = 3;
        end
        xpos = space + (i_row - 1) * (space + subplot_width);
        ypos = space + (i_col - 1) * (space + subplot_height);

         % Plot the trace
        subplot('position', [xpos, ypos, subplot_width, trace_height])
         for i_tfr_freq = 1:2
            plot(t, data.cond_alpha{i_chan, i_tfr_freq}.avg)
            hold on
        end
        hold off
        xlim([-1 1] * x_lim)
        xlabel('Time (s)')
        ylabel('Amp (T)')
        y_max = max(data.cond_alpha{i_chan,i_tfr_freq}.avg);
        scaler = floor(log10(y_max));
        ylims = ylim;

        % Plot the difference-over sum
        subplot('position', [xpos, ypos+trace_height, subplot_width, tfr_height])
        imagesc(t, f_tfr, x)
        hold on
        plot([-1 1], 63 * [1 1], '--w')
        plot([-1 1], 78 * [1 1], '--w')
        hold off
        set(gca, 'YDir', 'normal')
        xlim([-1 1] * x_lim)
        xticks([])
        ylabel('Freq (Hz)')
        h = colorbar('East');
        cb_pos = get(h, 'Position');
        set(h, 'Position',...
            [cb_pos(1)+0.05 cb_pos(2) cb_pos(3)/4 cb_pos(4)]);
%         title(chan_name)     
    end

    fig_width = 20;
    fig_height = 20;
    set(gcf,'units','centimeters')
    set(gcf,'paperunits','centimeters')
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'papersize', [fig_width fig_height])
    set(gcf,'paperposition', [0,0,fig_width,fig_height])
    set(gcf, 'renderer', 'painters');
    fn = sprintf('%splots/alpha_peaks/%s/%s', ...
        exp_dir, ...
        vers, ...
        strrep(fname, '/', '_'));  
    print('-dpng', fn)
    
        %{
        % Plot the difference-over-sum
        figure(1)
        subplot(3,3,i_chan)
        imagesc(t, f_tfr, x)
        hold on
        plot([-1 1], 63 * [1 1], '--w')
        plot([-1 1], 78 * [1 1], '--w')
        hold off
        set(gca, 'YDir', 'normal')
        xlim([-1 1] * x_lim)
        xlabel('Time (s)')
        ylabel('Freq (Hz)')
        colorbar('EastOutside')
        title(chan_name)
        
        % Get the spectrum of the difference-over-sum
        figure(2)
        subplot(3,3,i_chan)
        nfft = size(x, 2);
        f_fft = (1/sample_per) * (0:(nfft / 2)) / nfft;
        f_inx = 1 < f_fft & f_fft < 31;
        y = fft(x, nfft, 2);
        Pyy = 1 / (nfft * Fs) * abs(y(:,1:round(nfft/2)+1)) .^ 2; % Pow
        Pyy = Pyy(:,f_inx);
        f_fft = f_fft(f_inx);
        if i_subject == 0 %% Averages
            Pyy = squeeze(nanmean(all_spectra(:,i_chan,:,:), 1));
        elseif i_subject == inf
            disp('Simulated')
        else
            all_spectra(i_subject,i_chan,:,:) = Pyy;
        end
        imagesc(f_fft, f_tfr, db((Pyy), 'power'))
        hold on
        plot(minmax(f_fft), 63 * [1 1], '--w')
        plot(minmax(f_fft), 78 * [1 1], '--w')
        hold off
        set(gca, 'YDir', 'normal')
        xlabel('FFT Freq (Hz)')
        ylabel('TFR Freq (Hz)')
        xlim([1 30])
        colorbar;
        title(chan_name)
        
        % Plot the average alpha traces
        figure(3)
        subplot(3,3,i_chan)
        for i_tfr_freq = 1:2
            plot(t, data.cond_alpha{i_chan, i_tfr_freq}.avg)
            hold on
        end
        hold off
        xlim([-1 1] * x_lim)
        title(chan_name)
    end

    
    fignames = {'tfr' 'spectra' 'alpha'};
    for i_fig = 1:3
        figure(i_fig)
        fig_width = 30;
        fig_height = 15;
        set(gcf,'units','centimeters')
        set(gcf,'paperunits','centimeters')
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf,'papersize', [fig_width fig_height])
        set(gcf,'paperposition', [0,0,fig_width,fig_height])
        set(gcf, 'renderer', 'painters');
        fn = sprintf('%splots/alpha_peaks/%s/%s-%s', ...
            exp_dir, ...
            vers, ...
            fignames{i_fig}, ...
            strrep(fname, '/', '_'));  
        print('-dpng', fn)
    end
    
    %}

end


%% Plot alpha topography


clear variables
close all
rs_setup

tfr_dir = [exp_dir 'tfr/win_0.2s/'];

% Read in all data
% Array for all data: Subj x Chan x TFRfreq x Time x TargSide x Hit 
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end

    % Read in the data segmented around targets
    behav = rs_behavior(i_subject);
    fname = subject_info.meg{i_subject};
    fn = [tfr_dir 'trial/' fname '/low'];
    d = load(fn);
    d = d.low_freq_data;

    % Convert from complex fourier spec to power spec
    d.powspctrm = abs(d.fourierspctrm) .^ 2;
    d = rmfield(d, 'fourierspctrm');

    % Read in the grad structure (to be able to combine grads)
    fn = [exp_dir 'grad/' fname '/grad'];
    grad = load(fn);
    d.grad = grad.grad;
    % Combine planar gradiometers
    cfg = [];
    cfg.method = 'sum';
    cfg.updatesens = 'yes';
    d = ft_combineplanar(cfg, d);

    % Average in the alpha band
    cfg = [];
    cfg.trials = 'all';
    cfg.avgoverrpt = 'yes';
    
%     cfg.latency = [0.5 5];
%     cfg.avgovertime = 'yes';
%     cfg.nanmean = 'yes';
 
    cfg.frequency = [7 14];
    cfg.avgoverfreq = 'yes';
    cfg.nanmean = 'yes';
 
    cfg.channel = 'MEG*1';
    d_mag = ft_selectdata(cfg, d);
    
    cfg.channel = 'MEG*3';
    d_grad = ft_selectdata(cfg, d);
    
    % Plot them
    subplot(1,2,1)
    cfg = [];
    cfg.xlim = [0.5 4];
    cfg.style = 'straight';
    cfg.comment = 'no';
    cfg.shading = 'interp';
    cfg.markersymbol = '.';
    cfg.gridscale = 200;
    cfg.colorbar = 'South';
    cfg.layout = chan.mag.layout;
    ft_topoplotTFR(cfg, d_mag)
    title('Magnetometers')

    subplot(1,2,2)
    cfg.layout = chan.grad_cmb.layout;
    ft_topoplotTFR(cfg, d_grad)
    title('Gradiometers')
    
    fn = [exp_dir 'plots/alpha_topo/' strrep(fname, '/', '_')];
    print('-dpng', fn)
    
end

%% Identify channels where alpha decreases relative to baseline

clear variables
close all
rs_setup

%
for i_subject = [1:height(subject_info) 0]
    if i_subject == 0
%         figure
        d.tfr.powspctrm = squeeze(nanmean(overall_powspctrm, 1));
        fname = 'avg';
    elseif subject_info.exclude(i_subject)
        continue
    else
        % Read in the data segmented around targets
        fname = subject_info.meg{i_subject};
        fn = [exp_dir 'tfr/baseline_alpha/' fname '/baseline_alpha'];
        d = load(fn);
        x = d.tfr.powspctrm;
        overall_powspctrm(i_subject,:,:,:) = x;
    end

    high_alpha_chans = {'192x', '194x', '191x',...  % '204x', 
                        '234x', '232x', '231x',};   % '203x',
    chan_names = @(n) cellfun(@(c) ['MEG' c(1:3) num2str(n)], ...
        high_alpha_chans, ...
        'UniformOutput', false)
    high_alpha_chans = [chan_names(1) chan_names(2) chan_names(3)];
    cfg = [];
    cfg.channel = high_alpha_chans; %snr_roi;
    cfg.baseline = [-0.5 -0.2];
    cfg.baselinetype = 'relative';
    cfg.xlim = [-0.5 2];
    cfg.title = ' ';
    ft_singleplotTFR(cfg, d.tfr)
    
    fn = [exp_dir 'plots/alpha_baseline/' strrep(fname, '/', '_')];
    print('-dpng', fn)

end

%% Plot the topography of the avg response

grad = load([exp_dir 'grad/' subject_info.meg{1} '/grad']);
d.tfr.grad = grad.grad;

cfg = [];
cfg.method = 'sum';
x = ft_combineplanar(cfg, d.tfr);

cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relative';
x = ft_freqbaseline(cfg, x);

% ROI to highlight
high_alpha_chans = {'192x', '194x', '191x',...  % '204x', 
                    '234x', '232x', '231x',};   % '203x', 
mag_names = cellfun(@(s) ['MEG' s(1:(end-1)) '1'], ...
    high_alpha_chans, ...
    'UniformOutput', false);
cmb_grad_names = cellfun(...
    @(s) ['MEG' s(1:(end-1)) '2+' s(1:(end-1)) '3'], ...
    high_alpha_chans, ...
    'UniformOutput', false);

cfg = [];
cfg.xlim = [0.2 0.8];
cfg.ylim = [7 14];
cfg.layout = chan.grad_cmb.layout;
cfg.style = 'straight';
cfg.highlight = 'on';
cfg.highlightchannel = cmb_grad_names;
cfg.highlightsymbol = '.';x
cfg.highlightsize = 15;
cfg.comment = ' ';
cfg.colorbar = 'SouthOutside';
ft_topoplotTFR(cfg, x)
title('')

fn = [exp_dir 'plots/alpha_baseline/topo'];
print('-dpng', fn)


%% Plot filter characteristics for default bandpass filters

data = [];
data.fsample = 1000;
data.time = {-1:(1/data.fsample):1};
data.trial = {zeros(size(data.time{1}))};
data.trial{1}(1000) = 1;
data.label = {'e1'};

bands = [3 7; 7 14; 15 25];

for i_band = 1:size(bands, 1)
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = bands(i_band,:);
    cfg.bpfilttype = 'but';
    d = ft_preprocessing(cfg, data);
    
    subplot(2,1,1)
    plot(d.time{1}, d.trial{1}, 'LineWidth', 2)
    hold on
    
    % Plot the spectrum of the filter response
    cfg = [];
    cfg.taper = 'hanning';
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    f = ft_freqanalysis(cfg, d);
    
    subplot(2,1,2)
    plot(f.freq, f.powspctrm, 'LineWidth', 2)
    hold on
end

subplot(2,1,1)
hold off
xlim(0.5 * [-1 1])
xlabel('Time (s)')
ylabel('Amplitude')

subplot(2,1,2)
hold off
xlim([0 50])
ylabel('Power')
xlabel('Frequency (Hz)')

print('-dpng', [exp_dir 'plots/filter_char'])

%% Plot the EOG

clear variables
close all
rs_setup

% Array for all data: Subj x Chan x Time x TargSide x Hit 
agg_data_avg = nan([height(subject_info), 2, 2001, 2, 2]);
agg_data_var = agg_data_avg;

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    
%     % Be sure to comment out channel selection, and take EOG chans, and
%     % change the save filepath
%     rs_preproc(i_subject, 'target');

    fname = subject_info.meg{i_subject};
    fn = [exp_dir 'preproc/eog/target/' fname '/preproc'];    
    d = load(fn);
    d = d.data;
    behav = rs_behavior(i_subject);
    
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30;
    d = ft_preprocessing(cfg, d);
        
    % Information about each trial
    trial_numbers = d.trialinfo(:,2);
    targ_side_per_trial = behav.target_side(trial_numbers);

    % Put the data together
    side_labels = {'left' 'right'};
    for i_targ_side = 1:2
        targ_side_sel = strcmp(targ_side_per_trial, side_labels{i_targ_side});
        for hit = 0:1
            hit_sel = d.trialinfo(:,1) == hit;
            trial_sel = hit_sel & targ_side_sel;
            cfg = [];
            cfg.trials = trial_sel;
            cfg.keeptrials = 'no';
            cfg.vartrllength = 2;
            d_sub = ft_timelockanalysis(cfg, d);
            agg_data_avg(i_subject, :, :, i_targ_side, hit+1) = d_sub.avg;
            agg_data_var(i_subject, :, :, i_targ_side, hit+1) = d_sub.var;
        end
    end
end


% Plot it
side_label = {'left' 'right'};
hit_label = {'miss' 'hit'};
for field = {'avg' 'var'}
    field = field{1};
    switch field
        case 'avg'
            agg = agg_data_avg;
            ylims = [-1 1] * 3e-5;
        case 'var'
            agg = agg_data_var;
            ylims = [0 1] * 5e-9;
    end
    for i_targ_side = 1:2
        for hit = 0:1
            % Array for all data: Subj x Chan x Time x TargSide x Hit 
            x = agg(:, 2, :, i_targ_side, hit+1);
            x = squeeze(x);
            subplot(2, 2, (hit * 2) + i_targ_side)
            plot(d_sub.time, x)
            hold on
            plot(d_sub.time, nanmean(x, 1), '-k', 'LineWidth', 2)            
            hold off
            ylim(ylims)
            title(sprintf('%s %s', side_label{i_targ_side}, hit_label{hit+1}))
        end
    end
    print('-dpng', [exp_dir 'plots/eog/' field])
end

%% Coherence - Left hits vs right hits

close all
clear variables
rs_setup

% Load the data
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    d = load([exp_dir 'coherence/' fname '/coh']);
    % Dim: subject * chancmb * freq * time * RFT_freq_mapping
    x(i_subject,:,:,:,:) = d.x;
end

% Plot it
for i_subject = 0:height(subject_info)
    if i_subject == 0
        subj_sel = 1:height(subject_info);
        fname = 'avg';
    elseif subject_info.exclude(i_subject)
        continue
    else
        subj_sel = i_subject;
        fname = strrep(subject_info.meg{i_subject}, '/', '_');
    end
    
    for i_rft_map = 1:2 % Which side each tagged freq is on
        i_chancmb = 1; % Only using the RESS channel for stims on the right
        x_sub = x(subj_sel, i_chancmb, :, :, i_rft_map);
        x_sub = nanmean(x_sub, 1);
        x_sub = squeeze(x_sub);

        subplot(2, 1, i_rft_map)
        imagesc(d.coh.time, ...
            d.coh.freq, ...
            x_sub)
        set(gca, 'YDir', 'normal')
        title(sprintf('%s, %i HZ on left', ...
            d.coh.labelcmb{i_chancmb}, ...
            exp_params.tagged_freqs(i_rft_map)));
        colorbar;
    end

    fn = [exp_dir 'plots/coherence/left-hits_vs_right-hits/' fname];
    print('-dpng', fn)
end


%% Check whether the photodiode recordings can be approximated with a sine

clear variables
close all
rs_setup

freq_adjust = [-0.0535 -0.0648];

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    behav = rs_behavior(i_subject);
    d = load([exp_dir 'preproc/photodiode/trial/' fname '/preproc']);
    d = d.data;
    cfg = [];
    cfg.detrend = 'yes';
    d = ft_preprocessing(cfg, d);
    
    % Only keep data from when the stimulus was on
    trial_num = d.trialinfo(:,2);
    rt = behav.rt(trial_num);
    for i_trial = 1:length(d.trial)
        t = d.time{i_trial};
        after_stim_onset = t > 0 + 0.1;
        before_stim_end = t < 4 - 0.1;
        before_resp = t < rt(i_trial);
        t_sel = after_stim_onset & before_stim_end & before_resp;
        d.trial{i_trial} = d.trial{i_trial}(t_sel);
        d.time{i_trial} = d.time{i_trial}(t_sel);
    end
    % Toss empty trials
    empty_trial = cellfun(@isempty, d.time);
    cfg = [];
    cfg.trials = ~empty_trial;
    d = ft_selectdata(cfg, d);

    % Compute the phase of the real signal
    hilbphase = @(x) angle(hilbert(x));
    p_data = cellfun(hilbphase, d.trial, 'UniformOutput', false);

    % Check the actual frequency of the stimulus
    freq_measured = nan(size(d.trial));
    for i_trial = 1:length(d.trial)
        unwrapped = unwrap(p_data{i_trial});
        phase_change = unwrapped(end) - unwrapped(1);
        time_change = d.time{i_trial}(end) - d.time{i_trial}(1);
        freq_measured(i_trial) = phase_change / time_change / (2 * pi);
    end

    p_data_all = p_data;
    freq_measured_all = freq_measured;
    d_all = d;
    
    for i_freq = 1:2
        % Look at each frequency separately
        freq = exp_params.tagged_freqs(i_freq);
        keep_trials = abs(freq_measured_all - freq) < 1;
        p_data = p_data_all(keep_trials);
        freq_measured = freq_measured_all(keep_trials);
        cfg = [];
        cfg.trials = keep_trials;
        d = ft_selectdata(cfg, d_all);

        % Plot the distribution of measured frequency
        subplot(3, 2, i_freq)
        bins = 60:0.1:80;
        histogram(freq_measured, bins, 'DisplayStyle', 'stairs')
        hold on

        % Synthesize a sine wave at the frequency of this oscillation
        trial_num = d.trialinfo(:,2);
        freq_right = behav.freq_right(trial_num);
        d_synth = cell(size(d.trial));
        for i_trial = 1:length(d.trial)
            f = freq_right(i_trial);
            %%% Adjust freq using measured diff in freq due to the frame
            %%% rate being slightly slower than 120 Hz
            f = f + freq_adjust(i_freq);
            t = d.time{i_trial};
            amp = max(d.trial{i_trial});
            y = amp * sin(2 * pi * f * t);
            d_synth{i_trial} = y;
        end
        clear f t y

        % Check whether the synthesized signals match the real signals
        p_synth = cellfun(hilbphase, d_synth, 'UniformOutput', false);
        p_diff = cellfun(@minus, p_data, p_synth, 'UniformOutput', false);
        p_diff = cellfun(@wrapToPi, p_diff, 'UniformOutput', false);

        % Plot the phase difference between the real and synthesized
        % signals
        subplot(3, 2, 2 + i_freq)
        histogram(cat(2, p_diff{:}), 250, 'DisplayStyle', 'stairs')
        hold on
        
        % Plot phase difference as a function of time in the trial
        subplot(3, 2, 4 + i_freq)
        t_vec = cat(2, d.time{:});
        p_vec = cat(2, p_diff{:});
        s = randsample(length(t_vec), 1000);
        scatter(t_vec(s), p_vec(s))
        hold on
        
        % Fit a line to the phase lag over time
        % This will tell us what the real frequency is
        if i_freq == 1 % Toss weird values due to dropped frames
            keep_obs = (-2 < p_vec) & (p_vec < 1); 
        elseif i_freq == 2
            keep_obs = (-pi < p_vec) & (p_vec < -1); 
            p_vec = wrapTo2Pi(p_vec);
        end
        p_vec(~keep_obs) = [];
        t_vec(~keep_obs) = [];
        X = [ones(length(p_vec), 1) t_vec'];
        beta = X\p_vec';
        b(:,i_subject, i_freq) = beta';
        % Units of b(2,:) - rad/s  
        % So the drift for each subject is 
    end
end

for i_plot = 1:2
    subplot(3,2,i_plot)
    title('Measured frequency')
    xlabel('Frequency (Hz)')
    ylabel('Count')
    hold off
end

for i_plot = 3:4
    subplot(3,2,i_plot)
    title('Phase: measured vs synth')
    xlim([-pi pi])
    xlabel('Phase difference (rad)')
    ylabel('Count')
    hold off
end

for i_plot = 5:6
    subplot(3,2,i_plot)
    title('Phase diff by time')
    ylim([-pi pi])
    ylabel('Phase difference (rad)')
    xlabel('Time (s)')
    hold off
end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 8];
print('-dpng', [exp_dir 'plots/coherence/phase_test'])

% Get the difference in Hz between expected and displayed signals
b(b==0) = NaN;
beta = squeeze(b(2,:,:)); % Get slopes
beta = beta * (1 / (2 * pi)); % Convert to Hz
nanmean(beta)
