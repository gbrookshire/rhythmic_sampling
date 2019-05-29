% Compare power at the tagged frequency for hits on the L vs hits on the R.

clear variables
rs_setup
sides = {'left' 'right'};

segment_type = 'target'; % target | trial
w0in_size = 0.1;

% win_str = sprintf('win_%.1fs', win_size);
vers = 'raw_chans';
tfr_dir = [exp_dir 'tfr/' vers '/' segment_type '/'];

clear x
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        x(i_subject,:,:,:) = nan;
        continue
    end
    
    % Load target-segmented HF TFR
    % Read in the data segmented around targets
    fname = subject_info.meg{i_subject};
    fn = [tfr_dir fname '/high'];
    d = load(fn);
    d = d.high_freq_data;

    % Load the behavioral data
    behav = rs_behavior(i_subject);

    % Which trials were hits
    hit_sel = d.trialinfo(:,1) == 1;
    trial_num = d.trialinfo(:, 2);

    for i_tagged_freq = 1:2
        % Get the frequencies of the target and distractor
        targ_freq = exp_params.tagged_freqs(i_tagged_freq);
        dist_freq = exp_params.tagged_freqs(mod(i_tagged_freq, 2) + 1);
        
        % Which data to select
        cfg = [];
        cfg.frequency = targ_freq + [-0.1 0.1];
        cfg.avgoverrpt = 'yes';
        cfg.nanmean = 'yes';
        
        % Hits on the LEFT with the target at THIS freq
        stim_side_sel = strcmp(behav.target_side(trial_num), 'left');
        stim_freq_sel = behav.freq_left(trial_num) == targ_freq;
        cfg.trials = hit_sel & stim_side_sel & stim_freq_sel;
        d_l = ft_selectdata(cfg, d);
        x_l = d_l.powspctrm;
        
        % Hits on the RIGHT with the target at the OTHER freq
        % (To ensure that the sensors are contalateral to the same freq)
        stim_side_sel = strcmp(behav.target_side(trial_num), 'right');
        stim_freq_sel = behav.freq_right(trial_num) == dist_freq;
        cfg.trials = hit_sel & stim_side_sel & stim_freq_sel;
        d_r = ft_selectdata(cfg, d);
        x_r = d_r.powspctrm;

        % Compare them
        x(i_subject,:,i_tagged_freq,:) = (x_l - x_r) ./ (x_l + x_r);
    end
end


%% Plot it

for i_sensor_side = 1:2
    
    %{
    % RESS side is coded as _stimulus_ side. Flip around to get to the
    % _sensor_ side.
    i_ress_side = mod(i_sensor_side, 2) + 1;
    chan_inx = i_ress_side;
    %}
    
    % Select channels for the raw-channel (non-RESS) analysis
    if i_sensor_side == 1
        chan_inx = 3:4; % Sensors on the left side
    elseif i_sensor_side == 2
        chan_inx = 1:2;
    end
    
    for i_tagged_freq = 1:2
        subplot(2, 2, (i_tagged_freq - 1) * 2 + i_sensor_side)
        mi = squeeze(mean(x(:,chan_inx,i_tagged_freq,:), 2));
        mi_mean = nanmean(mi, 1);
        plot(d_l.time, mi, 'color', 0.7 * [1 1 1])
        hold on
        plot([-0.5 0.5], [0 0], '-k')
        plot([0 0], [min(min(mi)) max(max(mi))], '-k')
        plot(d_l.time, mi_mean, '-r', 'LineWidth', 2)
        title(sprintf('%s, %i Hz', ...
            sides{i_sensor_side}, exp_params.tagged_freqs(i_tagged_freq)));
        ylim([min(min(mi)) max(max(mi))])
        % Simple stats
        [h, p] = ttest(mi, 0, 'dim', 1);
        plot(d_l.time(p < 0.05), zeros([1 sum(h)]), '*b')
        hold off
    end
end

print('-dpng',...
    [exp_dir 'plots/accuracy/high_freq/left-hits_vs_right-hits'])

%% Plot the difference between the right and left sides

% % For selecting RESS channels
% left_inx = 2;
% right_inx = 1;

% For selecting raw channels
left_inx = 3:4;
right_inx = 1:2;

for i_tagged_freq = 1:2
    subplot(2, 1, i_tagged_freq)
    mi_left = squeeze(mean(x(:, left_inx, i_tagged_freq, :), 2));
    mi_right = squeeze(mean(x(:, right_inx, i_tagged_freq, :), 2));
    mi_diff = mi_right - mi_left;
    mi_diff_mean = nanmean(mi_diff, 1); 
    plot(d_l.time, mi_diff, 'color', 0.7 * [1 1 1])
    hold on
    plot([-0.5 0.5], [0 0], '-k')
    plot([0 0], [min(min(mi_diff)) max(max(mi_diff))], '-k')
    plot(d_l.time, mi_diff_mean, '-r', 'LineWidth', 2)
    title(sprintf('%i Hz', ...
        exp_params.tagged_freqs(i_tagged_freq)));
    ylim([min(min(mi_diff)) max(max(mi_diff))])
    % Simple stats
    [h, p] = ttest(mi_diff, 0, 'dim', 1, 'alpha', 0.01);
    plot(d_l.time(logical(h)), zeros([1 sum(h)]), '*b')
    hold off 
end

print('-dpng',...
    [exp_dir 'plots/accuracy/high_freq/left-hits_vs_right-hits_diff'])


