% Compare power at the tagged frequency for hits on the L vs hits on the R.

clear variables
rs_setup
sides = {'left' 'right'};

segment_type = 'target'; % target | trial
win_size = 0.1;

win_str = sprintf('win_%.1fs', win_size);
tfr_dir = [exp_dir 'tfr/' win_str '/' segment_type '/'];

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
        stim_side_sel = strcmp(behav.target_side(trial_num), 'right');
        stim_freq_sel = behav.target_side_freq(trial_num) == dist_freq;
        cfg.trials = hit_sel & stim_side_sel & stim_freq_sel;
        d_r = ft_selectdata(cfg, d);
        x_r = d_r.powspctrm;

        % Compare them
        x(i_subject,:,i_tagged_freq,:) = (x_l - x_r) ./ (x_l + x_r);
    end
end


%% Plot it

for i_sensor_side = 1:2
    % RESS side is coded as _stimulus_ side. Flip around to get to the
    % _sensor_ side.
    i_ress_side = mod(i_sensor_side, 2) + 1;
    for i_tagged_freq = 1:2
        subplot(2, 2, (i_tagged_freq - 1) * 2 + i_sensor_side)
        z = squeeze(x(:,i_ress_side,i_tagged_freq,:));
        z_mean = nanmean(z, 1);
        plot(d_l.time, z, 'color', 0.7 * [1 1 1])
        hold on
        plot([-0.5 0.5], [0 0], '-k')
        plot([0 0], [min(min(z)) max(max(z))], '-k')
        plot(d_l.time, z_mean, '-r', 'LineWidth', 2)
        title(sprintf('%s, %i Hz', ...
            sides{i_sensor_side}, exp_params.tagged_freqs(i_tagged_freq)));
        ylim([min(min(z)) max(max(z))])
        % Simple stats
        [h, p] = ttest(z, 0, 'dim', 1);
        plot(d_l.time(p < 0.05), zeros([1 sum(h)]), '*b')
        hold off
    end
end

print('-dpng',...
    [exp_dir 'plots/accuracy/high_freq/left-hits_vs_right-hits'])