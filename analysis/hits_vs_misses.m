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

    % Only look at hits and misses (no FAs or late responses)
    cfg = [];
    cfg.trials = d.trialinfo(:,1) == 1;
    hits = ft_selectdata(cfg, d);
    cfg.trials = d.trialinfo(:,1) == 0;
    misses = ft_selectdata(cfg, d);

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

%% plot it

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
        imagesc(d_sub.time, d_sub.freq, ...
            hit_miss(:,:,2) - hit_miss(:,:,1))
        set(gca, 'YDir', 'normal')
        title(sprintf('%s, %i Hz, Diff', ...
            targ_side, targ_freq))

        i_plot = i_plot + 1;
    end
end

%% Plot it, collapsing over side
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
    imagesc(d_sub.time, d_sub.freq, ...
        hit_miss(:,:,2) - hit_miss(:,:,1))
    set(gca, 'YDir', 'normal')
    title(sprintf('%i Hz, Diff', ...
        targ_freq))

    i_plot = i_plot + 1;
end

