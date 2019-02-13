rs_setup

i_subject = 1;

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


%% Plot it

i_plot = 1;
for targ_side = {'left' 'right'}
    sel_targ_side = strcmp(targ_side_per_trial, targ_side);
    for targ_freq = exp_params.tagged_freqs
        sel_targ_freq = targ_freq_per_trial == targ_freq;
        
        for hit = [0 1]
            sel_hit = d.trialinfo(:,1) == hit;
            trial_sel = sel_targ_side & sel_targ_freq & sel_hit;
            cfg = [];
            cfg.trials = trial_sel;
            cfg.channel = targ_side; %%% keep same or other side RESS channel?
            d_sub = ft_selectdata(cfg, d);
            
            subplot(4, 2, i_plot)
            imagesc(d_sub.time, d_sub.freq, ...
                squeeze(nanmean(d_sub.powspctrm, 1)))
            set(gca, 'YDir', 'normal')
            title(sprintf('%s, %i Hz, Hit:%i', ...
                targ_side{1}, targ_freq, hit))
            
            i_plot = i_plot + 1;
        end
    end
end

%% Plot them each
i_plot = 1;
for data = {hits misses}
    subplot(3, 1, i_plot)
    imagesc(data{1}.time, data{1}.freq, ...
        squeeze(mean(nanmean(data{1}.powspctrm, 1), 2)))
    set(gca, 'YDir', 'normal')
    i_plot = i_plot + 1;
end

diff = hits;
diff.powspctrm = nanmean(hits.powspctrm,1) - nanmean(misses.powspctrm,1);
subplot(3, 1, 3)
imagesc(diff.time, diff.freq, ...
    squeeze(mean(nanmean(diff.powspctrm, 1), 2)))
set(gca, 'YDir', 'normal')
