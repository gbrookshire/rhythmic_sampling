% For looking at effects of LF phase, should this be done on RESS virtual
% channels or on the full data? Probably on the full data.
%
% TODO
% - run on the raw channels, not RESS virtual channels
% - separate by hemisphere
%   - Or just compute separately for targets on each side of space?
%   - Maybe first, check correlation of LF phase across hemispheres?
% - Replace 'mean' with phase consistency (or something)


clear variables
close all
rs_setup


% Array for all data: Subj x Chan x Freq x Time x TargetSide x Hit
agg_data = nan([height(subject_info), 2, 11, 51, 2, 2]);

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end

    % Read in the data segmented around targets
    behav = rs_behavior(i_subject);
    fname = subject_info.meg{i_subject};
    fn = [exp_dir 'tfr/target/' fname '/low'];
    d = load(fn);
    d = d.low_freq_data;

    % Information about each trial
    trial_numbers = d.trialinfo(:,2);
    targ_side_per_trial = behav.target_side(trial_numbers);
    targ_freq_per_trial = behav.target_side_freq(trial_numbers);

    % Put the data together
    targ_side_labels = {'left' 'right'};
    for i_side = 1:2 % Which side did the target appear on?
        targ_side = targ_side_labels{i_side};
        sel_targ_side = strcmp(targ_side_per_trial, targ_side);
        for hit = [0 1]
            sel_hit = d.trialinfo(:,1) == hit;
            trial_sel = sel_hit & sel_targ_side;
            cfg = [];
            cfg.trials = trial_sel;
            d_sub = ft_selectdata(cfg, d);
            x = PHASECONSISTENCY(d_sub.fourierspctrm);
            agg_data(i_subject,:,:,:,i_side,hit+1) = x;

end
