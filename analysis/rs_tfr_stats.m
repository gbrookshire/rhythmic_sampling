function rs_tfr_stats(i_subject)

% Does accuracy correlate with power at the tagged freq?
% This function runs VERY slowly (~4.5 hrs / subj). Find a different way to
% do the analysis?

rs_setup

devcode = @(x) x - 0.5; % Move DV from 0/1 to +/-0.5 regression coding

% Load in a sample 
i_subject = 15;

% Read in the data segmented around targets
fname = subject_info.meg{i_subject};
fn = [exp_dir 'tfr/target/' fname '/high']; %%%% CHANGE TO TARGET
d = load(fn);
d = d.high_freq_data;

% Load the behavioral data
fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
behav = rs_behavior(fn);
behav = behav(225:end, :); % main blocks of trials

% Select only trials that are hits and misses
if any(isnan(behav.hit))
    trial_inx = ~isnan(behav.hit);
    behav = behav(trial_inx, :);
    cfg = [];
    cfg.trials = trial_inx;
    d = ft_selectdata(cfg, d);
    clear trial_inx cfg
end

% The following code is extremely slow (est 4.5 hrs for 1 subj)
% Run a regression for each (time,freq) pair
stats = cell([length(d.time), length(d.freq), length(d.label)]);
for i_t = 1:length(d.time)
    for i_f = 1:length(d.freq)
        % Extract power at this (time,freq) point
        % Make a matrix of Trial x Channel
        pwr = nan([length(d.trialinfo) length(d.label)]);
        for i_trial = 1:length(d.trialinfo)
            pwr(i_trial,:) = squeeze(d.powspctrm(i_trial,:,i_f,i_t));
        end
        % Run a regression on power and hit/miss, controlling for side of
        % screen and frequency.
        % This loop takes ~7 sec. Any way to speed it up?
        for i_chan = 1:length(d.label)
            % Set up variables for a regression
            tbl = table(normalize(pwr(:, i_chan)), ...
                devcode(behav.hit), ...
                devcode(strcmp(behav.target_side, 'right')), ...
                devcode(behav.target_side_freq == 78), ...
                'VariableNames', {'Power', 'Hit', 'Side', 'Freq'});
            m = fitglm(tbl, ...
                'Power ~ Hit + Freq + Side');
            stats{i_t, i_f, i_chan} = m.Coefficients;
        end
    end
end

save([exp_dir 'tfr/trial/' fname '/high_hit_lm'])

% Extract p-values like this:
% pvals = cellfun(@(c) c{'Hit', 'pValue'}, stats);
