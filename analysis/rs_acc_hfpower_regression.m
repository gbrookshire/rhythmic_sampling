function rs_acc_hfpower_regression(i_subject)

% Does accuracy correlate with power at the tagged freq?
% This function runs VERY slowly (~4.5 hrs / subj). Find a different way to
% do the analysis?

rs_setup

devcode = @(x) x - 0.5; % Move DV from 0/1 to +/-0.5 regression coding

% Read in the data segmented around targets
fname = subject_info.meg{i_subject};
fn = [exp_dir 'tfr/target/' fname '/high'];
d = load(fn);
d = d.high_freq_data;

% Load the behavioral data
behav = rs_behavior(i_subject);
behav = behav(225:end, :); % main blocks of trials

% Select only hits and misses
[hits, nans] = rs_resptype(i_subject);

% Exclude the trials with NaNs in the target trialdef
hits = hits(~nans);
behav = behav(~nans,:);

% Only look at hits and misses (no FAs or late responses)
only_hits_and_misses = ismember(hits, [0 1]);
hits = hits(only_hits_and_misses);
behav = behav(only_hits_and_misses,:);
keep_trials = ismember(d.trialinfo(:,2), behav.TrialNumber);
cfg = [];
cfg.trials = keep_trials;
d = ft_selectdata(cfg, d);

% Get whether each data segment was a hit or miss
hit_inx

% Keep only the hits and misses (no FAs or late responses)
hits_and_misses_inx = ismember(hits, [0 1]);

d.trialinfo(

cfg = [];
cfg.trials = hits_and_misses_inx;
d = ft_selectdata(cfg, d);
behav = behav(hits_and_misses_inx,:);
hit = hits(hits_and_misses_inx);
clear hits nan

% The following code is extremely slow (est 4.5 hrs for 1 subj)
% Run a regression for each (time,freq) pair
stats = cell([length(d.time), length(d.freq), length(d.label)]);
for i_t = 1:length(d.time)
    for i_f = 1:length(d.freq)
        % Extract power at this (time,freq) point
        % Make a matrix of Trial x Channel
        pwr = nan([size(d.powspctrm, 1) length(d.label)]);
        for i_trial = 1:size(d.powspctrm, 1)
            pwr(i_trial,:) = squeeze(d.powspctrm(i_trial,:,i_f,i_t));
        end
        % Run a regression on power and hit/miss, controlling for side of
        % screen and frequency.
        % This loop takes ~7 sec. Any way to speed it up?
        for i_chan = 1:length(d.label)
            % Set up variables for a regression
            tbl = table(normalize(pwr(:, i_chan)), ...
                devcode(hit'), ...
                devcode(strcmp(behav.target_side, 'right')), ...
                devcode(behav.target_side_freq == 78), ...
                'VariableNames', {'Power', 'Hit', 'Side', 'Freq'});
            m = fitglm(tbl, ...
                'Power ~ Hit + Freq + Side');
            stats{i_t, i_f, i_chan} = m.Coefficients;
        end
    end
end

time = d.time;
freq = d.freq;
label = d.label;
dimord = 'time_freq_chan';
save([exp_dir 'tfr/target/' fname '/hfpower_acc_stats'], ...
    'stats', 'time', 'freq', 'label', 'dimord')

% Extract p-values like this:
% pvals = cellfun(@(c) c{'Hit', 'pValue'}, stats);
