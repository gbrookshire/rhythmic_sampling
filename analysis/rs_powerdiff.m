function out = rs_powerdiff(i_subject, win_size, segment_type)

% Get the ratio of normalized power at each side
%
% Power time-courses
% HP filter
% Cut out transients
% Z-score
% Compare power in target-nontarget tagged frequencies as a func of time
%   Calculate the difference (subtraction)
%       Ratio leads to instability b/c power is 0-centered after HP filter
%   compare hits vs misses

rs_setup

% segment_type = 'target'; % target | trial
% win_size = 0.2;

win_str = sprintf('win_%.1fs', win_size);
tfr_dir = [exp_dir 'tfr/' win_str '/' segment_type '/'];

% Load target-segmented HF TFR
% Read in the data segmented around targets
fname = subject_info.meg{i_subject};
fn = [tfr_dir fname '/high'];
d = load(fn);
d = d.high_freq_data;

% Load the behavioral data
behav = rs_behavior(i_subject);

% Main analysis
sides = {'left' 'right'};
powdiff = nan(size(d.powspctrm, 1), length(d.time));
for i_targ_side = 1:2
        
    % Select trials with the target on this side
    targ_side = sides{i_targ_side};
    trial_num = d.trialinfo(:, 2);
    targ_side_sel = strcmp(targ_side, behav.target_side(trial_num));

    for i_targ_freq = 1:2

        % Select targets that appeared in stimuli tagged at this freq
        targ_freq = exp_params.tagged_freqs(i_targ_freq);
        targ_freq_sel = targ_freq == behav.target_side_freq(trial_num);

        % Get the power timecourse
        trial_sel = targ_side_sel & targ_freq_sel;
        cfg = [];
        cfg.trials = trial_sel;
        d_sub = ft_selectdata(cfg, d);

        % HP filter 
        hp_freq = 2; % Hz
        filt_order = 5;
        fsample = mean(diff(d.time)) ^ -1;
        nyq = fsample / 2;
        [b, a] = butter(filt_order, hp_freq / nyq, 'high');
        pwr = d_sub.powspctrm;
        pwr_filt = nan(size(d_sub.powspctrm));
        for i_trial = 1:size(pwr, 1)
            for i_channel = 1:length(d.label)
                for i_freq = 1:length(d.freq)
                    x = squeeze(pwr(i_trial,i_channel,i_freq,:));
                    x = filtfilt(b, a, x);
                    pwr_filt(i_trial,i_channel,i_freq,:) = x;
                end
            end
        end
        clear pwr
        
        % Cut out transients (don't have to worry about that for target-seg)

        % Z-score power in each frequency separately
        z = nan(size(pwr_filt));
        for i_freq = 1:length(d.freq)
            for i_chan = 1:length(d.label)
                x = pwr_filt(:,i_chan,i_freq,:);
                x_mean = nanmean(reshape(x, [1 numel(x)]));
                x_std = nanstd(reshape(x, [1 numel(x)]));
                z(:,i_chan,i_freq,:) = (x - x_mean) ./ x_std;
            end
        end
        clear pwr_filt
        
        % TESTING
        %{
        f = d_sub.freq;
        t = d_sub.time;
        plt = @(a) imagesc(t, f, squeeze(nanmean(a(:,1,:,:), 1)));
        subplot(3, 1, 1)
        plt(pwr), set(gca, 'YDir', 'normal'), colorbar;
        subplot(3, 1, 2)
        plt(pwr_filt), set(gca, 'YDir', 'normal'), colorbar;
        subplot(3, 1, 3)
        plt(z), set(gca, 'YDir', 'normal'), colorbar;
        %}

        % Compare power in the target and non-target stimuli
        % Get the timecourse of power at the target stimulus
        target_timecourse = z(:,i_targ_side,i_targ_freq,:);
        % Get the timecourse of power at the non-target stimulus
        k_targ_side = mod(i_targ_side, 2) + 1;
        k_targ_freq = mod(i_targ_freq, 2) + 1;
        nontarget_timecourse = z(:,k_targ_side,k_targ_freq,:);
        % Compute difference in power
        power_diff = target_timecourse - nontarget_timecourse;
        power_diff = squeeze(power_diff);
        %{
        % Don't use ratio -- because of HP filt and z-score, power is
        % centered around 0, so computing ratios leads to weird numbers
        power_ratio = target_timecourse ./ nontarget_timecourse;
        power_ratio = squeeze(power_ratio);
        %}
        
        % TESTING
        %{
        subplot(2, 1, 1)
        plot(t, power_ratio)
        subplot(2, 1, 2)
        plot(t, power_diff)
        %}

        % Put this back into the main data object
        powdiff(trial_sel,:) = power_diff;
    end
end

out = [];
out.time = d.time;
out.powdiff = powdiff;
out.trialinfo = d.trialinfo;
