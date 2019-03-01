function d = rs_powerratio(i_subject)

% Get the ratio of normalized power at each side
%
% Power time-courses
% HP filter
% Cut out transients
% Z-score
% Ratio of Power(63) / Power(78) as a function of time, compare hits vs misses

rs_setup

win_size = 0.2;
win_str = sprintf('win_%.1fs', win_size);
tfr_dir = [exp_dir 'tfr/' win_str '/target/'];

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
        hp_freq = 1; % Hz
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
                    pwr_filt(i_trial,i_channel,i_freq,:) = filtfilt(b,a,x);
                end
            end
        end 
        
        % Cut out transients (don't have to worry about that for target-seg)

        % Z-score power in each frequency separately
        z = nan(size(pwr));
        for i_freq = 1:length(d.freq)
            for i_chan = 1:length(d.label)
                x = pwr(:,i_chan,i_freq,:);
                x_mean = nanmean(reshape(x, [1 numel(x)]));
                x_std = nanstd(reshape(x, [1 numel(x)]));
                z(:,i_chan,i_freq,:) = (x - x_mean) / x_std;
            end
        end
        
        % TESTING
        %{
        f = d_sub.freq;
        t = d_sub.time;
        plt = @(x) imagesc(t, f, squeeze(nanmean(x(:,1,:,:), 1)));
        subplot(3, 1, 1)
        plt(pwr), set(gca, 'YDir', 'normal')
        subplot(3, 1, 2)
        plt(pwr_filt), set(gca, 'YDir', 'normal')
        subplot(3, 1, 3)
        plt(z), set(gca, 'YDir', 'normal')
        %}

        % Compute the ratio of Ipsi:Contra power
        warning('Make sure these ipsi/contra labels are correct')
        ipsi_timecourse = z(:,i_targ_side,i_targ_freq,:);
        k_targ_side = mod(i_targ_side, 2) + 1; % Get the non-target side
        k_targ_freq = mod(i_targ_freq, 2) + 1;
        contra_timecourse = z(:,k_targ_side,k_targ_freq,:);
        power_ratio = ipsi_timecourse ./ contra_timecourse;

        % Put this back into the main data object
        d.powerratio(trial_sel,:) = squeeze(power_ratio);

    end
end

% Ratio of Power(63) / Power(78)???


