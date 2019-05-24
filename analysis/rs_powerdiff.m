function out = rs_powerdiff(i_subject, win_size, segment_type, k_bootstrap)

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

% i_subject = 1;
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

% Z-score power in each channel and frequency separately
pwr = d.powspctrm;
z = nan(size(pwr));
for i_freq = 1:length(d.freq)
    for i_chan = 1:length(d.label)
        x = pwr(:,i_chan,i_freq,:);
        x_mean = nanmean(reshape(x, [1 numel(x)]));
        x_std = nanstd(reshape(x, [1 numel(x)]));
        z(:,i_chan,i_freq,:) = (x - x_mean) ./ x_std;
    end
end
clear pwr

% Load the behavioral data
behav = rs_behavior(i_subject);

% Main analysis
sides = {'left' 'right'};
powdiff_hit = nan(0, length(d.time));
powdiff_miss = nan(0, length(d.time));

for i_targ_side = 1:2
        
    % Select trials with the target on this side
    targ_side = sides{i_targ_side};
    trial_num = d.trialinfo(:, 2);
    targ_side_sel = strcmp(targ_side, behav.target_side(trial_num));

    for i_targ_freq = 1:2
        targ_freq = exp_params.tagged_freqs(i_targ_freq);
        
        for hit = 0:1
            hit_sel = d.trialinfo(:,1) == hit;

            % Select targets that appeared in stimuli tagged at this freq
            targ_freq_sel = targ_freq == behav.target_side_freq(trial_num);

            % Get the power timecourse
            trial_sel = targ_side_sel & targ_freq_sel & hit_sel;
            z_sub = z(trial_sel,:,:,:);
            
%             % HP filter
%             warning('HP filtering power time-courses at 2 Hz!')
%             hp_freq = 2; % Hz
%             filt_order = 5;
%             fsample = mean(diff(d.time)) ^ -1;
%             nyq = fsample / 2;
%             [b, a] = butter(filt_order, hp_freq / nyq, 'high');
%             pwr_filt = nan(size(d_sub.powspctrm));
%             for i_trial = 1:size(pwr, 1)
%                 for i_channel = 1:length(d.label)
%                     for i_freq = 1:length(d.freq)
%                         x = squeeze(pwr(i_trial,i_channel,i_freq,:));
%                         x = filtfilt(b, a, x);
%                         pwr_filt(i_trial,i_channel,i_freq,:) = x;
%                     end
%                 end
%             end
%             clear pwr
        
            % Cut out transients (don't have to worry about that for target-seg)
           
            % Make a bootstrap distribution of these values by randomly
            % shuffling the two channels between trials. Do this within
            % hits/misses.
            if k_bootstrap > 0
                size_z = size(z_sub);
                z_old = z_sub;
                perm1 = randsample(size_z(1), k_bootstrap, true);
                perm2 = randsample(size_z(1), k_bootstrap, true);
                z_sub = z_old(perm1, 1, :, :);
                z_sub(:,2,:,:) = z_old(perm2, 2, :, :);        
            end

            % Index the target and non-target frequencies
            rft_freqs = exp_params.tagged_freqs;
            f_t_inx = abs(d.freq - rft_freqs(i_targ_freq)) < 0.5;
            f_n_inx = abs(d.freq - rft_freqs(mod(i_targ_freq, 2) + 1)) < 0.5;
            
            % Compare power in the target and non-target stimuli
            % Get the timecourse of power at the target stimulus
            x_targ = z_sub(:,i_targ_side,f_t_inx,:);
            % Get the timecourse of power at the non-target stimulus
            k_targ_side = mod(i_targ_side, 2) + 1;
            x_nontarg = z_sub(:,k_targ_side,f_n_inx,:);
            % Compute difference in power
            power_diff = x_targ - x_nontarg;
            power_diff = squeeze(power_diff);

            if hit == 0
                powdiff_miss = cat(1, powdiff_miss, power_diff);
            elseif hit == 1
                powdiff_hit = cat(1, powdiff_hit, power_diff);
            end

        end
    end
end

out = [];
out.time = d.time;
out.powdiff_hit = powdiff_hit;
out.powdiff_miss = powdiff_miss;
