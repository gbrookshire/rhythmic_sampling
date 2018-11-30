function rs_accuracy_lfphase(i_subject)

% Replicate analysis from Fiebelkorn et al (2018, Neuron)

% Load LF TFR for one subject
%   - Complex Morlet wavelets
%   - 3-10 Hz: 2 cycles; 11-14 Hz: 3 cyc; 15-20 Hz: 4 cyc; >20 Hz: 5 cyc
% Select a time "just prior to target presentation".
%   - Target doesn't overlap with this window
% Take the angle of the complex output
% Make phase-bins that are 180 deg wide
% Calculate the hit rate in that 180 deg phase bin
% Shift the bin by 5 deg and do the same
% Take the ouput of this shifting phase-bin procedure
% Apply an FFT to it
% Take the 2nd component, which represents a 1-cycle sine wave
% The amplitude of that sine wave is your effect for this freq

rs_setup

% Read in the data segmented around targets
fname = subject_info.meg{i_subject};
fn = [exp_dir 'tfr/target/' fname '/low'];
d = load(fn);
d = d.low_freq_data;

% Select only hits and misses
[hits, nans] = rs_resptype(i_subject);

% Exclude the trials with NaNs in the target trialdef
hits = hits(~nans);
hits_and_misses_inx = ismember(hits, [0 1]);

% Keep only the hits and misses (no FAs or late responses)
cfg = [];
cfg.trials = hits_and_misses_inx;
d = ft_selectdata(cfg, d);
hit = hits(hits_and_misses_inx);
clear hits nan

% Compute the phase angle
phase = angle(d.fourierspctrm);

% Get the latest phase estimate that doesn't overlap with the target
tfr_window = d.cfg.previous.t_ftimwin;
t_choice = nan([1 length(d.freq)]); % Which TFR window to use for each freq
for i_freq = 1:length(d.freq)
    % Get the end times for each TFR bin
    window_end = d.time + (tfr_window(i_freq) / 2);
    % Get rid of windows that include the target (end after t=0)
    window_end(window_end >= 0) = [];
    % The latest window that doesn't overlap the target
    t_choice(i_freq) = length(window_end);
end

%{
% TEST - add hits only when 0 < phase < pi/2
x_freq = 5;
x_chan = 224;
for i_trial = 1:size(d.fourierspctrm, 1)
    phi = phase(i_trial, x_chan, x_freq, t_choice(x_freq));
    if 0 < phi && phi < pi/2
        hit(i_trial) = 1;
    else
        hit(i_trial) = 0;
    end
end
%}

% Get the binned hit rate
% Get the average hit rate in bins of LF phase
phase_bin_width = 180; % Degrees
phase_bin_step_size = 5; 
phase_bin_width = deg2rad(phase_bin_width);
phase_bin_step_size = deg2rad(phase_bin_step_size); 
n_bins =  2 * pi / phase_bin_step_size;
hit_rate = nan([... % Bin * Channel * Freq
    n_bins, ...
    length(d.label), ...
    length(d.freq)]);
n_trials = hit_rate;
for i_bin = 1:n_bins % Takes a minute or two
    % Get the difference between MEG phase and center of the bin
    bin_center = -pi + (i_bin * phase_bin_step_size);
    complex_diff = phase - bin_center;
    theta_diff = abs(angle(exp(1i * complex_diff)));
    for i_freq = 1:length(d.freq)
        tfr_window_inx = t_choice(i_freq); % TFR window just before target
        for i_chan = 1:length(d.label)
            % Keep trials that are within width/2 of the bin center
            theta = theta_diff(:, i_chan, i_freq, tfr_window_inx);
            trials_in_bin = theta < (phase_bin_width / 2);
            hr = mean(hit(trials_in_bin));
            hit_rate(i_bin, i_chan, i_freq) = hr;
            n_trials(i_bin, i_chan, i_freq) = sum(trials_in_bin);
        end
    end
end

% Run an FFT to get the strength with which it's modulated by phase
y = fft(hit_rate, n_bins, 1);
z = squeeze(abs(y(2,:,:))); % Get periodicity at the freq of interest

save([exp_dir 'tfr/target/' fname '/high_power_acc_stats'], ...
    'hit_rate', 'n_trials', 'z')

% % Plot it
% roi = ismember(d.label, occip_roi);
% plot(d.freq, z(roi,:))
