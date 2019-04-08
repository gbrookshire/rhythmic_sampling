function rs_acc_lfphase_fiebelkorn(i_subject)

% Replicate analysis from Fiebelkorn et al (2018, Neuron)
%
% Output:
%   hit_rate: Array of hit rate (Phase Bin * Channel * Freq)
%   n_trials: Number of trials in each cell (same dims as hit_rate)
%   z: Amplitude of oscillations at each freq (Channel * Freq)
%   n_bins: How many phase bins were used
%   label: Channel labels
%   freq: Frequencies at which accuracy was calculated
%   dimord: string describing dimensions

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
fn = [exp_dir 'tfr/win_0.1s/target/' fname '/low'];
d = load(fn);
d = d.low_freq_data;

% Read in CSV to get which side target was on
behav = rs_behavior(i_subject);

% Select only hits and misses
cfg = [];
cfg.trials = ismember(d.trialinfo(:,1), [0 1]);
d = ft_selectdata(cfg, d);

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

hit = d.trialinfo(:,1); % Was each trial a hit or miss
targ_right = strcmp('right', behav.target_side(d.trialinfo(:,2)));

%{
% TEST - add hits only when 0 < phase < pi/2
fname = 'SIMULATED';
x_freq = 5; % Add hits as a function of phase at this freq
x_chan = 200; % ... calculated at this MEG channel
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
hit_rate = nan([... % Phase Bin * Target Side * Channel * Freq
    n_bins, ...
    2, ...
    length(d.label), ...
    length(d.freq)]);
n_trials = hit_rate;

for i_bin = 1:n_bins % Takes a minute or two
    % Get the difference between MEG phase and center of the bin
    bin_center = -pi + (i_bin * phase_bin_step_size);
    complex_diff = phase - bin_center;
    theta_diff = abs(angle(exp(1i * complex_diff)));
    for i_freq = 1:length(d.freq)
        tfr_win_inx = t_choice(i_freq); % TFR wind just before targ
        for i_chan = 1:length(d.label)
            for targ_side = 1:2 % 1:left, 2:right
                trl_inx = targ_right == (targ_side - 1);
                % Keep trials that are within width/2 of the bin center
                theta = theta_diff(trl_inx, i_chan, i_freq, tfr_win_inx);
                trials_in_bin = theta < (phase_bin_width / 2);
                hr = mean(hit(trials_in_bin));
                hit_rate(i_bin, targ_side, i_chan, i_freq) = hr;
                trl_count = sum(trials_in_bin);
                n_trials(i_bin, targ_side, i_chan, i_freq) = trl_count;
            end
        end
    end
end

% Run an FFT to get the strength with which it's modulated by phase
y = fft(hit_rate, n_bins, 1);
z = squeeze(abs(y(2,:,:))); % Get periodicity at the freq of interest

label = d.label;
freq = d.freq;
dimord = 'phasebin_chan_freq';
save([exp_dir 'tfr/win_0.1s/target/' fname '/lfphase_acc_fieb'], ...
    'hit_rate', 'n_trials', 'z', 'n_bins', 'label', 'freq', 'dimord')