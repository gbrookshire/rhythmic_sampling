function data = rs_simulate_flicker()

% Simulate the kind of data we want to see, and make sure out analysis
% pipeline can reconstruct it.
%
% Simulation details
% 2 sources in L and R visual cortex
% One oscillating at 63 and the other at 78 Hz
% Both modulated at 4-5 Hz in antiphase

rs_setup

% Load a gradiometer structure
grad = load([exp_dir 'grad/181009_b46d/181009/grad']);

% Set the random seed
rng(1)

fsample = 1000;
t = -0.5:(1/fsample):4.5;
n_trials = 336;

% Simulate the flicker frequency modulated at 5 Hz
stim_on = zeros(size(t)); % When the stimulus is flickering
stim_on(t >= 0 & t <= 4) = 1;
mod_freq = 5;
amp = 1;
sig = nan([2 length(t)]); % Channel x Time x Trial
for i_trial = 1:n_trials
    mod_phase = rand(1) * 2 * pi; % Random phase of modulation
    for i_freq = 1:length(exp_params.tagged_freqs)
        % Carrier - flicker freq
        car_freq = exp_params.tagged_freqs(i_freq);
        car_sig = amp * sin(2 * pi * car_freq * t);
        % Modulation - Theta osc
        mod_phi = mod_phase + (pi * i_freq);
        mod_sig = amp * (1/2 + 1/2 * sin(2 * pi * mod_freq * t + mod_phi)); 
        sig(i_freq,:,i_trial) = car_sig .* mod_sig .* stim_on;
    end
end
% plot(t, sig(:,:,1))

%{
% Simulate a dipole
cfg = [];
%cfg.vol = vol; %%% Is this necessary?
cfg.grad = grad.grad;
cfg.dip.pos = [0 0 4];    % you can vary the location, here the dipole is along the z-axis
cfg.dip.mom = [1 0 0]';   % the dipole points along the x-axis
cfg.dip.signal = sig;
cfg.fsample = fsample;
cfg.ntrials = 20;
cfg.relnoise = 10;
cfg.randomseed = 1;
data = ft_dipolesimulation(cfg);
%}

% Construct a Fieldtrip object for the simulated data
% This is similar to performing a RESS/GEDb spatial filter
left_freq_inx = randsample(2, n_trials, true); % Which freq on the left?
left_freq = exp_params.tagged_freqs(left_freq_inx);
hit = randsample(2, n_trials, true) - 1;
trial_num = 225:560;
trialinfo = [hit trial_num' left_freq'];
data = [];
data.label = {'left' 'right'};
data.trialinfo = trialinfo; 
data.trial = cell([1 n_trials]);
data.time = cell([1 n_trials]);
for i_trial = 1:n_trials
    data.time{i_trial} = t;
    if left_freq(i_trial) == 63
        s = sig(:,:,i_trial);
    else
        s = flipud(sig(:,:,i_trial));
    end
    data.trial{i_trial} = s;
end

