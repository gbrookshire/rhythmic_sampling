function out = rs_simulate_alpha_peaks()

% Simulate the signal with alpha and noise
% Alpha is BP-filtered noise (6-14 Hz)
% Predicted signal, RFTs modulated by alpha (counterphase)


% Basic params and setup
fsample = 1000;
t = -0.5:(1/fsample):4.5;
n_trials = 336;
trials = cell(1, n_trials);
alpha_amp = 1;
rft_amp = 0.75;

for i_trial = 1:n_trials
    % Make alpha
    x_noise = rand(size(t));
    f_low = 6;
    f_high = 14;
    [b,a] = butter(5, [f_low, f_high] / (fsample / 2), 'bandpass');
    x_filt = filter(b, a, x_noise);
    x_filt = x_filt - min(x_filt);
    x_filt = alpha_amp * x_filt / max(x_filt);

    % Make RFT signals
    f1 = 63;
    f2 = 78;
    rft = @(f) rft_amp * (1/2) * (1 + sin(t * f * 2 * pi));
    rft1 = rft(f1);
    rft2 = rft(f2);

    % Combine them
    s = 4; % 'Sharpness' exponent - higher numbers reduce the duty cycle
    c1 = rft1 .* (x_filt .^ s);
    c2 = rft2 .* ((1 - x_filt) .^ s);
    c1 = detrend(c1);
    c2 = detrend(c2);
    c_comb = c1 + c2;

    % Make into three channels with a different mix of each signal
    y = [c1; c_comb; c2;];
    y = y + x_filt;
    
    % Add noise
    noise_amp = 0.0;
    y = y + (noise_amp * rand(size(y)));   
    
    trials{i_trial} = y;
end
    
% Put together a final group of signals in a fieldtrip data structure
labels = {'MEG2031' 'MEG2032' 'MEG2033'};
out = [];
out.label = labels;
out.trial = trials;
out.trialinfo = [(225:560)', round(rand(n_trials, 1))];
out.time = cell(1, n_trials);
for i_trial = 1:n_trials
    out.time{i_trial} = t;
end

end