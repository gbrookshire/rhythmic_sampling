function out = rs_simulate_alpha_peaks()

% Simulate the signal with alpha and noise
% Alpha is BP-filtered noise (6-14 Hz)
% Predicted signal, RFTs modulated by alpha (counterphase)


% Basic params and setup
fsample = 1000;
t = -0.5:(1/fsample):4.5;
n_trials = 336;
trials = cell(1, n_trials);

for i_trial = 1:n_trials
    % Make alpha
    x_noise = rand(size(t));
    f_low = 6;
    f_high = 14;
    [b,a] = butter(5, [f_low, f_high] / (fsample / 2), 'bandpass');
    x_filt = filter(b, a, x_noise);
    x_filt = x_filt - min(x_filt);
    x_filt = x_filt / max(x_filt);

%     plot(t, x_noise, '-k')
%     hold on
%     plot(t, x_filt, '-r')
%     hold off

    % Make RFT signals
    f1 = 63;
    f2 = 78;
    rft = @(f) (1/2) * (1 + sin(t * f * 2 * pi));
    rft1 = rft(f1);
    rft2 = rft(f2);

    % Combine them
    % Square alpha power to limit the duty-cycle of RFT oscillations
    x_comb1 = rft1 .* (x_filt .^ 2) + x_filt;
    x_comb2 = rft2 .* ((1 - x_filt) .^ 2) + x_filt;

%     plot(t, x_comb1, t, x_comb2);
%     xlim([0 0.5])

    % Add noise
    noise_amp = 0.25;
    y1 = x_comb1 + (noise_amp * rand(size(t)));
    y2 = x_comb2 + (noise_amp * rand(size(t)));

%     plot(t, y1, t, y2)
%     xlim([0 0.5])

    % Make into three channels with a different mix of each signal
    y = [y1; y1 + y2; y2];
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