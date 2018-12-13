function rs_tagged_spect(i_subject, segment_type)

% Does power at the RFT frequency vary rhythmically?

% Read in the pre-saved TFR
% Extract power at the tagged frequencies
% Power spectrum of RFT power
%   FFT on each trial
%   Normalize by dividing the area under the curve
%   Then average trial-wise FFTs
% To avoid the stimulus-driven response
%   Start at 0.5 s?
%   Subtract out the average response across trials?

% TODO:
% - Split into separate segments of equal length
% - Toss segments smaller than that length (or zero-pad?)
% - Toss segments close to the response (0.5 s from the end)
% - Divide by area under the curve

% i_subject = 1;
% segment_type = 'trial';

rs_setup
approx_eq = @(x,y) abs(x - y) < 0.1;
fname = subject_info.meg{i_subject};

% Load the TFR
d = load([exp_dir 'tfr/' segment_type '/' fname '/high']);
d = d.high_freq_data;

% Average over channels in the ROI
cfg = [];
cfg.channel = snr_roi;
cfg.avgoverchan = 'yes';
cfg.latency = [0.5 4]; %FIXME - does this keep post-response stuff?
d = ft_selectdata(cfg, d);

% Compute the spectrum of power at the tagged freq for each trial
fsample = 1 / mean(diff(d.time));
nfft = 2 ^ nextpow2(size(d.powspctrm, 4));
f = fsample / 2 * linspace(0, 1, nfft / 2 + 1); % Frequencies of the FFT
spectra = nan([... % Trial * Time-step * Freq (averaged over channel)
    size(d.powspctrm, 1) ...
    nfft / 2 + 1 ...
    length(d.freq)]);
%     length(exp_params.tagged_freqs)]);
for i_freq = 1:length(d.freq) %length(exp_params.tagged_freqs)
    %freq_inx = approx_eq(exp_params.tagged_freqs(i_freq), d.freq);
    freq_inx = i_freq;
    for i_trial = 1:size(d.powspctrm, 1)
        % Get rid of NaNs in the TFR (from trials shorter than the max time
        x = squeeze(d.powspctrm(i_trial, :, freq_inx, :));
        x(isnan(x)) = [];
        y = fft(x, nfft); % Compute FFT for each trial
        y = y / nfft; % Normalize the amplitude
        y_amp = 2 * abs(y(1:nfft/2+1)); % Single-sided amplitude spect
        spectra(i_trial,:,i_freq) = y_amp; % Store the amplitude spectra
    end
end
warning('Make sure that we dont get power at the trial length!')

plot(f, squeeze(mean(spectra, 1)))
xlabel('Frequency (Hz)')
ylabel('Power')
