function [cfc_data, mod_freq] = cfc2(data, freq, nfft, width)

% Perform phase-to-power coupling using coherence
%
% Append all data first to improve robustness of CFC calculation.
%
% ARGS
% data: Fieldtrip data structure
% freq: Vector of frequencies at which to look for modulation
% nfft: size of the FFT window
% width: How many cycles to include in the wavelet analysis (vec or scalar)
%
% OUTPUTS
% cfc_data: Coherence values
% mod_freq: The frequencies of modulation for the coherence

fsample = data.fsample;
width = width * ones(size(freq));

% Set up the wavelets
wavelets = cell([1 length(freq)]);
for i_freq = 1:length(freq)
    f = freq(i_freq);
    N = floor(width(i_freq) * fsample / f);
    taper = hanning(N)';
    wavelet = taper .* exp(1i * 2 * pi * f .* (1:N) / fsample);
    wavelets{i_freq} = wavelet;
end
clear wavelet

% Compute CFC
cfc_data = nan(length(freq), ... % CarFreq x ModFreq x Chan
    nfft, ...
    length(data.label));
% 
for i_channel = 1:length(data.label)
    for i_freq = 1:length(freq)
        % Each run through this takes ~1.6 sec
        % 1.6 s * 304 (channels) * 61 (freqs) = 8.25 hours

        % Compute the power time-course
        % This loop takes ~0.2 s
        pwr = cell([1 length(data.trial)]);
        for i_trial = 1:length(data.trial)
            s = data.trial{i_trial}(i_channel, :);
            c = conv(s, wavelets{i_freq}, 'same');
            sP = abs(c) .^ 2; % Fluctuations in power over time
            pwr{i_trial} = sP;
        end
                
        % Append all the trials into one big trial separated by NaNs
        x_padded = cell([length(data.trial) 1]); % Raw data
        y_padded = cell([length(data.trial) 1]); % Power time-course
        % This loop takes ~1.3 s
        for i_trial = 1:length(data.trial)
            s = data.trial{i_trial}(i_channel, :); % Raw data
            % Detrend the segment
            s = detrend(s);
            % Pad out to the next multiple of NFFT
            trial_length = size(data.trial{i_trial}, 2);
            pad_length = (ceil(trial_length / nfft) * nfft) - trial_length;
            x_padded{i_trial} = nanpad(s, pad_length);
            y_padded{i_trial} = nanpad(pwr{i_trial}, pad_length);
        end
        x = cat(2, x_padded{:}); % Raw data
        y = cat(2, y_padded{:}); % Power time-course
        clear x_padded y_padded

        % Compute the CFC
        
        % Split each signal into overlapping segments
        seg_length = nfft;
        seg_overlap = 0.5;
        seg_overlap = round(nfft * seg_overlap);
        buf = @(z) buffer(z, seg_length, seg_overlap, 'nodelay');
        x_s = buf(x);
        y_s = buf(y);

        % Apply a hanning taper to each segment
        taper = @(z) bsxfun(@times, z, hann(nfft));
        x_s = taper(x_s);
        y_s = taper(y_s);

        % Compute FFT of each segment
        X_s = fft(x_s, [], 1);
        Y_s = fft(y_s, [], 1);

        % Compute the cross-spectra
        xspec = X_s .* conj(Y_s);

        % Compute cross-frequency coupling
        nansum = @(a, dim) sum(a, dim, 'omitnan');
        num = abs(nansum(xspec, 2));
        denom = sqrt(nansum(abs(X_s).^2, 2) .* nansum(abs(Y_s).^2, 2));
        cfc = num ./ denom;
        cfc_data(i_freq, :, i_channel) = cfc;

    end
end

% Only keep the real mod freqs
n_keep_freqs = floor(size(cfc_data, 2) / 2);
cfc_data = cfc_data(:, 1:n_keep_freqs, :); 

mod_freq = (0:(nfft - 1)) * fsample / nfft;
mod_freq = mod_freq(1:n_keep_freqs);
end

function b = nanpad(a, n)
% Pad an array <a> with <n> nans, distributed across the front and back
nrows = size(a, 1);
pre = nan(nrows, floor(n / 2));
post = nan(nrows, ceil(n / 2));
b = cat(2, pre, a, post);
end