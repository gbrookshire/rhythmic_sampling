function [cfc_data, mod_freq] = cfc2(data, freq, nfft, width)

% Load data
data = rs_simulate_flicker();
% d.fsample = mean(diff(d.time{1}));
cfg = [];
cfg.channel = {'left' 'right'};
data = ft_selectdata(cfg, data);
% For testing, make each trial a random length from 1.5 - 4.5
for i_trial = 1:length(data.trial)
    len = (4.5 - 1.5) .* rand(1) + 1.5;
    keep_samps = data.time{i_trial} < len;
    data.trial{i_trial} = data.trial{i_trial}(:,keep_samps);
    data.time{i_trial} = data.time{i_trial}(:,keep_samps);
end

% Parameters
freq = [63 78];
width = 6;
nfft = 2^10;

fsample = data.fsample;
width = width * ones(size(freq));

% Compute the power time-course
pwr = cell([1 length(data.trial)]);
for i_trial = 1:length(data.trial)
%     % Within that cell: Chan x Time x Freq matrices
%     p = nan(length(d.label), length(d.time{i_trial}), length(freq));
    p = nan(1, length(data.time{i_trial})); % Only 1 chan/freq at a time
    for i_freq = 1:length(freq)
        f = freq(i_freq);
        N = floor(width(i_freq) * fsample / f);
        taper = hanning(N)';
        wavelet = taper .* exp(1i * 2 * pi * f .* (1:N) / fsample);
        for i_channel = 1:length(data.label)
            s = data.trial{i_trial}(i_channel, :);
            c = conv(s, wavelet, 'same');
            sP = abs(c) .^ 2;
            p(i_channel, :, i_freq) = sP;
        end
    end
    pwr{i_trial} = p;
end

% Append all the data into a matrix
%   Raw data: Chan x Time
%   Power: Chan x Time x Freq
% Pad each trial with NaNs to an integer mult of nfft
% This ensures that FFTs only cover data from a contiguous segment of data
append_nan = @(a,n) cat(2, a, nan(size(a, 1), n, size(a, 3)));
x_padded = cell([length(data.trial) 1]);
y_padded = cell([length(data.trial) 1]);
for i_trial = 1:length(data.trial)
    trial_length = size(data.trial{i_trial}, 2);
    pad_length = (ceil(trial_length / nfft) * nfft) - trial_length;
    x_padded{i_trial} = append_nan(data.trial{i_trial}, pad_length);
    y_padded{i_trial} = append_nan(pwr{i_trial}, pad_length);
end
x = cat(2, x_padded{:}); % Raw data
y = cat(2, y_padded{:}); % Power time-course
clear x_padded y_padded

% Compute CFC on each channel & frequency
cfc_data = nan(length(data.label), ... % Chan x CarFreq x ModFreq
    length(freq), ...
    nfft);
for i_channel = 1:length(data.label)
    for i_freq = 1:length(freq)
        
        % Select data from this chan and freq
        x_sel = x(i_channel, :);
        y_sel = y(i_channel, :, i_freq);
        
        % Split each signal into overlapping segments
        seg_length = nfft;
        seg_overlap = 0.5;
        seg_overlap = round(nfft * seg_overlap);
        buf = @(z) buffer(z, seg_length, seg_overlap, 'nodelay');
        x_s = buf(x_sel);
        y_s = buf(y_sel);

        % Apply a hanning taper to each segment
        taper = @(z) bsxfun(@times, z, hann(nfft));
        x_s = taper(x_s);
        y_s = taper(y_s);

        % Compute FFT of each segment
        X_s = fft(x_s, [], 1);
        Y_s = fft(y_s, [], 1);
        mod_freq = (0:(nfft - 1)) * fsample / nfft;

        % Compute the cross-spectra
        xspec = X_s .* conj(Y_s);

        % Compute cross-frequency coupling
        nansum = @(a, dim) sum(a, dim, 'omitnan');
        num = abs(nansum(xspec, 2));
        denom = sqrt(nansum(abs(X_s).^2, 2) .* nansum(abs(Y_s).^2, 2));
        cfc = num ./ denom;
        cfc_data(i_channel, i_freq, :) = cfc;
    end
end

