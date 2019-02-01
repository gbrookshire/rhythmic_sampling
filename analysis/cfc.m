function [cfc_data, mod_freq] = cfc(data, freq, nfft, width)

% Perform phase-to-power coupling using coherence
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
%
% This will only work if there are no NaNs in the data.
% Use ft_rejectartifacts with 'partial' before this point to avoid NaNs.

% Author: Geoff Brookshire
% Adapted from a script provided by Ole Jensen

% TODO
% This is pretty slow - find a way to speed it up or parallelize it

if length(width) == 1
    width = ones(size(freq)) * width;
elseif length(width) ~= length(freq)
    error('<width> must have one element per frequency')
end

if isfield(data, 'fsample')
    fsample = data.fsample;
else
    fsample = 1 / mean(diff(data.time{1}));
end 

% Initialize data structure for holding coherence
% CarrierFreq x ModFreq x Channel x Trial
cfc_data = nan(length(freq), ...
    nfft / 2 + 1, ...
    length(data.label), ...
    length(data.trial));

for i_freq = 1:length(freq) 
    % Set up the wavelet to get power over time
    f = freq(i_freq);
    N = floor(width(i_freq) * fsample / f);
    taper = hanning(N)';
    wavelet = taper .* exp(1i * 2 * pi * f .* (1:N) / fsample);
    % Compute coherence
    for i_trial = 1:length(data.trial) 
        for i_channel = 1:length(data.label)
            s = data.trial{i_trial}(i_channel, :);
            sP = abs(conv(s, wavelet)) .^ 2;
            sP = sP(ceil(N/2):length(sP)-floor(N/2));
            [coh, mod_freq] = mscohere(s ,sP, ...
                hanning(nfft), nfft / 2, nfft, fsample);
            cfc_data(i_freq, :, i_channel, i_trial) = coh;
        end 
    end 
end

% Only keep the real mod freqs
cfc_data = cfc_data(:, 1:floor(end/2), :, :); 
mod_freq = mod_freq(1:floor(end/2));

% Weighted average over trials
for i_trial = 1:length(data.trial)
    n = length(data.time{i_trial});
    cfc_data(:,:,:,i_trial) = cfc_data(:,:,:,i_trial) .* n;
end
total_length = sum(cellfun(@length, data.time));
cfc_data = sum(cfc_data, 4) ./ total_length;
