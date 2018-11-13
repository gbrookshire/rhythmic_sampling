function data_out = rs_snr(data_in, freq)

% Compute the SNR on a Fieldtrip data structure
% data: Fieldtrip data structure
% freq: Frequency to compute the SNR for

% Get power at the tagged freq
cfg = [];
cfg.frequency = [-0.5 0.5] + freq;
cfg.avgoverfreq = 'yes';
p_center = ft_selectdata(cfg, data_in);

% Get power at a lower frequency
cfg.frequency = [-1 1] + freq - 2;
p_low = ft_selectdata(cfg, data_in);

% Get power at a higher frequency
cfg.frequency = [-1 1] + freq + 2;
p_high = ft_selectdata(cfg, data_in);

% Get the field that holds freq data - pow or powspctrm
if isfield(data_in, 'pow')
    param = 'pow';
elseif isfield(data_in, 'powspctrm')
    param = 'powspctrm';
else
    error('No ''pow'' or ''powspctrm'' field found in the data')
end

% Calculate the SNR
s = p_center.(param) ./ ((p_low.(param) + p_high.(param)) / 2);

% Put it together into an output object
data_out = p_center;
data_out.(param) = s;