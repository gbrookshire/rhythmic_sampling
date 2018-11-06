function tagged = rs_tagged(data)

% Return a Fieldtrip structure with power at the two tagged frequencies

rs_setup

%{
% With the Hilbert transform
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [-1 1] + exp_params.tagged_freqs;
cfg.hilbert = 'abs'; % Get the power at tagged freq
cfg.padding = 15; % Pad the data to reduce filtering artifacts
cfg.padtype = 'data';
tagged = ft_preprocessing(cfg, dat);
%}

% %{
% With a standard TFR procedure
% Think about the best way to choose these parameters
window_length = 0.1; % Length of FFT window
window_interval = 0.01; % How much time between centers of FFT windows
cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = exp_params.tagged_freqs;
cfg.t_ftimwin = ones(size(cfg.foi)) .* window_length;
cfg.toi = -0.5:window_interval:4; % Will this add NaNs for short trials?
cfg.keeptrials = 'yes';
tagged = ft_freqanalysis(cfg, data);
%}