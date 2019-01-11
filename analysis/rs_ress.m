function [data_out, maps] = rs_ress(data_in, f, fwhm)

% RESS (Gulbinaite & Cohen, 2017 NeuroImage)
% Spatial filter for extracting periodic activity
%
% INPUTS
%   data_in: A fieldtrip data structure (from preprocessing)
%   f: peak frequency of filter
%   fhwm: standard deviation of filter, 
%         defined as full-width at half-maximum in Hz
%
% OUTPUTS
%   data_out: A fieldtrip data structure of activity at the RESS filter
%   maps: Scalp topography of the filters

% Append all trials into one big matrix (Channel * Time)
x = cat(2, data_in.trial{:});

% Check whether the input data have full rank. If not, eigenvalues will not
% be stable. This will happen if components are rejected from ICA during
% preprocessing.
if rank(x') ~= size(x,1)
    warning('Input data are rank-deficient; eigenvalues may be unstable')
end

% Mean-center each row before getting covariance matrix
% This only affects the covariance matrices very weakly (1e-38)
center = @(x) bsxfun(@minus, x, mean(x, 2));

% Get covariance matrix of the broadband activity
cov_bb = cov(center(x)', 1);

% Get covariance matrix of narrowband-filtered activity
fsample = 1 / mean(diff(data_in.time{1}));
x_filt = filterFGx(x, fsample, f, fwhm);
cov_filt = cov(center(x_filt)', 1);

% perform generalized eigendecomposition
[evecs, evals] = eig(cov_filt, cov_bb);

if ~isreal(evecs)
    warning('GED found complex-valued eigenvectors')
end
% find maximum component
[~, comp2plot] = max(diag(evals));
% normalize vectors (not really necessary, but OK)
evecs = bsxfun(@rdivide, evecs, sqrt(sum(evecs .^ 2, 1)));

% extract components and force sign
% maps = inv(evecs'); % get maps (this is fine for full-rank matrices)
% Get maps in a way that's robust to rank-deficiency
maps = cov_filt * evecs / (evecs' * cov_filt * evecs);
[~,idx] = max(abs(maps(:, comp2plot))); % find biggest component
maps = maps * sign(maps(idx, comp2plot)); % force to positive sign

% % reconstruct RESS component time series
% ress_ts1 = zeros(EEG.pnts, size(data,3));
% for ti=1:size(data,3)
%     ress_ts1(:,ti) = evecs(:,comp2plot)' * squeeze(data(:,:,ti));
% end

% Make Fieldtrip-style object to return out
data_out = [];
data_out.fsample = fsample;
data_out.label = {'RESS'};
data_out.trialinfo = data_in.trialinfo;
data_out.time = data_in.time;
data_out.trial = cell(size(data_in.trial));
for i_trial = 1:length(data_in.trial)
    trial_data = data_in.trial{i_trial};
    data_out.trial{i_trial} = evecs(:, comp2plot)' * squeeze(trial_data);
end

end

%% Testing
%{

% Try looking at neighboring frequencies instead of BB activity

rs_setup
i_subject = 1;
fname = subject_info.meg{i_subject};
data_preproc = rs_preproc(i_subject, 'trial');
grad = load([exp_dir 'grad/' fname '/grad'], 'grad');
behav = rs_behavior(i_subject);

filter_freq = exp_params.tagged_freqs(2);

% Select trials with consistent freq/side mapping
keep_trials = ismember(...
    data_preproc.trialinfo(:,2), ...
    find(behav.freq_right == filter_freq));
cfg = [];
cfg.trials = keep_trials;
data_sub = ft_selectdata(cfg, data_preproc);

% Compute RESS components
[data_ress, maps] = rs_ress(data_sub, filter_freq, 0.5);

% Plot the RESS component maps
% Make aggregated data structure
n_maps = size(maps, 2);
d_maps = [];
d_maps.label = data_sub.label;
d_maps.time = 1:n_maps; % Actually not time, but component number
d_maps.avg = real(maps);
d_maps.dimord = 'chan_time';
d_maps.grad = grad.grad;
% Combine planar gradiometers
cfg = [];
cfg.method = 'sum';
d_maps = ft_combineplanar(cfg, d_maps);
% Plot the maps
figure
for i_map = 1:6
    subplot(4, 3, i_map)
    cfg = [];
    cfg.marker = 'off';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.layout = chan.grad_cmb.layout;
    cfg.xlim = [-0.1 0.1] + i_map;
    ft_topoplotER(cfg, d_maps)
end
subplot(4,3,1)
title('Retained comp')

% Check the spectra of the time-course
% Compute the spectra
cfg = [];
cfg.channel = {'RESS'};
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'hanning';
cfg.pad = 'nextpow2';
cfg.padtype = 'zero';
cfg.polyremoval = 1; % Remove linear trends
spec = ft_freqanalysis(cfg, data_ress);
subplot(2,1,2)
pow = 20 * log10(spec.powspctrm);
plot(spec.freq, pow)
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xlim([0 100])
%ylim([-580 -550])
hold on
plot(filter_freq, max(pow), 'vr')
hold off
%}