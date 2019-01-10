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

% % Setup script for testing
% cd('C:\Users\brookshg\Documents\rhythmic_sampling\rhythmic_sampling\analysis')
% addpath('C:\Users\brookshg\Documents\fieldtrip-20180805');
% load('C:\Users\brookshg\Documents\rhythmic_sampling\sample_data\preproc\trial\181009_b46d\181009\preproc.mat')

% Append all trials into one big matrix (Channel * Time)
x = cat(2, data_in.trial{:});

if rank(x') ~= size(x,1)
    warning('The input data are rank-deficient')
end
% The full data are rank-deficient, but I'm not sure why. When I select
% only the magnetometers, or only the gradiometers, the data have full
% rank. Why does this happen?

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
data_out = data_in;
data_out.label = {'RESS'};
data_out.trial = cell(size(data_in.trial));
for i_trial = 1:length(data_in.trial)
    trial_data = data_in.trial{i_trial};
    data_out.trial{i_trial} = evecs(:, comp2plot)' * squeeze(trial_data);
end

end

%% Testing

% Plot the maps -- right now they make no sense

n_maps = size(maps,1);
d = [];
d.fsample = data_out.fsample;
d.label = data_in.label;
d.time = 1:n_maps; % Actually not time, but component number
d.avg = maps;
d.dimord = 'chan_time';

for i_map = 1:9
    subplot(3, 3, i_map)
    cfg = [];
    cfg.marker = 'off';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.layout = chan.grad.layout;
    cfg.xlim = [-0.1 0.1] + i_map;
    ft_topoplotER(cfg, d)
end

% Check the spectra of the time-course