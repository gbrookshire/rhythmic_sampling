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

% Mean-center each row before getting covariance matrix
center = @(x) bsxfun(@minus, x, mean(x, 2));

% Get covariance matrix of the broadband activity
cov_bb = cov(center(x)', 1);

% Get covariance matrix of narrowband-filtered activity
% signal_band = [f_low, f_high];
% [b,a] = butter(filter_order, signal_band / (data.fsample / 2));
% x_filt = filtfilt(b, a, x); % Does this need to be transposed?
x_filt = filterFGx(x, data_in.fsample, f, fwhm);
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

% reconstruct RESS component time series
ress_ts1 = zeros(EEG.pnts, size(data,3));
for ti=1:size(data,3)
    ress_ts1(:,ti) = evecs(:,comp2plot)' * squeeze(data(:,:,ti));
end

% Make Fieldtrip-style object to return out
data_out = data_in;
data_out.label = {'RESS'};
data_out.trial = cell(size(data_in.trial));
for i_trial = 1:length(data_in.trial)
    trial_data = data_in.trial{i_trial};
    data_out.trial{i_trial} = evecs(:, comp2plot)' * squeeze(trial_data);
end

end