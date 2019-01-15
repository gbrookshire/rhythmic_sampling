function [data_out, maps, ress_filt] = rs_ress(data_in, f, fwhm)

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
%   ress_filt: The spatial filter to apply to MEEG data

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

% Sort the eigenvalues and -vectors
[srtd, srt_ind] = sort(diag(evals), 'descend');
evals = diag(srtd); % diag matrix of sorted eigenvalues
evecs = evecs(:, srt_ind); % sort the eigenvectors

% normalize vectors (not really necessary, but OK)
evecs = bsxfun(@rdivide, evecs, sqrt(sum(evecs .^ 2, 1)));

% Get maps in a way that's robust to rank-deficiency
% For full-rank matrices, it's fine to do: maps = inv(evecs');
maps = cov_filt * evecs / (evecs' * cov_filt * evecs);

% Make component maps have positive sign
[~,peak_ind] = max(abs(maps), [], 1); % find peak activation in each map
peak_sign = nan(size(peak_ind)); % Init
for i_map = 1:length(peak_ind)
    peak_sign(i_map) = sign(maps(peak_ind(i_map), i_map));
end
maps = maps * diag(peak_sign);

% Matrix to apply the spatial filter to MEEG data
ress_filt = evecs(:,1)';

% % Make Fieldtrip-style object to return out
% data_out = [];
% data_out.fsample = fsample;
% data_out.label = {'RESS'};
% data_out.trialinfo = data_in.trialinfo;
% data_out.time = data_in.time;
% data_out.trial = cell(size(data_in.trial));
% for i_trial = 1:length(data_in.trial)
%     trial_data = data_in.trial{i_trial};
%     % Keep only the top component
%     data_out.trial{i_trial} = ress_filt * squeeze(trial_data);
% end

data_out = rs_applyressfilt(data_in, ress_filt);

end