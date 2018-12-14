function [y, A] = rs_gedb(data, f_low, f_high)

% Generalized eigendecomposition, broadband
% Adapted from
%   https://github.com/svendaehne/matlab_SPoC/blob/master/SSD/ssd.m

% Append all trials into one big matrix
x = cat(2, data.trial{:});

% Get covariance matrix of the broadband activity
cov_bb = cov(x', 1);

% Get covariance matrix of narrowband-filtered activity
signal_band = [f_low, f_high];
[b,a] = butter(filter_order, signal_band / (data.fsample / 2));
x_filt = filtfilt(b, a, x); % Does this need to be transposed?
cov_filt = cov(x_filt', 1);

% Get the joint eigenvalues of cov_bb' * cov_filt
[V, D] = eig(cov_bb' * cov_filt); % GEDb
%%%% Should be the same as the following...? (from Gulbanaite & Cohen
[V, D] = eig(cov_filt, cov_bb); % RESS

% Sort the eigenvalues/vectors
[~, sort_idx] = sort(diag(D), 'descend');
W = V(:, sort_idx);

% Get the timecourse of the GEDb spatial filter
y = W(:,1) * x;

% Get the matrix of the patterns (in columns)
A = cov_filt * W / (W' * cov_filt * W);

end