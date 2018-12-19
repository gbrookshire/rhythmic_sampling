function [y, A] = rs_gedb(data, f_low, f_high)

% Generalized eigendecomposition, broadband
% Adapted from
%   https://github.com/svendaehne/matlab_SPoC/blob/master/SSD/ssd.m
% --- This script mixes together RESS and GEDb. Switch to just RESS

% Append all trials into one big matrix
x = cat(2, data.trial{:});

% Get covariance matrix of the broadband activity
cov_bb = cov(x', 1);

% Get covariance matrix of narrowband-filtered activity
signal_band = [f_low, f_high];
[b,a] = butter(filter_order, signal_band / (data.fsample / 2));
x_filt = filtfilt(b, a, x); % Does this need to be transposed?
cov_filt = cov(x_filt', 1);

% Run the generalized eigendecomposition
[V, D] = eig(inv(cov_bb) * cov_filt); % GEDb
[V, D] = eig(cov_filt, cov_bb); % RESS

% % Two equivalent ways to compute generalized eigendecomposition:
% % Derivation:
% % [V,D] = eig(A, B)
% % A*V = B*V*D
% % inv(B)*A*V = V*D
% % Which is equivalent to normal eigendecomposition if inv(B)*A is treated
% % as a single matrix C: C*V = V*D
% % [V,D] = eig(inv(B), A)
% % Example:
% A = rand(10,2);
% B = rand(10,2);
% A = cov(A);
% B = cov(B);
% [V1,D1] = eig(A, B); % First way
% [V2,D2] = eig(inv(B) * A); % Second way
% rounding_error = 1e-10;
% all(D1 - D2 < rounding_error) % Eigenvalues are almost identical
% V1 ./ V2 % Eigenvectors are the same after multiplying by a scalar
% % The first method (eig(A, B)) is better for 'numerical stability'
% according to Gulbinaite and Cohen (2017).

% Sort the eigenvalues/vectors
[~, sort_idx] = sort(diag(D), 'descend');
W = V(:, sort_idx);

% Get the timecourse of the spatial filter
y = W(:,1) * x;

% Get the matrix of the patterns (in columns)
A = cov_filt * W / (W' * cov_filt * W);

end