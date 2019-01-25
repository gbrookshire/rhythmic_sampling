function [data_out, maps, ress_filt] = rs_ress(data_in, f, fwhm)

% Spatial filter for extracting periodic activity
% RESS (Gulbinaite & Cohen, 2017 NeuroImage)
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
%
% g.brookshire@bham.ac.uk
% Adapted from mikexcohen@gmail.com (and including his code for the
% gaussian filter)

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

data_out = rs_applyressfilt(data_in, ress_filt);

end


function [filtdat,empVals] = filterFGx(data,srate,f,fwhm,showplot)
% filterFGx   Narrow-band filter via frequency-domain Gaussian
%  [filtdat,empVals] = filterFGx(data,srate,f,fwhm,showplot)
% 
% 
%    INPUTS
%       data : 1 X time or chans X time
%      srate : sampling rate in Hz
%          f : peak frequency of filter
%       fhwm : standard deviation of filter, 
%              defined as full-width at half-maximum in Hz
%   showplot : set to true to show the frequency-domain filter shape
% 
%    OUTPUTS
%    filtdat : filtered data
%    empVals : the empirical frequency and FWHM
% 
% Empirical frequency and FWHM depend on the sampling rate and the
% number of time points, and may thus be slightly different from
% the requested values.
% 
% mikexcohen@gmail.com

%% input check

if size(data,1)>size(data,2)
    help filterFGx
    error('Check data size')
end

if (f-fwhm)<0
%     help filterFGx
%     error('increase frequency or decrease FWHM')
end

if nargin<4
    help filterFGx
    error('Not enough inputs')
end

if fwhm<=0
    error('FWHM must be greater than 0')
end

if nargin<5
    showplot=false;
end

%% compute and apply filter

% frequencies
hz = linspace(0,srate,size(data,2));

% create Gaussian
s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-f;                 % shifted frequencies
fx = exp(-.5*(x/s).^2);    % gaussian
fx = fx./max(fx);          % gain-normalized

%% filter

filtdat = 2*real( ifft( bsxfun(@times,fft(data,[],2),fx) ,[],2) );

%% compute empirical frequency and standard deviation

idx = dsearchn(hz',f);
empVals(1) = hz(idx);

% find values closest to .5 after MINUS before the peak
empVals(2) = hz(idx-1+dsearchn(fx(idx:end)',.5)) - hz(dsearchn(fx(1:idx)',.5));

%% inspect the Gaussian (turned off by default)

if showplot
    figure(10001+showplot),clf
    plot(hz,fx,'o-')
    hold on
    plot([hz(dsearchn(fx(1:idx)',.5)) hz(idx-1+dsearchn(fx(idx:end)',.5))],[fx(dsearchn(fx(1:idx)',.5)) fx(idx-1+dsearchn(fx(idx:end)',.5))],'k--')
    set(gca,'xlim',[max(f-10,0) f+10]);
    
    title([ 'Requested: ' num2str(f) ', ' num2str(fwhm) ' Hz; Empirical: ' num2str(empVals(1)) ', ' num2str(empVals(2)) ' Hz' ])
    xlabel('Frequency (Hz)'), ylabel('Amplitude gain')
end

%% done.
end