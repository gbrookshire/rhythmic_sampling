function rs_tagged_spect(i_subject, segment_type)

% Does power at the RFT frequency vary rhythmically?

% Read in the pre-saved TFR
% Extract power at the tagged frequencies
% Power spectrum of RFT power
%   FFT on each trial
%   Normalize by dividing the area under the curve
%   Then average trial-wise FFTs
% To avoid the stimulus-driven response
%   Start at 0.5 s?
%   Subtract out the average response across trials?

% i_subject = 1;
% segment_type = 'trial';

rs_setup
approx_eq = @(x,y) abs(x - y) < 0.1;
fname = subject_info.meg{i_subject};

% Load the TFR
d = load([exp_dir 'tfr/' segment_type '/' fname '/high']);
d = d.high_freq_data;

% Average over channels in the ROI
cfg = [];
cfg.channel = snr_roi;
cfg.avgoverchan = 'yes';
cfg.latency = [0.5 4];
d = ft_selectdata(cfg, d);

% Make into a data structure for getting FFTs
d_new = [];
d_new.fsample = 1 / mean(diff(d.time));
d_new.trialinfo = d.trialinfo;
d_new.grad = d.grad;
d_new.label = cellfun(@(x) num2str(x), num2cell(d.freq)'); % 'chan' = freq
d_new.trial = cell([1 size(d.powspctrm, 1)]);
d_new.time = cell([1 size(d.powspctrm, 1)]);
for i_rpt = 1:size(d.powspctrm, 1)
    curr_rpt = squeeze(d.powspctrm(i_rpt,1,:,:)); % Data from current rpt
    active_samples = all(~isnan(curr_rpt), 1); % CHECK: Right dimension?
    d_new.time{i_rpt} = d.time(active_samples);;
    d_new.trial{i_rpt} = curr_rpt(:, active_samples);
end
d = d_new;
clear d_new

% Split into separate segments
cfg = [];
cfg.minlength = 1; % Toss segments smaller than 1 s
cfg.length = 1; % Split into n-second segments
cfg.overlap = 0.8; % Segments overlap by this prop
d = ft_redefinetrial(cfg, d);

% Find trials that overlap with the response, or occur after the response
behav = rs_behavior(i_subject);
includes_resp = nan(size(d.time));
for i_rpt = 1:length(d.time)
    n_trial = d.trialinfo(i_rpt, 2);
    resp_time = behav.rt(behav.TrialNumber == n_trial);
    t = d.time{i_rpt};
    includes_resp(i_rpt) = (min(t) < resp_time);
end

% Compute the spectra
cfg = [];
cfg.trials = ~includes_resp; % Exclude segments that overlap with a resp
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'hanning';
cfg.polyremoval = 1; % Remove linear trends
spectra = ft_freqanalysis(cfg, d);

save([exp_dir 'tfr/' segment_type '/' '/' fname '/spect'], 'spectra')

% plot(f, squeeze(mean(spectra, 1)))
% xlabel('Frequency (Hz)')
% ylabel('Power')
