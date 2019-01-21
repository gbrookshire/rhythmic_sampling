function rs_tagged_spect(i_subject)

% Does power at the RFT frequency vary rhythmically?

% Read in the pre-saved TFR
% Extract power at the tagged frequencies
% Power spectrum of RFT power
%   FFT on each trial
%   Normalize by dividing the area under the curve
%   Then average trial-wise FFTs
% To avoid the stimulus-driven response
%   Subtract out the average response across trials?

rs_setup
fname = subject_info.meg{i_subject};

% Load the data
behav = rs_behavior(i_subject);
d = load([exp_dir 'tfr/trial/' fname '/high']);
d = d.high_freq_data;

%{
%%% TESTING WITH WHITE NOISE
nan_inx = isnan(d.powspctrm);
d.powspctrm = rand(size(d.powspctrm));
d.powspctrm(nan_inx) = nan;
%}

% To facilitate FFT calculation, make a FT struct with a separate 'channel'
% for each frequency at each RESS filter
data_ress = [];
for side = {'left' 'right'}
    % Make a FT data struc of RESS filter for one side
    chan_inx = strcmp(side{1}, d.label);
    d_side = [];
    d_side.trialinfo = d.trialinfo;
    % Set channel names as the frequencies
    d_side.label = cellfun(@(x) sprintf('%.6f', x), ...
        num2cell(d.freq)', ...
        'UniformOutput', false);
    d_side.trial = cell([1 size(d.powspctrm, 1)]);
    d_side.time = cell([1 size(d.powspctrm, 1)]);
    for i_rpt = 1:size(d.powspctrm, 1)
        curr_rpt = squeeze(d.powspctrm(i_rpt,chan_inx,:,:)); % Current rpt
        active_samples = all(~isnan(curr_rpt), 1); % CHECK: Right dimension?
        d_side.time{i_rpt} = d.time(active_samples);;
        d_side.trial{i_rpt} = curr_rpt(:, active_samples);
    end
    fprintf('\nTRIALS: %i\n', length(d_side.trial))
    
    % Toss trials with no samples
    % This happens when calculating TFR on short trials
    cfg = [];
    cfg.trials = cellfun('length', d_side.time) > 1;
    d_side = ft_selectdata(cfg, d_side);
    fprintf('\nTRIALS: %i\n', length(d_side.trial))

    % Split into short segments
    cfg = [];
    cfg.minlength = 1; % Toss segments smaller than 1 s
    cfg.length = 1; % Split into n-second segments
    cfg.overlap = 0.8; % Segments overlap by this prop
    d_side = ft_redefinetrial(cfg, d_side);
    fprintf('\nTRIALS: %i\n', length(d_side.trial))

    % Toss segments that overlap with or occur after the response
    % Or that include the transient response at the beginning of the trial
    includes_resp = nan(size(d_side.time));
    beginning_of_trial = nan(size(d_side.time));
    for i_rpt = 1:length(d_side.time)
        n_trial = d_side.trialinfo(i_rpt, 2);
        resp_time = behav.rt(behav.TrialNumber == n_trial);
        t = d_side.time{i_rpt};
        includes_resp(i_rpt) = resp_time < max(t);
        beginning_of_trial(i_rpt) = min(t) < 0.5;
    end
    cfg = [];
    cfg.trials = ~includes_resp & ~beginning_of_trial;
    d_side = ft_selectdata(cfg, d_side);
    fprintf('\nTRIALS: %i\n', length(d_side.trial))

    % Compute the spectra
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.taper = 'hanning';
    cfg.polyremoval = 1; % Remove linear trends
    cfg.keeptrials = 'yes';
    spectra = ft_freqanalysis(cfg, d_side);
    
    data_ress.(side{1}) = spectra;
    clear d_side spectra;
end    

save([exp_dir 'tfr/trial/' fname '/spect'], 'data_ress')
