function rs_lga(i_subject)

% Reproduce (something similar to) the analysis from Landau et al 2015

% Load the power at each tagged frequency
% Lateralized Gamma Amplitude
%   Difference b/w left and right power
%   Source-projected gamma-band activity contralateral to the target minus
%   the homologous ipsilateral activity was referred to as lateralized
%   gamma-band activity (LGA).

rs_setup

% Load target-segmented HF TFR
% Read in the data segmented around targets
fname = subject_info.meg{i_subject};
fn = [exp_dir 'tfr/target/' fname '/high'];
d = load(fn);
d = d.high_freq_data;

% Load the behavioral data
behav = rs_behavior(i_subject);

% Compute LGA
d_all = {};
sides = {'left' 'right'};
for side_63Hz = sides
    
    % Select trials with 63 Hz oscillation on this side
    if strcmp(side_63Hz, 'left')
        trial_freq_sel = d.trialinfo(:,3) == 63;
    elseif strcmp(side_63Hz, 'right')
        trial_freq_sel = d.trialinfo(:,3) == 78;
    end
    
    for side_target = sides
        
        % Select trials with the target on this side
        trial_num = d.trialinfo(:, 2);
        targ_side_sel = strcmp(side_target, behav.target_side(trial_num));
        trial_sel = trial_freq_sel & targ_side_sel;

        % Find which side is ipsi/contralateral to the target
        if strcmp(side_target, 'left')
            side_ipsi = 'left';
            side_contra = 'right';
        elseif strcmp(side_target, 'right')
            side_ipsi = 'right';
            side_contra = 'left';
        else
            error('could not find the ipsi/contra side')
        end

        % Find which frequency is ipsi/contralateral to the target
        if strcmp(side_ipsi, side_63Hz)
            f_ipsi = 63;
            f_contra = 78;
        elseif strcmp(side_contra, side_63Hz)
            f_ipsi = 78;
            f_contra = 63;
        else
            error('could not find freq on each side')
        end
        
        % Get the power at the tagged frequencies
        
        % Ipsi
        cfg = [];
        cfg.trials = trial_sel;
        cfg.channel = side_ipsi;
        cfg.frequency = f_ipsi + [-0.1 0.1];
        d_ipsi = ft_selectdata(cfg, d);
        
        % Contra
        cfg = [];
        cfg.trials = trial_sel;
        cfg.channel = side_contra;
        cfg.frequency = f_contra + [-0.1 0.1];
        d_contra = ft_selectdata(cfg, d);
        
        % Compute LGA
        d_lga = d_ipsi;
        d_lga.label = {'lga'};
        d_lga.freq = 0;
        % TODO
        % Does plain subtraction cut it, since the power is so different?
        d_lga.powspctrm = d_contra.powspctrm - d_ipsi.powspctrm;
        
        % Select samples immediately before the target
        % TODO: Redo analyses with shorter TFR window
        tfr_window_length = max(d.cfg.t_ftimwin);
        cfg = [];
        cfg.latency = [-0.5 -tfr_window_length];
        d_lga = ft_selectdata(cfg, d_lga);
        
        d_all{end+1} = d_lga;
    end
end

% Put all the data together
cfg = [];
d_lga = ft_appenddata(cfg, d_all{:}); % TEST: d_all{1}, d_all{3}
d_lga.label = {'lga'};

% Toss all incomplete trials
trial_length = cellfun(@length, d_lga.time);
cfg = [];
cfg.trials = trial_length == max(trial_length);
d_lga = ft_selectdata(cfg, d_lga);


% Run FFTs on LGA for each trial
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper = 'hanning';
cfg.keeptrials = 'yes';
cfg.pad = 'nextpow2';
d_spect = ft_freqanalysis(cfg, d_lga);

% Phase Consistency Metric - PCM
% Compare cosine of the phase difference for
% every possible pairing of hits/misses
hit_inx = find(d_spect.trialinfo(:,1) == 1);
miss_inx = find(d_spect.trialinfo(:,1) == 0);
cosine = nan([length(d_spect.freq), length(hit_inx) * length(miss_inx)]);
k = 1;
for i_hit = hit_inx'
    for i_miss = miss_inx'
        fourier_hit = d_spect.fourierspctrm(i_hit,:,:);
        fourier_miss = d_spect.fourierspctrm(i_miss,:,:);
        phase_diff = angle(fourier_hit .* conj(fourier_miss)); % In radians
        c = cos(phase_diff);
        cosine(:,k) = c;
        k = k + 1;
    end
end

%% Plot for testing
plot(d_spect.freq, mean(cosine, 2))