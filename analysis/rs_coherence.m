function rs_coherence(i_subject)

% Compute coherence between the MEG recordings and the photodiode
% 
% Get coherence separately for hits-on-the-left, and hits-on-the-right
% Separate by frequency contralateral to the stimulus

rs_setup

fname = subject_info.meg{i_subject};

% Load preprocessed data
behav = rs_behavior(i_subject);
d_meg = rs_preproc_ress(i_subject, 'target');
d_photo = load([exp_dir 'preproc/photodiode/target/' fname '/preproc']);
d_photo = d_photo.data;
data = ft_appenddata([], d_meg, d_photo);
data.trialinfo = d_meg.trialinfo;
clear d_meg d_photo1

% Compute TFR 
cfg = [];
cfg.output = 'fourier'; % FUll Fourier representation needed for coherence
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 55:100;
cfg.toi = -0.5:0.05:0.5;
cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.1;
cfg.keeptrials = 'yes';
fourier = ft_freqanalysis(cfg,data);

% Get rid of trials with nans
% And keep only hits
nan_trials = any(squeeze(isnan(fourier.fourierspctrm(:,1,1,:))), 2);
hits = fourier.trialinfo(:,1) == 1;
cfg = [];
cfg.trials = hits & ~nan_trials;
fourier = ft_selectdata(cfg, fourier);

% Find which side each RFT frequency was on in each trial
trial_num = fourier.trialinfo(:,2);
freq_left = behav.freq_left(trial_num);
target_side = behav.target_side(trial_num);

% Get coherence separately for each mapping of freq to each side
sides = {'left' 'right'};
clear cohspctrm % chancmb * freq * time * RFTfreq * target_side
clear trial_counts
for i_freq = 1:2
    % i_freq == 1: left side at 63 Hz, right side at 78 Hz
    % i_freq == 2: left side at 78 Hz, right side at 63 Hz
    freq_sel = freq_left == exp_params.tagged_freqs(i_freq);
    for i_targ_side = 1:2
        % i_targ_side == 1: target is on the left (left-hits)
        % i_targ_side == 2: target is on the right (right-hits)
        targside_sel = strcmp(target_side, sides{i_targ_side});
        % Combine selection vectors to get trial selection
        trial_sel = freq_sel & targside_sel;
        trial_counts(i_freq, i_targ_side) = sum(trial_sel);
        % Compute coherence
        cfg = [];
        cfg.trials = freq_sel & targside_sel;
        cfg.method = 'coh';
        cfg.channelcmb = {'left' 'MISC004'
                          'right' 'MISC004'};
        coh = ft_connectivityanalysis(cfg, fourier);
        cohspctrm(:,:,:,i_freq,i_targ_side) = coh.cohspctrm;
    end
end

% Get left-hits vs right-hits
x_l = cohspctrm(:,:,:,:,1);
x_r = cohspctrm(:,:,:,:,2);
x = (x_l - x_r) ./ (x_l + x_r);
% Dimensions of x: chancmb * freq * time * RFT_freq_mapping

% Clean out the raw coherence spectrum so we don't get confused
coh = rmfield(coh, 'cohspctrm');

% Save the data
save_dir = [exp_dir 'coherence/'];
[~,~,~] = mkdir(save_dir, fname);
save([save_dir '/' fname '/coh'], 'coh', 'x', 'cohspctrm')
