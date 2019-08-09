function rs_coherence(i_subject)

% Compute coherence between the MEG recordings and the photodiode
% 
% Get coherence separately for hits-on-the-left, and hits-on-the-right

rs_setup

fname = subject_info.meg{i_subject};

% Load preprocessed data
d_meg = rs_preproc_ress(i_subject, 'target');
d_photo = load([exp_dir 'preproc/photodiode/target/' fname '/preproc']);
d_photo = d_photo.data;
data = ft_appenddata([], d_meg, d_photo);
data.trialinfo = d_meg.trialinfo;
clear d_meg d_photo1

% Compute TFR 
cfg = [];
cfg.output = 'fourier'; % FUll Fourier representation needed for coherence
% cfg.channel = {'MEGGRAD','MISC004'}; % Gradiometers and the photodiode
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 55:100;
cfg.toi = -0.5:0.05:0.5;
cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.1;
cfg.keeptrials = 'yes';
fourier = ft_freqanalysis(cfg,data);
% There are a small number of nans here
nan_trials = any(squeeze(isnan(fourier.fourierspctrm(:,1,1,:))), 2);

% Compute coherence
cfg = [];
cfg.trials = ~nan_trials;
cfg.method = 'coh';
cfg.channelcmb = {'left' 'MISC004'
                  'right' 'MISC004'};
coh = ft_connectivityanalysis(cfg, fourier);

% % Combine planar gradiometers
% cfg = [];
% cfg.method = 'sum';
% coh.powspctrm = coh.cohspctrm;
% coh.label = fourier.label(1:length(coh.labelcmb));
% coh_cmb = ft_combineplanar(cfg, coh); %cfg.method = 'sum' only works for frequency data with powspctrm
% coh_cmb = rmfield(coh_cmb,'cohspctrm');
% Save the data

save_dir = [exp_dir 'coherence/'];
[~,~,~] = mkdir(save_dir, fname);
save([save_dir '/' fname '/coh'], 'coh', '-v7.3')






