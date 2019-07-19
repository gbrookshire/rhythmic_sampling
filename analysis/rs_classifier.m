function rs_classifier(i_subject, data_type)

% Try to classify Left-Hits vs Right-Hits in each sample
% Run on
%   - Raw MEG
%   - Alpha power
%   - RFT power
%
% data_type: alpha, raw, rft

tic
rs_setup

cfg_sel = []; % for ft_selectdata later on

behav = rs_behavior(i_subject);
fname = subject_info.meg{i_subject};

switch data_type

    case 'alpha'
        % Read in the data segmented around targets
        data_dir = [exp_dir 'tfr/lf_4cyc/'];
        fn = [data_dir fname '/low_standard'];
        d = load(fn);
        d = d.low_freq_data;
        % Average over power in the alpha band
        cfg_sel.frequency = [7 14];
        cfg_sel.avgoverfreq = 'yes';
        % Make a function to grab the data for this timepoint
        get_data = @(d, i_t) d.powspctrm(:,:,:,i_t);
        % Function to get the time vector
        get_time = @(d) d.time;

    case 'raw'
        % Read in the data segmented around targets
        data_dir = [exp_dir 'preproc/target/'];
        fn = [data_dir fname '/preproc'];
        d = load(fn);
        d = d.data;
        % Only keep full trials (where no samples have been excluded)
        trial_lengths = cellfun(@length, d.time);
        full_trials = trial_lengths == max(trial_lengths);
        % Filter data
        cfg = [];
        cfg.trials = full_trials;
        cfg.detrend = 'yes';
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [0.5 40];
        cfg.padding = 2;
        d = ft_preprocessing(cfg, d);
        % Downsample the data
        cfg = [];
        cfg.resamplefs = 100;
        cfg.detrend = 'no';
        d = ft_resampledata(cfg, d);
        % Function to grab the data for this timepoint
        get_data = @get_data_raw;
        % Function to get the time vector
        get_time = @(d) d.time{1};
                
    case 'rft'
        error('Not implemented for RFT data yet')
        % Read in the data
        xxxxx
        % Make a function to grab the data for this timepoint
        get_data = @get_data_rft;
        % Function to get the time vector
        get_time = @(d) d.time;
    
    otherwise
        error('data_type %s not implemented', data_t)
end

rf_dir = [exp_dir 'forest/' data_type '/'];

% Select hit trials
hits = d.trialinfo(:,1) == 1;
cfg_sel.trials = hits;
cfg_sel.avgoverrept = 'no';

% Only look at gradiometers
cfg_sel.channel = {'MEG*2' 'MEG*3'};

% Make the selections
d = ft_selectdata(cfg_sel, d);

% Determine whether each trial is a left-hit or right-hit
trial_numbers = d.trialinfo(:,2);
targ_side_per_trial = behav.target_side(trial_numbers);
targ_left = strcmp(targ_side_per_trial, 'left');

% Parameters for the random forest analysis
n_trees = 500;
n_features = round(sqrt(length(d.label)));
in_bag_fraction = 2/3;

t = get_time(d);
forests = cell(size(t));
for i_t = 1:length(t)
    fprintf('%i: ', i_t)
    toc
    x = get_data(d, i_t);
    if ~all(all(isnan(x)))
        mdl = TreeBagger(n_trees, x, targ_left, ...
            'InBagFraction', in_bag_fraction, ...
            'NumPredictorsToSample', n_features, ...
            'OOBPrediction', 'on', ...
            'OOBPredictorImportance', 'on');
        forests{i_t} = mdl;
    end
end

[~,~,~] = mkdir(rf_dir, fname);
save([rf_dir fname '/rf'], 'forests', 't')
end

% Function to get a timepoint of data from raw (non-TFR) data
function dout = get_data_raw(d, i_t)
data_matrix = cat(3, d.trial{:});
dout = data_matrix(:, i_t, :);
dout = squeeze(dout)';
end

% Function to get a timepoint of data from RFT data
function dout = get_data_rft(d, i_t)
%%%%%%% ----- This should be moved, so it doesn't get computed on each loop
% Select data at each tagged frequency +/- 1 Hz
d_freq = cell(1, 2);
for i_freq = 1:2
    freq = exp_params.tagged_freqs(i_freq);
    cfg = [];
    cfg.frequency = [-1 1] + freq;
    cfg.avgoverfreq = 'yes';
    d_freq{i_freq} = ft_selectdata(cfg, d); 
end
d = ft_appendfreq([], d_freq{:});
data = d.powspctrm(:,:,:,i_t);
%
% Separate features for power at each tagged freq
% Include a feature for "rft_mapping": which side has 63 Hz
% Random forests discover interactions automatically
%
end


