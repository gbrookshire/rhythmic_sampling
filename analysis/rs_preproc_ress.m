function data_ress = rs_preproc_ress(i_subject, segment_type)

% Preprocess the data with RESS spatial filters

rs_setup

% Load data
fname = subject_info.meg{i_subject};
data_preproc = load([exp_dir 'preproc/' segment_type '/' fname '/preproc']);
data_preproc = data_preproc.data;
behav = rs_behavior(i_subject);
ress = load([exp_dir 'ress/' fname '/ress']);
ress = ress.ress_maps;

% Apply the RESS filters
% Derive filters using oscillations at 63 Hz, and then apply those to
% oscillations at both tagged frequencies
data_l = rs_applyressfilt(data_preproc, ress.left.f63.ress);
data_l.label = {'left'};
data_r = rs_applyressfilt(data_preproc, ress.right.f63.ress);
data_r.label = {'right'};

% Combine the two RESS filters
data_ress = ft_appenddata([], data_l, data_r);
clear data_l data_r

% Add a column to trialinfo for which frequency was shown on the left
data_ress.trialinfo(:,3) = behav.freq_left(data_ress.trialinfo(:,2));