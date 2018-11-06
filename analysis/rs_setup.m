exp_dir = 'C:\Users\brookshg\Documents\rhythmic_sampling\data\';
addpath([exp_dir '..\analysis'])

% Genearal experiment parameters
exp_params = [];
exp_params.tagged_freqs = [63 78]; % Hz
exp_params.max_trial_dur = 4; % seconds

% Triggers
triggers = [];
triggers.trial_start = 'STI002';
triggers.target = 'STI003';
triggers.response = 'STI015';

% Channel details
chan = [];
chan.bad.names = {'-MEG1843' '-MEG2122'};
chan.all.names = ['MEG*' chan.bad.names];
chan.all.layout = 'neuromag306all.lay';
chan.mag.names = ['MEG*1' chan.bad.names];
chan.mag.layout = 'neuromag306mag.lay';
chan.grad.names = [{'MEG*2' 'MEG*3'} chan.bad.names];
chan.grad.layout = 'neuromag306planar.lay';
chan.grad_cmb.layout = 'neuromag306cmb.lay';
chan.eog.names = {'BIO001' 'BIO002' 'BIO003'}; % EOG & EKG
chan.eyetracker.names = {'MISC001' 'MISC002'};

% Block info
block_info = [];
block_info.all = 1:5;
block_info.quest = 1:2;
block_info.main = 3:5;

% Read in subject info
subject_info = readtable([exp_dir 'subject_info.csv'], 'Delimiter', ',');