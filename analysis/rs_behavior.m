function behav = rs_behavior(i_subject)

% Load the behavioral data
rs_setup
fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
behav = readtable(fn);

% Add a column for recording number
trials_per_rec = size(behav,1) / 5;
rec_num = repelem((1:5)', trials_per_rec, 1);
behav.rec_num = rec_num;

% Convert logical variables to Booleans
behav.quest = strcmp(behav.hit, 'True');
behav.hit = strcmp(behav.hit, 'True');
behav.false_alarm = strcmp(behav.false_alarm, 'True');

% Convert RTs to numeric
behav.rt = str2double(behav.rt);