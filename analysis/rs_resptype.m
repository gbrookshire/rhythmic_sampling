function [hits, nans] = rs_resptype(i_subject)

% Return a vector of hits/misses/etc

rs_setup
fname = subject_info.meg{i_subject};

hits = [];
nans = [];
for i_rec = block_info.main
    % Read in the trial definition
    fn = [exp_dir 'trialdef/' fname '/' num2str(i_rec) '.mat'];
    trialdef = load(fn);
    hits = [hits trialdef.trl.trial(:,4)'];
    nans = [nans isnan(trialdef.trl.target(:,1))'];
end