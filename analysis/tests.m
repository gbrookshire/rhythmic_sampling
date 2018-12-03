% Make sure the trialdef is giving sensible trial lengths

rs_setup
reclen = 112; % Trials per recording
for i_subject = 1:height(subject_info)
    lens = nan([ 1 336*5]);
    if subject_info.exclude(i_subject)
        continue
    end
    for i_rec = 1:5
        n = num2str(i_rec);
        x = load([exp_dir 'trialdef/' subject_info.meg{i_subject} '/' n]);
        trllen = x.trl.trial(:,2) - x.trl.trial(:,1);
        lens((1:reclen) + ((i_rec - 1) * reclen)) = trllen;
    end
    subplot(4,4,i_subject)
    hist(lens, 40)
end

subplot(4,4,13)
ylabel('Count')
xlabel('Length (samples)')


%% Cross-check hits & FAs in the behavioral & trl structs
% The presentation script coded late responses as FAs, but they should be
% noted as a distinct kind of response. Don't use the behavioral data for
% FAs -- instead, use the 4th column of the trial-based trl

clear variables
rs_setup
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    fprintf('\n%i - %s: ', i_subject, fname);
    
    % Load the behavioral data
    fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
    behav = rs_behavior(fn);
    
    for i_rec = block_info.all
        fn = [exp_dir 'trialdef/' fname '/' num2str(i_rec) '.mat'];
        trialdef = load(fn);
        trial_inx = ((112 * (i_rec - 1)) + 1):(112 * i_rec);
        behav_rec = behav(trial_inx,:);
        
        behav_fa = behav_rec.false_alarm == 1;
        trl_fa = trialdef.trl.trial(:, 4) == -1;
        fa_mismatch = trl_fa ~= behav_fa;
        
        trl_not_hit = ismember(trialdef.trl.trial(:,4)', [-1 2]);
        behav_hit = behav_rec.hit(~trl_not_hit) == 1;
        trl_hit = trialdef.trl.trial(~trl_not_hit, 4) == 1;
        hit_mismatch = trl_hit ~= behav_hit;
        
        behav_rec.trl = trialdef.trl.trial(:,4);
        behav_rec.fa_mm = fa_mismatch;
        behav_rec.hit_mm(~trl_not_hit) = hit_mismatch;

%         fas = behav.false_alarm(trial_inx);
%         hit_mismatch = behav_hit(~fas) ~= ;
%         trl_fa = trialdef.trl.trial(:,4) == -1;
%         fa_mismatch = fas ~= trl_fa ;

        % Print the number of mismatching hits/misses and false alarms
        fprintf('%i/%i, ', sum(hit_mismatch), sum(fa_mismatch))
    end
end
fprintf('\n')


