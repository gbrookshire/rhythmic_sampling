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


%% Check the behavioral thresholds

close all
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    % Load the behavioral data
    fn = [exp_dir 'logfiles/' subject_info.behav{i_subject} '.csv'];
    behav = rs_behavior(fn);
    subplot(4,4,i_subject)
    hold on
    for s = {'left' 'right'}
        for f = [63 78]
            inx = strcmp(behav.target_side, s) & (behav.target_side_freq == f);
            plot(behav.target_opacity(inx))
        end
    end
end

