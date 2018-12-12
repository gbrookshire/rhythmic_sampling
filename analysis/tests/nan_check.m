function nan_check(data)
% Check for NaNs in a fieldtrip data structure

trial_len = cellfun(@(x) size(x,2), data.trial);
nan_counts = zeros([1 max(trial_len)]);
trial_counts = nan_counts;
for i_trial = 1:length(data.trial)
    x = data.trial{i_trial}(1,:); % Look for NaNs in one channel
    nan_counts(1:length(x)) = nan_counts(1:length(x)) + isnan(x);
    trial_counts(1:length(x)) = trial_counts(1:length(x)) + ~isnan(x);
end
subplot(2,1,1)
plot(nan_counts), ylabel('NaN count')
subplot(2,1,2)
plot(trial_counts), ylabel('Retained trial count')
end