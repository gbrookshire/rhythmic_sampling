clear variables
close all
rs_setup

data_type = 'alpha';

%% Compute and save the Out-of-bag error

data_dir = [exp_dir 'forest/' data_type '/'];
switch data_type
    case 'alpha';
        % Load one TFR to give some info
        d_tfr = load([exp_dir 'tfr/lf_4cyc/' subject_info.meg{1} '/low_standard']);
        d_tfr = d_tfr.low_freq_data;
        t = d_tfr.time;

    case 'raw';
        % Time variable for downsampled raw data
        t = -1:0.01:1;
        
    otherwise
        error('data_type %s not implemented',data_type)
end

ooberr = nan(height(subject_info), length(t));
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};
    fprintf('Loading data: ')
    disp(fname)
    fn = [data_dir fname '/rf'];
    d = load(fn);

    for i_t = 1:length(d.forests)
        fprintf('%0.3f,', t(i_t))
        if isempty(d.forests{i_t})
            continue
        end
        ooberr(i_subject,i_t) = oobError(d.forests{i_t}, 'Mode', 'ensemble');
    end
    fprintf('\n')
end

save([data_dir 'ooberr'], 'ooberr', 't')

%% Plot the results

load([data_dir 'ooberr'])

subplot(2, 1, 1)

plot(t, ooberr, '-', 'color', [1 1 1] * 0.7)
hold on
plot(t, nanmean(ooberr, 1), '-r', 'LineWidth', 2)
plot([-1 1], [0.5 0.5], '--k')
hold off
xlim([-0.75 0.75])
xlabel('Time (s)')
ylabel('Accuracy')

subplot(2, 1, 2)

pval = nan(size(t));
for i_t = 1:length(t)
    x = ooberr(:,i_t);
    if all(isnan(x))
        continue
    end
    %[~, p] = ttest(x, 0.5);
    p = myBinomTest(sum(x < 0.5), length(x), 0.5);
    pval(i_t) = p;
end
plot(t, log10(pval))
hold on
plot([-1 1], log10(0.05) * [1 1], '--k')
plot([-1 1], log10(0.01) * [1 1], '--k')
hold off
xlim([-0.75 0.75])
xlabel('Time (s)')
ylabel('log_{10} p-value')

print([exp_dir 'plots/forest/' data_type], '-dpng')