% Compare power at the tagged frequency for hits on the L vs hits on the R.

clear variables
rs_setup
sides = {'left' 'right'};

segment_type = 'target'; % target | trial
% win_size = 0.1;

% win_str = sprintf('win_%.1fs', win_size);
% vers = 'raw_chans/all';
vers = 'ress/win_0.1s'

tfr_dir = [exp_dir 'tfr/' vers '/' segment_type '/'];

occip_roi = {'MEG2042' 'MEG2043' 'MEG2032' 'MEG2033'};

clear x % Subject x Channel X TaggedFreq x Time

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        x(i_subject,:,:,:) = nan;
        continue
    end
    
    % Load target-segmented HF TFR
    % Read in the data segmented around targets
    fname = subject_info.meg{i_subject};
    disp(fname)
    fn = [tfr_dir fname '/high'];
    d = load(fn);
    d = d.high_freq_data;

    % Load the behavioral data
    behav = rs_behavior(i_subject);

    % Which trials were hits
    hit_sel = d.trialinfo(:,1) == 1;
    trial_num = d.trialinfo(:, 2);

    for i_tagged_freq = 1:2
        % Get the frequencies of the target and distractor
        targ_freq = exp_params.tagged_freqs(i_tagged_freq);
        dist_freq = exp_params.tagged_freqs(mod(i_tagged_freq, 2) + 1);
        
        % Which data to select
        cfg = [];
        cfg.frequency = targ_freq + [-0.1 0.1];
        cfg.avgoverrpt = 'yes';
        cfg.nanmean = 'yes';
        
        % Hits on the LEFT with the target at THIS freq
        stim_side_sel = strcmp(behav.target_side(trial_num), 'left');
        stim_freq_sel = behav.freq_left(trial_num) == targ_freq;
        cfg.trials = hit_sel & stim_side_sel & stim_freq_sel;
        d_l = ft_selectdata(cfg, d);
        x_l = d_l.powspctrm;
        
        % Hits on the RIGHT with the target at the OTHER freq
        % (To ensure that the sensors are contalateral to the same freq)
        stim_side_sel = strcmp(behav.target_side(trial_num), 'right');
        stim_freq_sel = behav.freq_right(trial_num) == dist_freq;
        cfg.trials = hit_sel & stim_side_sel & stim_freq_sel;
        d_r = ft_selectdata(cfg, d);
        x_r = d_r.powspctrm;

        % Compare them
        x(i_subject,:,i_tagged_freq,:) = (x_l - x_r) ./ (x_l + x_r);
    end
end

for fieldname = {'powspctrm' 'cumtapcnt' 'trialinfo' 'cfg' 'dimord'}
    d = rmfield(d, fieldname{1});
end
save([tfr_dir 'agg'], 'x', 'd')

%% Plot it

for i_subject = 0:height(subject_info)
    
if i_subject == 0
    subj_sel = 1:height(subject_info);
    save_fname = 'avg';
elseif subject_info.exclude(i_subject)
    continue
else
    subj_sel = i_subject;
    save_fname = strrep(subject_info.meg{i_subject}, '/', '_');
end

for i_sensor_side = 1:2

    switch vers
        case 'ress/win_0.1s'
            % RESS side is coded as _stimulus_ side. 
            % Flip around to get to the sensor_ side.
            i_ress_side = mod(i_sensor_side, 2) + 1;
            chan_inx = i_ress_side;
        case 'raw_chans/occip'
            % Select channels for the raw-channel (non-RESS) analysis
            if i_sensor_side == 1
                chan_inx = 3:4; % Sensors on the left side
            elseif i_sensor_side == 2
                chan_inx = 1:2;
            end
        case 'raw_chans/all'
            if i_sensor_side == 1 % Right
                chan_inx = 3:4; 
            elseif i_sensor_side == 2 % Left
                chan_inx = 1:2;
            end
            chan_inx = ismember(d.label, occip_roi(chan_inx));
        otherwise
            error('Not implemented for TFR version %s', vers)
    end

    for i_tagged_freq = 1:2
        subplot(2, 2, (i_tagged_freq - 1) * 2 + i_sensor_side)
        mi = mean(x(subj_sel, chan_inx, i_tagged_freq,:), 2);
        mi_mean = nanmean(mi, 1);
        mi = squeeze(mi);
        mi_mean = squeeze(mi_mean);
        plot(d_l.time, mi, 'color', 0.7 * [1 1 1])
        hold on
        plot([-0.5 0.5], [0 0], '-k')
        plot([0 0], [min(min(mi)) max(max(mi))], '-k')
        plot(d_l.time, mi_mean, '-r', 'LineWidth', 2)
        title(sprintf('%s, %i Hz', ...
            sides{i_sensor_side}, exp_params.tagged_freqs(i_tagged_freq)));
        ylim([min(min(mi)) max(max(mi))])
        if i_subject == 0
            % Simple stats
            [h, p] = ttest(mi, 0, 'dim', 1);
            plot(d_l.time(p < 0.05), zeros([1 sum(h)]), '*b')
        end
        hold off
    end
end

save_dir = [exp_dir 'plots/accuracy/high_freq/'];
print('-dpng', [save_dir 'left-hits_vs_right-hits_' save_fname])
% break
end


%% Plot the difference between the right and left sides

switch vers
    case 'ress/win_0.1s'
        % For selecting RESS channels
        left_inx = 2;
        right_inx = 1;
    case 'raw_chans/occip'
        % For selecting raw channels
        left_inx = 3:4;
        right_inx = 1:2;
    case 'raw_chans/all'
        left_inx = ismember(d.label, occip_roi(3:4));
        right_inx = ismember(d.label, occip_roi(1:2));
    otherwise
        error('Not implemented for TFR version %s', vers)
end

for i_subject = 0:height(subject_info)
    
if i_subject == 0
    subj_sel = 1:height(subject_info);
    save_fname = 'avg';
elseif subject_info.exclude(i_subject)
    continue
else
    subj_sel = i_subject;
    save_fname = strrep(subject_info.meg{i_subject}, '/', '_');
end

% Difference between left and right sensors
for i_tagged_freq = 1:2
    subplot(2, 1, i_tagged_freq)
    mi_left = mean(x(subj_sel, left_inx, i_tagged_freq, :), 2);
    mi_right = mean(x(subj_sel, right_inx, i_tagged_freq, :), 2);
    mi_diff = mi_right - mi_left;
    mi_diff_mean = nanmean(mi_diff, 1); 
    mi_diff = squeeze(mi_diff);
    mi_diff_mean = squeeze(mi_diff_mean);
    plot(d_l.time, mi_diff, 'color', 0.7 * [1 1 1])
    hold on
    plot([-0.5 0.5], [0 0], '-k')
    plot([0 0], [min(min(mi_diff)) max(max(mi_diff))], '-k')
    plot(d_l.time, mi_diff_mean, '-r', 'LineWidth', 2)
    title(sprintf('%i Hz', ...
        exp_params.tagged_freqs(i_tagged_freq)));
    ylim([min(min(mi_diff)) max(max(mi_diff))])
    if i_subject == 0
        % Simple stats
        [h, p] = ttest(mi_diff, 0, 'dim', 1, 'alpha', 0.01);
        plot(d_l.time(logical(h)), zeros([1 sum(h)]), '*b')
    end
    hold off 
end

save_dir = [exp_dir 'plots/accuracy/high_freq/'];
print('-dpng',...
    [save_dir 'left-hits_vs_right-hits_diff_' save_fname])


% Collapse over RFT frequency for each sensor side
for i_sensor_side = 1:2
    subplot(2, 1, i_sensor_side)
    if i_sensor_side == 1
        chan_inx = left_inx;
    elseif i_sensor_side == 2
        chan_inx = right_inx;
    end
        
    mi_sub = mean(mean(x(subj_sel, chan_inx, :, :), 2), 3);
    mi_sub_mean = nanmean(mi_sub, 1);
    mi_sub = squeeze(mi_sub);
    mi_sub_mean = squeeze(mi_sub_mean);
    
    plot(d_l.time, mi_sub, 'color', 0.7 * [1 1 1])
    hold on
    plot([-0.5 0.5], [0 0], '-k')
    plot([0 0], [min(min(mi_sub)) max(max(mi_sub))], '-k')
    plot(d_l.time, mi_sub_mean, '-r', 'LineWidth', 2)
    title([sides{i_sensor_side} ' sensors']);
    ylim([min(min(mi_sub)) max(max(mi_sub))])
    if i_subject == 0
        % Simple stats
        [h, p] = ttest(mi_sub, 0, 'dim', 1, 'alpha', 0.01);
        plot(d_l.time(logical(h)), zeros([1 sum(h)]), '*b')
    end
    hold off 
end


save_dir = [exp_dir 'plots/accuracy/high_freq/'];
print('-dpng',...
    [save_dir 'left-hits_vs_right-hits_by-side_' save_fname])

end
%% Collapse over the two RFT frequencies

switch vers
    case 'ress/win_0.1.s'
        % For selecting RESS channels
        left_inx = 2;
        right_inx = 1;
    case 'raw_chans/occip'
        % For selecting raw channels
        left_inx = 3:4;
        right_inx = 1:2;
    case 'raw_chans/all'
        left_inx = ismember(d.label, occip_roi(3:4));
        right_inx = ismember(d.label, occip_roi(1:2));
    otherwise
        error('Not implemented for TFR version %s', vers)
end

mi_left = squeeze(mean(mean(x(:, left_inx, :, :), 2), 3));
mi_right = squeeze(mean(mean(x(:, right_inx, :, :), 2), 3));
mi_diff = mi_right - mi_left;
mi_diff_mean = nanmean(mi_diff, 1); 
plot(d_l.time, mi_diff, 'color', 0.7 * [1 1 1])
hold on
plot([-0.5 0.5], [0 0], '-k')
plot([0 0], [min(min(mi_diff)) max(max(mi_diff))], '-k')
plot(d_l.time, mi_diff_mean, '-r', 'LineWidth', 2)
title(sprintf('%i Hz', ...
    exp_params.tagged_freqs(i_tagged_freq)));
ylim([min(min(mi_diff)) max(max(mi_diff))])
% Simple stats
[h, p] = ttest(mi_diff, 0, 'dim', 1, 'alpha', 0.05);
plot(d_l.time(logical(h)), zeros([1 sum(h)]), '*b')
hold off 

print('-dpng',...
    [exp_dir 'plots/accuracy/high_freq/left-hits_vs_right-hits_diff-overall'])


%% Plot the topography of this effect

% x: Subject x Channel X TaggedFreq x Time
x_m = squeeze(nanmean(mean(x, 3), 1));

% Set up FT object
d_mi = [];
d_mi.label = d.label;
d_mi.time = d.time;
d_mi.avg = x_m;

% Combine planar gradiometers
fn = [exp_dir 'grad/' fname '/grad'];
grad = load(fn);
d_mi.grad = grad.grad;
cfg = [];
cfg.method = 'sum';
cfg.updatesens = 'yes';
d_mi = ft_combineplanar(cfg, d_mi);

% Plot it
t_centers = (-1 * [0.4 0.3 0.2 0.1 0]');
for i_plot = 1:length(t_centers)
    subplot(1, length(t_centers), i_plot)
    cfg = [];
    cfg.layout = chan.grad_cmb.layout;
    cfg.xlim = t_centers(i_plot) + [-0.01 0.01];
    cfg.zlim = 'maxabs';
    cfg.colorbar = 'no';
    cfg.style = 'straight';
    cfg.comment = 'no';
    cfg.shading = 'interp';
    cfg.markersymbol = '.';
    cfg.gridscale = 200;
    cfg.colorbar = 'yes';
    cfg.highlight = 'on';
    cfg.highlightchannel = {'MEG2042+2043' 'MEG2032+2033'};
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 10;
    ft_topoplotER(cfg, d_mi)
    title(sprintf('t = [%.2f, %.2f]', cfg.xlim(1), cfg.xlim(2)));
end

fn = [exp_dir 'plots/accuracy/high_freq/left-hits_vs_right-hits_topo'];
print('-dpng', fn)

