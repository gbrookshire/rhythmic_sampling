% Plot differences in LF power between sides as a function of hit/miss

clear variables
close all
rs_setup

tfr_dir = [exp_dir 'tfr/lf_4cyc/'];
%%
% Read in all data
% Array for all data: Subj x Chan x TFRfreq x Time x TargSide x Hit 
agg_data = nan([height(subject_info), 202, 28, 51, 2, 2]);
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end

    % Read in the data segmented around targets
    behav = rs_behavior(i_subject);
    fname = subject_info.meg{i_subject};
    fn = [tfr_dir fname '/low_standard'];
    d = load(fn);
    d = d.low_freq_data;
    
%     % Combine planar gradiometers
%     fn = [exp_dir 'grad/' fname '/grad'];
%     grad = load(fn);
%     d.grad = grad.grad;
%     cfg = [];
%     cfg.method = 'sum';
%     cfg.updatesens = 'yes';
%     d = ft_combineplanar(cfg, d);

    % Select only gradiometers
    cfg = [];
    cfg.channel = {'MEG*2' 'MEG*3'};
    d = ft_selectdata(cfg, d);

    % Information about each trial
    trial_numbers = d.trialinfo(:,2);
    targ_side_per_trial = behav.target_side(trial_numbers);

    % Put the data together
    side_labels = {'left' 'right'};
    for i_targ_side = 1:2
        targ_side_sel = strcmp(targ_side_per_trial, side_labels{i_targ_side});
        for hit = 0:1
            hit_sel = d.trialinfo(:,1) == hit;
            trial_sel = hit_sel & targ_side_sel;
            cfg = [];
            cfg.trials = trial_sel;
            d_sub = ft_selectdata(cfg, d);
            x = d_sub.powspctrm; % Extract data
            x = nanmean(x, 1); % Avg over trials
            agg_data(i_subject, :, :, :, i_targ_side, hit+1) = x;
        end
    end
    clear d_sub
end

% Save a version of the original data
agg_data_orig = agg_data;
for fieldname = {'powspctrm' 'cumtapcnt' 'trialinfo' 'cfg' 'dimord'}
    d = rmfield(d, fieldname{1});
end
save([tfr_dir 'agg'], 'agg_data_orig', 'd');

% Find which channels are on the left/right side

hmlgs = homologous_chans(); % map left/right homologous channels
close all

% Compare channels
% Get indices of left channels along with their homologous right channels
left_inx = [];
right_inx = [];
for i_chan = 1:length(d.label)
    chan_label = d.label(i_chan);
    chan_pair_inx = find(strcmp(hmlgs(:,1), chan_label));
    if isempty(chan_pair_inx) % This channel isn't on the left side
        continue
    end
    homologous_chan_label = hmlgs(chan_pair_inx, 2);
    homologous_chan_inx = find(strcmp(d.label, homologous_chan_label));
    left_inx(end+1) = i_chan;
    right_inx(end+1) = homologous_chan_inx;
end

% Select channels for an ROI
high_alpha_chans = {'192x', '194x', '191x',...
                    '234x', '232x', '231x',};
chan_names = {};
for chan_num = 2:3 % Only keep gradiometers
    chan_names = [chan_names ...
        cellfun(@(s) ['MEG' s(1:(end-1)) num2str(chan_num)], ...
        high_alpha_chans, 'UniformOutput', false)];
end
roi = ismember(d.label, chan_names);

%% Plot the results
cm = flipud(lbmap(100, 'RedBlue'));

for i_subject = 0:height(subject_info)
    if i_subject == 0 % Plot all subjects
        i_subject = true([1 height(subject_info)]);
        fn = 'avg';
    elseif subject_info.exclude(i_subject)
        continue
    else
        fn = strrep(subject_info.meg{i_subject}, '/', '_');
    end
    fn = [exp_dir 'plots/alpha_power/' fn];

    % Compare left-hits and right-hits
    % For each side of the head
    clear d_agg
    clear d_comp
    figure(1)
    for i_chan_side = 1:2
        % Select channels on this side
        if i_chan_side == 1
            chan_side_inx = left_inx;
        elseif i_chan_side == 2
            chan_side_inx = right_inx;
        else
            error('oops')
        end
        % Convert chan selection to logical inx to combine w/ ROI
        chan_inx = false([length(d.label) 1]);
        chan_inx(chan_side_inx) = true;
        chan_inx = chan_inx & roi;
        % Extract the data for hits and misses
        d_a = agg_data_orig(i_subject,chan_inx,:,:,:,2);
        d_l = d_a(:,:,:,:,1,:);
        d_r = d_a(:,:,:,:,2,:);
        % Compare left-hits vs right-hits
        d_x = (d_l - d_r) ./ (d_l + d_r);
        % Average over channels
        d_x = nanmean(d_x, 2);
        % Hold onto this for doing simple stats
        d_comp(i_chan_side, :, :, :) = squeeze(d_x);
        % Average over subjects
        d_x = nanmean(d_x, 1);
        d_x = squeeze(d_x);
        d_agg(:,:,i_chan_side) = d_x;

        % Plot the results
        i_plot = i_chan_side;
        subplot(3,1,i_plot)
        imagesc(d.time, d.freq, d_x)
        xlim(0.7 * [-1 1])
        ylabel('Freq (Hz)')
        colorbar('EastOutside')
        colormap(cm)
        set(gca, 'YDir', 'normal')
        caxis([-1 1] * max(max(abs(d_x))) * 0.75) % Center color axis
        title(sprintf('Chans: %s', ...
            side_labels{i_chan_side}))
    end
    % Difference between L & R channels
    subplot(3,1,3)
    imagesc(d.time, d.freq, d_agg(:,:,1) - d_agg(:,:,2))
    xlim(0.7 * [-1 1])
    xlabel('Time (s)')
    ylabel('Freq (Hz)')
    colorbar('EastOutside')
    colormap(cm)
    set(gca, 'YDir', 'normal')
    caxis([-1 1] * max(max(abs(d_x))) * 0.75) % Center color axis
    title('Difference')
    print('-dpng', fn)
    
    % Simple stats: a t-value map
    if endsWith(fn, 'avg')
        d_p = nan(size(d_x));
        for i_time = 1:size(d_comp, 4)
            for i_freq = 1:size(d_comp, 3)
                a = squeeze(d_comp(1,:,i_freq,i_time));
                b = squeeze(d_comp(2,:,i_freq,i_time));
                if all(isnan(a))
                    continue
                end
                [~,p] = ttest(a,b);
                d_p(i_freq, i_time) = p;
            end
        end
        figure(2)
        imagesc(d.time, d.freq, log10(d_p))
        hold on
        contour(d.time, d.freq, d_p, [0 0.05], ...
            'color', 'r', 'LineWidth', 1.5);
        hold off
        xlim(0.7 * [-1 1])
        xlabel('Time (s)')
        ylabel('Freq (Hz)')
        colorbar('EastOutside')
        set(gca, 'YDir', 'normal')
        print('-dpng', [fn '-pvals'])
    end
    
end


%% Look at raw alpha power in a scatterplot
% Compare left-hits -- right-hits, avg'd over sensors and over time during
% a pre-stimulus window

time_win = [0.0 0.2] - 0.1;
freq_win = [7 14];

% Compare left-hits and right-hits
% For each side of the head
for i_chan_side = 1:2
    subplot(1,2,i_chan_side)
    % Select channels on this side
    if i_chan_side == 1
        chan_side_inx = left_inx;
    elseif i_chan_side == 2
        chan_side_inx = right_inx;
    else
        error('oops')
    end
    % Convert chan selection to logical inx to combine w/ ROI
    chan_inx = false([length(d.label) 1]);
    chan_inx(chan_side_inx) = true;
    chan_inx = chan_inx & roi;
    % Extract the data for left-hits and right-hits
    d_x = agg_data_orig(:,chan_inx,:,:,:,2); % Get hits
    % Average over channels
    d_x = nanmean(d_x, 2);
    % Average over the time window
    t_sel = (time_win(1) < d.time) & (d.time < time_win(2));
    d_x = d_x(:,:,:,t_sel,:,:);
    d_x = mean(d_x, 4);
    % Average over the frequencies of interest
    f_sel = (freq_win(1) < d.freq) & (d.freq < freq_win(2));
    d_x = d_x(:,:,f_sel,:,:,:);
    d_x = mean(d_x, 3);
    d_x = squeeze(d_x);
    
    % Make the plot
    lims = [min(min(d_x)) * 0.9, max(max(d_x)) * 1.1];
    loglog(d_x(:,1), d_x(:,2), 'o')
    hold on
    plot(lims, lims, '--k')
    hold off
    xlabel('Power (left hits)')
    ylabel('Power (right hits)')
    xlim(lims)
    ylim(lims)
    if i_chan_side == 1
        title('Left hemi')
    elseif i_chan_side == 2
        title('Right hemi')
    end
end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
print('-dpng', [exp_dir 'plots/coherence/phase_test'])

fn = [exp_dir 'plots/alpha_power/win_4cyc'];
fn = sprintf('%s/each_side_%.1f-%.1f.png', fn, time_win(1), time_win(2));
print('-dpng', fn)


%% Correlate pre- and post-stimulus effects

time_windows = [-0.5 -0.3; 0.4 0.6];
mi = {};
for t_win = time_windows'
    % Array for all data: Subj x Chan x TFRfreq x Time x TargSide x Hit 
    % Compute modulation index for all data simultaneously
    d_a = agg_data_orig(:,:,:,:,:,2); % Get all hits
    d_l = d_a(:,:,:,:,1,:); % Targets on the left
    d_r = d_a(:,:,:,:,2,:); % Targets on the right
    % Compare left-hits vs right-hits
    d_x = (d_l - d_r) ./ (d_l + d_r);
    % Average pre-stimulus activity
    t_sel = (d.time > t_win(1)) & (d.time < t_win(2));
    d_x = mean(d_x(:,:,:,t_sel), 4);
    % Snip out the alpha band
    f_sel = (d.freq >= 8) & (d.freq <= 14);
    d_x = mean(d_x(:,:,f_sel), 3);
    mi{end+1} = d_x;
end

r = nan([1 length(d.label)]);
p = nan([1 length(d.label)]);
for i_chan = 1:length(d.label)
    [rho,pval] = corr(mi{1}(:,i_chan), mi{2}(:,i_chan), 'rows', 'pairwise');
    r(i_chan) = rho;
    p(i_chan) = pval;
end

% Plot correlation for individual channels
i_plot = 1;
for chan_name_stem = high_alpha_chans
    for i_grad = 2:3
        subplot(4, 3, i_plot)
        chan_id = sprintf('MEG%s%i', chan_name_stem{1}(1:3), i_grad);
        chan_sel = strcmp(d.label, chan_id);
        pre_mi = mi{1}(:,chan_sel);
        post_mi = mi{2}(:,chan_sel);
        plot(pre_mi, post_mi, 'ob')
        hold on
%         plot(nanmean(pre_mi), nanmean(post_mi), 'or', 'MarkerSize', 12)
        xl = xlim;
        yl = ylim;
        plot(xl, [0 0], '--k')
        plot([0 0], yl, '--k')
        xlim(xl)
        ylim(yl)
        hold off
        text(0.1, 0.9, chan_id, 'Units', 'normalized')
        if i_plot == 1
            xlabel('Pre-stimulus MI')
            ylabel('Post-stimulus MI')
        end
        i_plot = i_plot + 1;
    end
end
fn = [exp_dir 'plots/alpha_power/corr-pre-post-MI'];
print('-dpng', fn)

% Plot topo-plots of p and r
data_types = {r p};
for i_type = 1:2
    d_mi = [];
    d_mi.label = d.label;
    d_mi.time = 0;
    d_mi.avg = data_types{i_type}';

    % Combine planar gradiometers
    fn = [exp_dir 'grad/' fname '/grad'];
    grad = load(fn);
    d_mi.grad = grad.grad;
    cfg = [];
    cfg.method = 'sum';
    cfg.updatesens = 'yes';
    d_mi = ft_combineplanar(cfg, d_mi);
    d_mi.avg = d_mi.avg / 2; % Make it average instead of sum

    subplot(1,2,i_type)
    cfg = [];
    cfg.layout = chan.grad_cmb.layout;
    if i_type == 1
        cfg.zlim = 'maxabs';
    else
        d_mi.avg = log10(d_mi.avg);
        cfg.highlight = 'on';
        cfg.highlightchannel
    end
    cfg.colorbar = 'no';
    cfg.style = 'straight';
    cfg.comment = 'no';
    cfg.shading = 'interp';
    cfg.markersymbol = '.';
    cfg.gridscale = 200;
    cfg.colorbar = 'yes';
    ft_topoplotTFR(cfg, d_mi)
end
subplot(1,2,1), title('rho')
subplot(1,2,2), title('p')

fn = [exp_dir 'plots/alpha_power/topo-corr-pre-post-MI'];
print('-dpng', fn)



%% Topoplot of left-hits vs right-hits
% Array for all data: Subj x Chan x TFRfreq x Time x TargSide x Hit 

% Compute modulation index for all data simultaneously
d_a = agg_data_orig(:,:,:,:,:,2); % Get all hits
d_l = d_a(:,:,:,:,1,:); % Targets on the left
d_r = d_a(:,:,:,:,2,:); % Targets on the right
% Compare left-hits vs right-hits
d_x = (d_l - d_r) ./ (d_l + d_r);
% Average over subjects
d_x = nanmean(d_x, 1);
d_x = squeeze(d_x);

d_mi = [];
d_mi.label = d.label;
d_mi.time = d.time;
d_mi.freq = d.freq;
d_mi.powspctrm = d_x;
d_mi.dimord = 'chan_freq_time';

% Combine planar gradiometers
fn = [exp_dir 'grad/' fname '/grad'];
grad = load(fn);
d_mi.grad = grad.grad;
cfg = [];
cfg.method = 'sum';
cfg.updatesens = 'yes';
d_mi = ft_combineplanar(cfg, d_mi);

% Plot it
cfg = [];
cfg.layout = chan.grad_cmb.layout;
cfg.xlim = [0.4 0.6]; % Post-stimulus
% cfg.xlim = [-0.4 0]; % Pre-stimulus
cfg.ylim = [7 14];
cfg.zlim = 'maxabs';
cfg.colorbar = 'no';
cfg.style = 'straight';
cfg.comment = 'no';
cfg.shading = 'interp';
cfg.markersymbol = '.';
cfg.gridscale = 200;
cfg.colorbar = 'yes';
ft_topoplotTFR(cfg, d_mi)

fn = [exp_dir 'plots/alpha_power/topo'];
print('-dpng', fn)

%% Look at hits vs misses regardless of side of space

% Compare hits vs misses
d_h = agg_data_orig(:,:,:,:,:,2); % hits
d_m = agg_data_orig(:,:,:,:,:,1); % misses
d_x = (d_h - d_m) ./ (d_h + d_m); % compare them

% TFR plot
i_plot = 1;
for d_sub = {d_h, d_m, d_x}
    d_sub = d_sub{1}(:,roi,:,:,:);
    d_sub = squeeze(nanmean(mean(mean(d_sub, 5), 2), 1));
    subplot(3,1,i_plot)
    imagesc(d.time, d.freq, d_sub)
    xlim(0.7 * [-1 1])
    ylabel('Freq (Hz)')
    colorbar('EastOutside')
    colormap(cm)
    set(gca, 'YDir', 'normal')
    switch i_plot
        case 1
            title('Hit')
        case 2
            title('Miss')
        case 3
            title('Hit vs. Miss')
    end
    i_plot = i_plot + 1;
end
caxis([-1 1] * max(max(abs(d_sub))) * 0.75) % Center color axis
xlabel('Time (s)')
fn = [exp_dir 'plots/alpha_power/hit_vs_miss_occip_chans'];
print('-dpng', fn)
close all

% Topoplot
d_x = nanmean(nanmean(d_x, 5), 1); % Average over subjects & target side
d_x = squeeze(d_x);

% Set up FT object
d_mi = [];
d_mi.label = d.label;
d_mi.time = d.time;
d_mi.freq = d.freq;
d_mi.powspctrm = d_x;
d_mi.dimord = 'chan_freq_time';

% Combine planar gradiometers
fn = [exp_dir 'grad/' fname '/grad'];
grad = load(fn);
d_mi.grad = grad.grad;
cfg = [];
cfg.method = 'sum';
cfg.updatesens = 'yes';
d_mi = ft_combineplanar(cfg, d_mi);

% Plot it
cfg = [];
cfg.layout = chan.grad_cmb.layout;
cfg.ylim = [7 14];
cfg.zlim = 'maxabs';
cfg.colorbar = 'no';
cfg.style = 'straight';
cfg.comment = 'no';
cfg.shading = 'interp';
cfg.markersymbol = '.';
cfg.gridscale = 200;
cfg.colorbar = 'yes';

cfg.xlim = [-0.4 0]; % Pre-stimulus
ft_topoplotTFR(cfg, d_mi)
fn = [exp_dir 'plots/alpha_power/hit_vs_miss_topo_pre-stim'];
print('-dpng', fn)

cfg.xlim = [0.4 0.6]; % Post-stimulus
ft_topoplotTFR(cfg, d_mi)
fn = [exp_dir 'plots/alpha_power/hit_vs_miss_topo_post-stim'];
print('-dpng', fn)

