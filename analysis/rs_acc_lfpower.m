% Plot differences in LF power between sides as a function of hit/miss

clear variables
close all
rs_setup

tfr_dir = [exp_dir 'tfr/target/'];

% Read in all data
% Array for all data: Subj x Chan x TFRfreq x Time x TargSide x Hit 
agg_data = nan([height(subject_info), 304, 28, 41, 2, 2]);
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


%% Find which channels are on the left/right side

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

%% Select channels for an ROI
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

for i_subject = 0%0:height(subject_info)
    if i_subject == 0 % Plot all subjects
        i_subject = true([1 height(subject_info)]);
        fn = 'avg';
    elseif subject_info.exclude(i_subject)
        continue
    else
        fn = strrep(subject_info.meg{i_subject}, '/', '_');
    end
    fn = [exp_dir 'plots/alpha_power/' fn];

    % Compare hits and misses
    %{
    % For each side of the head
    for i_chan_side = 1:2
        for i_targ_side = 1:2
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
            d_a = agg_data_orig(i_subject,chan_inx,:,:,i_targ_side,:);
            d_h = d_a(:,:,:,:,:,2);
            d_m = d_a(:,:,:,:,:,1);
            % Compare hits and misses
            d_x = (d_h - d_m) ./ (d_h + d_m);
            % Average over channels
            d_x = nanmean(d_x, 2);
            % Average over subjects
            d_x = nanmean(d_x, 1);
            d_x = squeeze(d_x);

            % Plot the results
            i_plot = (2 * (i_targ_side - 1)) + i_chan_side;
            subplot(2,2,i_plot)
            imagesc(d.time, d.freq, d_x)
            xlim(0.7 * [-1 1])
            xlabel('Time (s)')
            ylabel('Freq (Hz)')
            colorbar('EastOutside')
            colormap(cm)
            set(gca, 'YDir', 'normal')
            caxis([-1 1] * max(max(abs(d_x))) * 0.75) % Center color axis
            title(sprintf('Chans: %s, Targ: %s', ...
                side_labels{i_chan_side}, ...
                side_labels{i_targ_side}))
        end
    end
    %}
    
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
        contour(d.time, d.freq, d_p, [0 0.01], ...
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
