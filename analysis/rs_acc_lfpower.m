% Plot differences in LF power between sides as a function of hit/miss

clear variables
close all
rs_setup

tfr_dir = [exp_dir 'tfr/target/'];

% Read in all data
% Array for all data: Subj x Chan x TFRfreq x Time x TargSide x Hit 
agg_data = nan([height(subject_info), 204, 28, 41, 2, 2]);
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
    
    % Combine planar gradiometers
    fn = [exp_dir 'grad/' fname '/grad'];
    grad = load(fn);
    d.grad = grad.grad;
    cfg = [];
    cfg.method = 'sum';
    cfg.updatesens = 'yes';
    d = ft_combineplanar(cfg, d);

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

hmlgs = homologous_chans(true); % map left/right homologous channels
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
high_alpha_chans = {'204x', '192x', '194x', '191x', '172x' ...
                    '203x', '234x', '232x', '231x', '252x'};
mag_names = cellfun(@(s) ['MEG' s(1:(end-1)) '1'], ...
    high_alpha_chans, ...
    'UniformOutput', false);
cmb_grad_names = cellfun(...
    @(s) ['MEG' s(1:(end-1)) '2+' s(1:(end-1)) '3'], ...
    high_alpha_chans, ...
    'UniformOutput', false);
roi = ismember(d.label, [mag_names cmb_grad_names]);

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

    % For each side of the head
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
        % Extract the data for hits
        d_a = agg_data_orig(i_subject,chan_inx,:,:,:,2);
        d_l = d_a(:,:,:,:,1,:);
        d_r = d_a(:,:,:,:,2,:);
        % Compare hits on the left and right
        d_x = (d_l - d_r) ./ (d_l + d_r);
        % Average over channels
        d_x = nanmean(d_x, 2);
        % Average over subjects
        d_x = nanmean(d_x, 1);
        d_x = squeeze(d_x);

        % Plot the results
        i_plot = i_chan_side;
        subplot(2,1,i_plot)
        imagesc(d.time, d.freq, d_x)
        xlim(0.7 * [-1 1])
        xlabel('Time (s)')
        ylabel('Freq (Hz)')
        colorbar('EastOutside')
        colormap(cm)
        set(gca, 'YDir', 'normal')
        caxis([-1 1] * max(max(abs(d_x))) * 0.75) % Center color axis
        title(sprintf('Chans: %s', side_labels{i_chan_side}))
    end
    
    print('-dpng', fn)
end
