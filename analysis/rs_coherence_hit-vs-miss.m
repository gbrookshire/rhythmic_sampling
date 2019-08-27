% Compute coherence between the MEG recordings and the photodiode
% 
% Overall -- just check how coherence is working
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We only have photodiode recordings for the stim on the right %
% To get around this, we can synthesize stim signals on the left
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
close all
rs_setup

% How should the hits and misses be compared
% cmp_fnc = @(a,b) a - b;
cmp_fnc = @(a,b) (a - b) ./ (a + b);

% The stimuli were not shown at exactly the intended frequencies, due to
% the fact that the screen was not refreshed at exactly 120 Hz. These
% adjustments here fix the difference in frequencies for the synth signals.
freq_adjust = [-0.0535 -0.0648]; % For 63 and 78 Hz, respectively

clear cohspctrm_overall % Dim: subject * [cohspctrm]

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    fname = subject_info.meg{i_subject};

    % Load preprocessed data
    behav = rs_behavior(i_subject);
    d_meg = rs_preproc_ress(i_subject, 'target');
    fn_photo = [exp_dir 'preproc/photodiode/target/' fname '/preproc'];
    d_photo = load(fn_photo);
    d_photo = d_photo.data;
    d_photo.label = {'stim_right'};

    % Synthesize surrogate photodiode response for the
    % non-measured left stim
    d_synth = d_photo;
    d_synth.label = {'stim_left'};
    trial_num = d_meg.trialinfo(:,2);
    for i_trial = 1:length(d_photo.trial)
        % Adjust freq using measured diff in freq due to the frame
        % rate being slightly slower than 120 Hz
        f = behav.freq_left(trial_num(i_trial)); % Freq on left in this trial
        f = f + freq_adjust(exp_params.tagged_freqs == f);
        t = d_photo.time{i_trial};
        y = sin(2 * pi * f * t);
        d_synth.trial{i_trial} = y;
    end
    clear f t y
    
    % Combine the data structures
    data = ft_appenddata([], d_meg, d_photo, d_synth);
    data.trialinfo = d_meg.trialinfo;

    % Compute TFR 
    cfg = [];
    cfg.output = 'fourier'; % Full Fourier representation
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.foi = 55:100;
    cfg.toi = -0.5:0.05:0.5;
    cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.1;
    cfg.keeptrials = 'yes';
    fourier = ft_freqanalysis(cfg,data);

    % Get rid of trials with nans
    % And keep only hits
    nan_trials = any(squeeze(isnan(fourier.fourierspctrm(:,1,1,:))), 2);
    % hits = fourier.trialinfo(:,1) == 1;
    cfg = [];
    cfg.trials = ~nan_trials;
    fourier = ft_selectdata(cfg, fourier);

    % Find which side each RFT frequency was on in each trial
    trial_num = fourier.trialinfo(:,2);
    freq_left = behav.freq_left(trial_num);

    % Get coherence separately for each mapping of freq to each side
    clear cohspctrm % chancmb * freq * time * RFT_mapping * hit
    clear trial_counts
    for i_rft_map = 1:2
        for hit = 0:1
            % i_rft_map == 1: left side at 63 Hz, right side at 78 Hz
            % i_rft_map == 2: left side at 78 Hz, right side at 63 Hz
            hit_sel = fourier.trialinfo(:,1) == hit;
            freq_sel = freq_left == exp_params.tagged_freqs(i_rft_map);
            trial_sel = hit_sel & freq_sel;

            trial_counts(i_rft_map,hit+1) = sum(trial_sel);
            % Compute coherence
            cfg = [];
            cfg.trials = trial_sel;
            cfg.method = 'coh';
            cfg.channelcmb = {'left' 'stim_left';
                              'right' 'stim_right'};
            coh = ft_connectivityanalysis(cfg, fourier);
            cohspctrm(:,:,:,i_rft_map,hit+1) = coh.cohspctrm;
            cohspctrm_overall(i_subject,:,:,:,i_rft_map,hit+1) = coh.cohspctrm;
        end
    end

    % Clean out the raw coherence spectrum so we don't get confused
    coh = rmfield(coh, 'cohspctrm');

    % Remember that the RESS filters are named after the stimulus that 
    % they're designed to select. So 'left' refers to a filter that 
    % selects stimuli on the left side -- not the channels on the 
    % left side of the head.

    % Plot it
    for i_chancmb = 1:2
        figure(i_chancmb)
        for i_rft_map = 1:2
            for hit = 0:2
                subplot(3, 2, (hit * 2) + i_rft_map)

                if hit == 2 % Compare Hits vs Misses
                    x_m = cohspctrm(i_chancmb, :, :, i_rft_map,1);
                    x_h = cohspctrm(i_chancmb, :, :, i_rft_map,2);
                    x = cmp_fnc(x_h, x_m); % Set at the top of the script
                else % For hits or misses alone
                    x = cohspctrm(i_chancmb, :, :, i_rft_map,hit+1);
                end
                x = squeeze(x);

                imagesc(coh.time, ...
                    coh.freq, ...
                    x)
                set(gca, 'YDir', 'normal')
                colorbar;
            end
        end

        subplot(3,2,1)
        hold on
        t = text(-1, 80, 'Miss');
        t.FontSize = 12;
        title('Right stim: 78 Hz')
        hold off

        subplot(3,2,2)
        title('Right stim: 63 Hz')

        subplot(3,2,3)
        hold on
        t = text(-1, 80, 'Hit');
        t.FontSize = 12;
        hold off

        subplot(3,2,5)
        hold on
        t = text(-1, 80, 'Diff');
        t.FontSize = 12;
        hold off

        % Save a plot of this
        fn = strrep(fname, '/', '_');
        fn = [fn '-' coh.labelcmb{i_chancmb} '_stim'];
        save_fname = [exp_dir 'plots/coherence/hit_miss/' fn];
        print('-dpng', save_fname)
    end
end

%% Plot the average over subjects

for i_chancmb = 1:2
    figure(i_chancmb)
    for i_rft_map = 1:2
        for hit = 0:2
            subplot(3, 2, (hit * 2) + i_rft_map)

            if hit == 2 % Hits vs Misses
                x_m = cohspctrm_overall(:,i_chancmb, :, :, i_rft_map,1);
                x_h = cohspctrm_overall(:,i_chancmb, :, :, i_rft_map,2);
                x = cmp_fnc(x_h, x_m);
            else
                x = cohspctrm_overall(:,i_chancmb, :, :, i_rft_map,hit+1);
            end
            x = nanmean(x, 1);
            x = squeeze(x);

            imagesc(coh.time, ...
                coh.freq, ...
                x)
            set(gca, 'YDir', 'normal')
            colorbar;
        end
    end

    subplot(3,2,1)
    hold on
    t = text(-1, 80, 'Miss');
    t.FontSize = 12;
    title('Right stim: 78 Hz')
    hold off

    subplot(3,2,2)
    title('Right stim: 63 Hz')

    subplot(3,2,3)
    hold on
    t = text(-1, 80, 'Hit');
    t.FontSize = 12;
    hold off

    subplot(3,2,5)
    hold on
    t = text(-1, 80, 'Diff');
    t.FontSize = 12;
    hold off

    % Save a plot of this
    fn = ['avg' '-' coh.labelcmb{i_chancmb} '_stim'];
    save_fname = [exp_dir 'plots/coherence/hit_miss/' fn];
    print('-dpng', save_fname)
end