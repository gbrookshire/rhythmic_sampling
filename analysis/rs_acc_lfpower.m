% Plot differences in LF power between sides as a function of hit/miss

clear variables
close all
rs_setup

% Which TFR window to use?
win_size = 0.2;
win_str = sprintf('win_%.1fs', win_size);
tfr_dir = [exp_dir 'tfr/' win_str '/'];

% Array for all data: Subj x Chan x TFRfreq x Time x TargSide x Hit 
agg_data = nan([height(subject_info), 204, 11, 101, 2, 2]);

for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end

    % Read in the data segmented around targets
    behav = rs_behavior(i_subject);
    fname = subject_info.meg{i_subject};
    fn = [tfr_dir 'target/' fname '/low'];
    d = load(fn);
    d = d.low_freq_data;
    
    % Convert from complex fourier spec to power spec
    d.powspctrm = abs(d.fourierspctrm) .^ 2;
    d = rmfield(d, 'fourierspctrm');
    
    % Read in the grad structure (to be able to combine grads)
    fn = [exp_dir 'grad/' fname '/grad'];
    grad = load(fn);
    d.grad = grad.grad;
    % Combine planar gradiometers
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

% Save a version of the original data
agg_data_orig = agg_data;


%% Normalize the power within subjects on each side
%{
agg_data = agg_data_orig; % Set to the original data
for i_subject = 1:height(subject_info)
    for i_chan = 1:length(d.label)
        for i_freq = 1:length(d.freq)
            x = agg_data(i_subject,i_chan,i_freq,:,:,:);
            sd_x = std(reshape(x, [1 numel(x)]), 'omitnan');
            m_x = nanmean(reshape(x, [1 numel(x)]));
            z_x = (x - m_x) / sd_x;
            agg_data(i_subject,i_chan,i_freq,:,:,:) = z_x;
        end
    end
end
%}

%% Select channels for an occipital ROI
occip_chans = {'191x' '231x' '201x' '202x' '203x' '204x'...
    '194x' '192x' '234x' '232x'};
mag_names = cellfun(@(s) ['MEG' s(1:(end-1)) '1'], ...
    occip_chans, ...
    'UniformOutput', false);
cmb_grad_names = cellfun(...
    @(s) ['MEG' s(1:(end-1)) '2+' s(1:(end-1)) '3'], ...
    occip_chans, ...
    'UniformOutput', false);
roi = ismember(d.label, [mag_names cmb_grad_names]);


% Plot the results
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

    i_plot = 1;
    for hit = [0 1] 
        for i_targ_side = 1:2 

            % Find which channels are ipsi/contra-lateral to the target
            switch side_labels{i_targ_side}
                case 'left'
                    ipsi_inx = left_inx;
                    contra_inx = right_inx;
                case 'right'
                    ipsi_inx = right_inx;
                    contra_inx = left_inx;
                otherwise
                    error('oops')
            end
            
            % Convert chan selection to logical inx to combine w/ ROI
            i_inx = false([1 length(d.label)]);
            i_inx(ipsi_inx) = true;
            i_inx = i_inx' & roi;
            c_inx = false([1 length(d.label)]);
            c_inx(contra_inx) = true;
            c_inx = c_inx' & roi;
            clear ipsi_inx contra_inx

            % Compute comparison: (Ipsi - Contra) / (Ipsi + Contra)
            x_ipsi = agg_data(:,i_inx,:,:,i_targ_side,hit+1);
            x_contra = agg_data(:,c_inx,:,:,i_targ_side,hit+1);
            x_comp = (x_ipsi - x_contra) ./ (x_ipsi + x_contra);
            % x_comp = (x_ipsi - x_contra);
            

            % Plot the results
            subplot(2,2,i_plot)
            x = squeeze(mean(nanmean(x_comp(i_subject,:,:,:), 1), 2));
            imagesc(d.time, d.freq, x)
            set(gca, 'YDir', 'normal')
            title(sprintf('Targ: %s, Hit: %i', side_labels{i_targ_side}, hit))
%             % Fix the color scale so it's constant between all plots
%             caxis([min(reshape(x_comp, [1 numel(x_comp)])) ...
%                 max(reshape(x_comp, [1 numel(x_comp)]))])
            %caxis([-1 1] * 0.5)
            i_plot = i_plot + 1;
        end 
    end
    
    print('-dpng', fn)
end












