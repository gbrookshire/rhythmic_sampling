% Compute the TFR around alpha peaks

%{
Following Spaak et al (2012)
- Only select segments in which alpha power is > 60th percentile for >= 800 ms
    - This selected ~3% of the data
- Trimmed out 1 s of data centered on these alpha bursts
- BP-filter the alpha signal
    - two-pass FIR-LS filter b/w 7 and 14 Hz (order: 426)
- TFRs of the signal
    - Hanning taper
    - 20-300 Hz in steps of 2 Hz
    - 7 cycles
%}


i_subject = 14;

rs_setup

% % Load data
% fname = subject_info.meg{i_subject};
% data_preproc = load([exp_dir 'preproc/trial/' fname '/preproc']);
% data_preproc = data_preproc.data;

data_preproc = rs_preproc_ress(i_subject, 'trial');
fsample = 1 / mean(diff(data_preproc.time{1}));

% Get the band-passed alpha oscillations
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [7 14];
cfg.bpfilttype = 'but'; %FIRLS get a warning: not recommended for neural signals
data_alpha = ft_preprocessing(cfg, data_preproc);

% Get alpha power in order to trim out high-power segments
cfg = [];
cfg.hilbert = 'abs';
data_alpha_pow = ft_preprocessing(cfg, data_alpha);

% Convenience function to concatenate all trials
conc = @(x) cat(2, x{:});

% Get the alpha power cutoffs at each channel 
thrsh_alpha = 60; % Percentile
cutoffs = prctile(conc(data_alpha_pow.trial), thrsh_alpha, 2);

% Look for contiguous segments of data with alpha power above that threshold
thrsh_time = 0.8 % seconds
thrsh_samp = thrsh_time * fsample; % minimum number of samples to keep
high_power = cellfun(@(x) x > cutoffs, ...
    data_alpha_pow.trial, ...
    'UniformOutput', false);
blobs = cellfun(@(x) bwconncomp(x), high_power, 'UniformOutput', false);
blobs_inx = cell(size(blobs)); % Cell of cells with blob indices
blobs_info = {};
for i_trial = 1:length(data_preproc.trial)
    trial_size = size(data_preproc.trial{i_trial});
    trial_time = data_preproc.time{i_trial};
    for i_blob = 1:blobs{i_trial}.NumObjects
        % Convert from element-wise to row/column-indexing
        [r_inx c_inx] = ind2sub(trial_size, blobs{i_trial}.PixelIdxList{i_blob});
        blob_subinx = [r_inx c_inx];
        % Find the blobs separately for each channel
        contains_blob = zeros(trial_size); 
        contains_blob(blobs{i_trial}.PixelIdxList{i_blob}) = 1;
        blob_length_per_channel = sum(contains_blob, 2);
        for i_chan = 1:length(data_preproc.label)
            % Cut out blobs that are too short
            if blob_length_per_channel(i_chan) < thrsh_samp
                blob_subinx(blob_subinx(:,1) == i_chan, :) = [];
                continue
            end
            % If the blob is long enough, save it 
            % Save independently for each channel
            bb = blob_subinx(blob_subinx(:,1) == i_chan, :); % Blob at one channel
            % if numel(bb) > 0
            %     keyboard
            % end
            b = [];
            b.center_time_inx = round(median(bb(:,2)));
            b.center_time = trial_time(b.center_time_inx);
            b.chan_inx = i_chan;
            b.trial = i_trial;
            b.width = blob_length_per_channel(i_chan);
            blobs_info{end+1} = b;
        end
        % % Save the remaining blobs of data
        % blobs_inx{i_trial}{i_blob} = blob_subinx;
    end
    % % Cut out empty blobs
    % blobs_inx{i_trial}( cellfun(@(c) isempty(c), blobs_inx{i_trial})) = [];
end

% % Trim out data centered on these blobs
% % Make sure to only take segments >0.5 s after stimulus onset, to avoid transients
% % Find new 0-times 
% zero_times = nan(0, 2); % Segment * TrialNum * ChanNum
% for i_trial = 1:length(blobs_inx)
%     for i_seg = 1:length(blobs_inx{i_trial})
%         trial
%     end
% end
