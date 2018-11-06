function [trl, event] = rs_trialfun(cfg)

%%%%% 
warning('This doesn''t catch all the targets!')
%%%%%

% Segment trials from stimulus onset to response
% Adapted from the example on the Fieldtrip website

rs_setup

% read the header information and the events from the data
hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% Trial duration in samples
trial_dur = exp_params.max_trial_dur * hdr.Fs;

% Samples of stimulus onset
onset_samp = [event(strcmp(triggers.trial_start, {event.type})).sample]';

% Samples of responses
resp_samp = [event(strcmp(triggers.response, {event.type})).sample]';
% Exclude button presses before the first trial
resp_samp(resp_samp < onset_samp(1)) = [];

% For each trial, end at either the response or the end of the trial
end_samp = nan(size(onset_samp));
i_resp = 1;
i_onset = 1;
while i_onset <= length(onset_samp) % Go through each trial
    trial_start = onset_samp(i_onset);
    try % Get the time of the next response
        resp_time = resp_samp(i_resp);
    catch err % No more responses in this block
        if strcmp(err.identifier, 'MATLAB:badsubscript')
            resp_time = Inf;
        else
            rethrow(err)
        end
    end
    try % Get the start time of the following trial
        next_trial_start = onset_samp(i_onset + 1);
    catch err % No more trials
        if strcmp(err.identifier, 'MATLAB:badsubscript')
            next_trial_start = trial_start + 10;
        else
            rethrow(err)
        end
    end
    % Button-press occurred between trials
    if resp_time < trial_start
        i_resp = i_resp + 1;
    % Button press occurred during this trial
    elseif (trial_start < resp_time) && (resp_time < next_trial_start)
        end_samp(i_onset) = resp_time;
        i_resp = i_resp + 1;
        i_onset = i_onset + 1;
    % No button press on this trial
    else
        end_samp(i_onset) = trial_start + trial_dur;
        i_onset = i_onset + 1;
    end
end

% determine the number of samples before and after the trigger
pretrig = round(cfg.trialdef.prestim * hdr.Fs);
posttrig = round(cfg.trialdef.poststim * hdr.Fs);

% Offset of trial start from trigger
offset = -pretrig * ones(size(onset_samp));

trl = [onset_samp - pretrig, end_samp + posttrig, offset];