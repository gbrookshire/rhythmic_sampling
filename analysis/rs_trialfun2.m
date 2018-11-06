function trl = rs_trialfun2(dataset)

% For some reason some target triggers are not registered. No idea why...

rs_setup

% %% for testing
% i_subject = 5;
% i_block = 5;
% fname = subject_info.meg{i_subject};
% dataset = [exp_dir 'raw\' fname '\' num2str(i_block) '.fif'];
% %%

% read the header information and the events from the data
hdr = ft_read_header(dataset);
event = ft_read_event(dataset);

% Only keep relevant triggers
keep_trigs = {triggers.trial_start triggers.target triggers.response};
keep_trigs_inx = ismember({event.type}, keep_trigs);
event = event(keep_trigs_inx);

trials_per_block = 112;
trl_trial = nan([trials_per_block, 4]);
trl_target = nan([trials_per_block, 3]);
trl_response = nan([trials_per_block, 3]);

% Window within which responses are considered hits (in seconds)
hit_window = 1;

% Maximum trial duration in samples
max_trial_dur = exp_params.max_trial_dur * hdr.Fs;

% How many samples to keep before and after the trigger
pretrig = round(1 * hdr.Fs);
posttrig = round(1 * hdr.Fs);
% Offset of trial start from trigger
offset = -pretrig;

% Make sure there are the expected number of trials
trial_onset_inx = find(strcmp({event.type}, triggers.trial_start));
if length(trial_onset_inx) ~= trials_per_block
    error('Only %i trials found instead of expected %i', ...
        length(trial_onset_inx), trials_per_block)
end

% Check for expected number of targets
if sum(strcmp({event.type}, triggers.target)) ~= 112
    msg = 'Only found %i targets in file %s';
    warning(msg, ...
        sum(strcmp({event.type}, triggers.target)), ...
        dataset)
end

% Get information for each trial
for i_trial = 1:length(trial_onset_inx)
    inx = trial_onset_inx(i_trial);
    onset_t = event(inx).sample;
    % The next trigger is a target (expected outcome)
    if strcmp(event(inx + 1).type, triggers.target)
        target_t = event(inx + 1).sample;
        % The next trigger (i+2) is a response
        if length(event) > (inx + 2) && ...
                strcmp(event(inx + 2).type, triggers.response)
            % Get the response time and set it as trial offset time
            resp_t = event(inx + 2).sample;
            offset_t = resp_t;
            % Make sure it falls within the window to be considered a hit
            if (resp_t - target_t) / hdr.Fs < hit_window
                hit = 1;
            else
                hit = 0;
            end
        else % No response was given
            resp_t = NaN;
            offset_t = onset_t + max_trial_dur;
            hit = 0;
        end
    else % Target trigger is missing from the MEG data!
        target_t = NaN;
        resp_t = NaN;
        hit = NaN;
        offset_t = onset_t + max_trial_dur;
    end

    % Make the lines for the trl arrays
    % Trial
    trl_trial(i_trial,:) = [onset_t-pretrig offset_t+posttrig offset hit];
    % Target
    trl_target(i_trial,:) = [target_t-pretrig target_t+posttrig offset];
    % Response
    trl_response(i_trial,:) = [resp_t-pretrig resp_t+posttrig offset];
    
    clear inx onset_t offset_t target_t resp_t hit
end

% Make a structure of all the trl matrices
trl = [];
trl.trial = trl_trial;
trl.target = trl_target;
trl.response = trl_response;
trl.event = event;