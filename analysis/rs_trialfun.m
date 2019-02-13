function trl = rs_trialfun(dataset)

rs_setup

% for testing
%{
i_subject = 1;
i_block = 3;
fname = subject_info.meg{i_subject};
dataset = [exp_dir 'raw/' fname '/' num2str(i_block) '.fif'];
%}

% read the header information and the events from the data
hdr = ft_read_header(dataset);
event = ft_read_event(dataset);

% Only keep relevant triggers
keep_trigs = {triggers.trial_start triggers.target triggers.response};
keep_trigs_inx = ismember({event.type}, keep_trigs);
event = event(keep_trigs_inx);

trials_per_block = 112;
trl_trial = nan([trials_per_block, 4]);
trl_target = nan([trials_per_block, 4]);
trl_response = nan([trials_per_block, 4]);

% Time window within which responses are considered hits (seconds)
hit_window = 1;

% Maximum trial duration (samples)
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

% Get sample information for each trial
for i_trial = 1:length(trial_onset_inx)
    inx = trial_onset_inx(i_trial); % Index of the trigger for trial onset
    begin_t = event(inx).sample; % Trial beginning sample
    
    % If this trial was a false alarm
    if strcmp(event(inx + 1).type, triggers.response)
        % Exclude these trials by setting to NaNs in the target trl
        target_t = NaN;
        resp_t = event(inx + 1).sample;
        hit = -1;
        end_t = event(inx + 1).sample; % Trial end sample
    
    % The next trigger is a target
    elseif strcmp(event(inx + 1).type, triggers.target)
        target_t = event(inx + 1).sample;
    
        % The trigger after the target (i+2) is a response
        if length(event) >= (inx + 2) && ...
                strcmp(event(inx + 2).type, triggers.response)
            % Get the response time and set it as trial offset time
            resp_t = event(inx + 2).sample;
            end_t = resp_t;
        
            % Make sure it falls within the window to be considered a hit
            rt = (resp_t - target_t) / hdr.Fs;
            if rt < hit_window % It's a hit
                hit = 1;
            else % There was a response, but it was too slow
                hit = 2;
                if rt > exp_params.max_trial_dur %Resp wasn't in this trial
                    resp_t = NaN;
                    end_t = begin_t + 1000;
                end
            end
            
        else % No response was given -- miss
            hit = 0;
            resp_t = NaN;
            end_t = begin_t + max_trial_dur;
        end
        
    else % No target or response in this trial (Shouldn't happen)
        msg = sprintf(...
            'Missing target or response trigger in file %s', ...
            dataset);
        warning(msg);
 
        target_t = NaN;
        resp_t = NaN;
        hit = NaN;
        end_t = begin_t + 1; % Make the trial 1 sample long
    end

    keyboard
    % Make the lines for the trl arrays
    % Trial
    trl_trial(i_trial,:) = [begin_t-pretrig end_t+posttrig offset hit];
    % Target
    trl_target(i_trial,:) = [target_t-pretrig target_t+posttrig offset hit];
    % Response
    trl_response(i_trial,:) = [resp_t-pretrig resp_t+posttrig offset hit];
    
    clear inx onset_t offset_t target_t resp_t hit
end

% Make a structure of all the trl matrices
trl = [];
trl.trial = trl_trial;
trl.target = trl_target;
trl.response = trl_response;
trl.event = event;

%% Test the trl objects -- make sure they're the expected length
%{
rs_setup
for i_subject = 1:height(subject_info)
    disp(['Subject: ' num2str(i_subject)])
    if subject_info.exclude(i_subject)
        continue
    end
    
    fname = subject_info.meg{i_subject};
    trl_dir = [exp_dir 'trialdef\' fname '\'];

    for i_block = block_info.all
        x = load([trl_dir num2str(i_block)]);
        trial_length = x.trl.trial(:,2) - x.trl.trial(:,1);
        too_big = trial_length > 7000;
        if sum(too_big) > 0
            disp(['Weird: ' num2str(find(too_big))])
            keyboard
        end
    end

end
%}