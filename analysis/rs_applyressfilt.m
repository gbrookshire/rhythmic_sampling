function data_out = rs_applyressfilt(data_in, ress_filter)

% Apply a RESS/GEDb spatial filter to a Fieldtrip data structure

% Make Fieldtrip-style object to return out
data_out = [];
data_out.label = {'RESS'};
data_out.trialinfo = data_in.trialinfo;
data_out.time = data_in.time;
data_out.trial = cell(size(data_in.trial));
for i_trial = 1:length(data_in.trial)
    trial_data = data_in.trial{i_trial};
    data_out.trial{i_trial} = ress_filter * squeeze(trial_data);
end