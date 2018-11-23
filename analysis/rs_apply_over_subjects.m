function rs_apply_over_subjects(func_handle, parallel, subj_selection)

% Run a function for every subject
%
% func_handle: @-handle of which function to run
% parallel: true to run with `parfor`, or false to use a serial for-loop
% subj_selection: Which subjects to run (default: all)

rs_setup

% If no subject selection is provided, use all non-excluded subjects
if nargin < 3
    subj_selection = 1:height(subject_info);
    subj_selection(logical(subject_info.exclude)) = [];
end

if parallel
    nworkers = length(subj_selection);
else
    nworkers = 0;
end

% Run the function for each subject
parfor (subj_inx = 1:length(subj_selection), nworkers)
    % Indices to parfor need to be sequential integers, so we have to make
    % an indexing index here.
    i_subject = subj_selection(subj_inx);
    feval(func_handle, i_subject);
end