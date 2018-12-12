function data = rs_preproc(i_subject, segment_type)

% fname: Base filename for this subject
% evt: Which event type to segment out (trials|targets|responses)
% reject: Which artifacts to reject: visual, photo, eye

rs_setup
fname = subject_info.meg{i_subject};
art_path = [exp_dir 'artifacts/' fname '/'];

disp('Loading artifact definitions...')
art = [];
art.ica = load([art_path 'ica']); % Artifact defs
art.visual = load([art_path 'visual']);
art.photo = load([art_path 'photodiode']);
art.eye = load([art_path 'eye']);

data_by_block = cell(size(block_info.all));
for i_block = block_info.main

    % Read in the trial definition
    fn = [exp_dir 'trialdef/' fname '/' num2str(i_block) '.mat'];
    if ~exist(fn, 'file')
        warning('No trialdef for sub %s, block %d', fname, i_block)
        continue
    end
    trialdef = load(fn);
    
    % When segmenting by targets, exclude trials in which there was no targ
    if strcmp(segment_type, 'target')
        trl = trialdef.trl.target;
        trl = trl(~isnan(trialdef.trl.target(:,1)),:);
    else
        trl = trialdef.trl.(segment_type);
    end

    % Preprocess the data
    fn = [exp_dir 'raw/' fname '/' num2str(i_block) '.fif'];
    cfg = [];
    cfg.dataset = fn;
    cfg.event = trialdef.trl.event;
    cfg.trl = trl;   

    % Preprocess the data
    cfg.channel = art.ica.comp.cfg.channel; % Select the good channels
    cfg.bsfilter = 'yes'; % Band-stop filter to get rid of line noise
    cfg.bsfreq = [48 52];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 100; % Changed from 80 in FLUX pipeline (too close to 78)
    cfg.padding = 8; % Pad the data to reduce filtering artifacts
    cfg.padtype = 'data';
    d = ft_preprocessing(cfg);
    %nan_check(d) %% No NaNs at this stage
    
    % Reject artifacts visually-identified and photodiode artifacts
    a = []; % Make the artfctdef object
    a.reject = 'nan';
    a.minaccepttim = 0.5;
    a.visual = art.visual.cfg_art.grad{i_block}.artfctdef.visual;
    a.photodiode.artifact = art.photo.photo_artfctdef{i_block};
    a.eye.artifact = art.eye.eyes_artfctdef{i_block};
    cfg = []; % Put it into a cfg
    cfg.artfctdef = a;
    % Check whether artifacts are as expected
    %{
    cfg.viewmode = 'vertical';
    cfg.layout = chan.grad.layout;
    ft_databrowser(cfg, d);
    %}
    d = ft_rejectartifact(cfg, d);
    %nan_check(d) %% More NaNs than expected
    
    % Downsample
    cfg = [];
    cfg.resamplefs = 250;
    d = ft_resampledata(cfg, d);
    %nan_check(d) %% No big change from after rejecting artifacts

    % Reject artifact ICs
    cfg = [];
    cfg.component = art.ica.reject_comp;
    d = ft_rejectcomponent(cfg, art.ica.comp, d);
    %nan_check(d) %% Same as after downsampling
    
    % Save the preproc data in a cell obj
    data_by_block{i_block} = d;
    clear d
end

% Combine all the data
data = ft_appenddata([], data_by_block{block_info.main});

% Save the data
save_dir = [exp_dir 'preproc/' segment_type '/'];
[~,~,~] = mkdir(save_dir, fname);
save([save_dir '/' fname '/preproc'], 'data')

end
