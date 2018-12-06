function data = rs_preproc(fname, evt)

% fname: Base filename for this subject
% evt: Which event type to segment out (trials|targets|responses)
% reject: Which artifacts to reject: visual, photo, eye

rs_setup
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
    if strcmp(evt, 'target')
        trl = trialdef.trl.target;
        trl = trl(~isnan(trialdef.trl.target(:,1)),:);
    else
        trl = trialdef.trl.(evt);
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
%     cfg.polyremoval = 'yes';
%     cfg.polyorder = 1;
    d = ft_preprocessing(cfg);
    
    % Reject artifacts visually-identified and photodiode artifacts
    a = []; % Make the artfctdef object
    a.reject = 'nan';
    a.minaccepttim = 0.5;
    a.visual = art.visual.cfg_art.grad{i_block}.artfctdef.visual;
    a.photodiode.artifact = art.photo.photo_artfctdef{i_block};
    a.eye.artifact = art.eye.eyes_artfctdef{i_block};
    cfg = []; % Put it into a cfg
    cfg.artfctdef = a;
    d = ft_rejectartifact(cfg, d);
    
    % Downsample
    cfg = [];
    cfg.resamplefs = 250;
    d = ft_resampledata(cfg, d);

    % Reject artifact ICs
    cfg = [];
    cfg.component = art.ica.reject_comp;
    d = ft_rejectcomponent(cfg, art.ica.comp, d);
    
    % Save the preproc data in a cell obj
    data_by_block{i_block} = d;
    clear d
end

% Combine all the data
data = ft_appenddata([], data_by_block{block_info.main});
