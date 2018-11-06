function data = rs_preproc(fname, evt)

% fname: Base filename for this subject
% evt: Which event type to segment out (trials|targets|responses)

rs_setup
art_path = [exp_dir 'artifacts\' fname '\'];
trialdef_path = [exp_dir 'trialdef\' fname '\' evt '_'];

disp('Loading artifact definitions...')
art = [];
art.ica = load([art_path 'ica']); % Artifact defs
art.visual = load([art_path 'visual']);
art.photo = load([art_path 'photodiode']);

data_by_block = cell(size(block_info.all));
for i_block = block_info.main
    % Read in the trial definition
    fn = sprintf('%s%i.mat', trialdef_path, i_block);
    if ~exist(fn, 'file')
        warning('No trialdef for sub %s, block %d', fname, i_block)
        continue
    end
    cfg = load(fn);
    cfg = cfg.cfg;

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
    cfg = [];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.minaccepttim = 0.5;
    cfg.artfctdef.mag = art.visual.cfg_art.mag{i_block}.artfctdef.visual;
    cfg.artfctdef.grad = art.visual.cfg_art.grad{i_block}.artfctdef.visual;
    %%%cfg.artfctdef.photodiode.artifact = art.photo.photo_artfctdef{i_block};
    d = ft_rejectartifact(cfg, d);
    warning('Make sure artifacts are being rejected correctly!')
    
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
