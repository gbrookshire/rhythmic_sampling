function rs_preproc(i_subject, segment_type)

% Make data structures of just the photodiode recordings, with matching
% samples to the MEG data

rs_setup
fname = subject_info.meg{i_subject};
art_path = [exp_dir 'artifacts/' fname '/'];

disp('Loading artifact definitions...')
art = [];
art.ica = load([art_path 'ica']); % Artifact defs
art.visual = load([art_path 'visual']);
art.photo = load([art_path 'photodiode']);
%art.eye = load([art_path 'eye']); % eye-tracker
art.eog = load([art_path 'eog']);

data_by_block = cell(size(block_info.all));
for i_block = block_info.main

    % Read in the trial definition
    fn = [exp_dir 'trialdef/' fname '/' num2str(i_block) '.mat'];
    if ~exist(fn, 'file')
        warning('No trialdef for sub %s, block %d', fname, i_block)
        continue
    end
    trialdef = load(fn);
    trl = trialdef.trl.(segment_type);

    % Append the trial number to the trialinfo field
    trialnum = ((i_block - 1) * 112) + (1:112);
    trl(:,end+1) = trialnum;

    % When segmenting by targets, exclude trials in which there was no targ
    if strcmp(segment_type, 'target')
        trl = trl(~isnan(trl(:,1)),:);
    end

    % Preprocess the data
    fn = [exp_dir 'raw/' fname '/' num2str(i_block) '.fif'];
    cfg = [];
    cfg.dataset = fn;
    cfg.event = trialdef.trl.event;
    cfg.trl = trl;   

    % Preprocess the data
    cfg.channel = {'MISC004'}; % Only the photodiode
    d = ft_preprocessing(cfg);
    
    % Reject artifacts visually-identified and photodiode artifacts
    a = []; % Make the artfctdef object
    a.reject = 'partial';
    a.minaccepttim = 0.5;
    a.visual = art.visual.cfg_art.grad{i_block}.artfctdef.visual;
    a.photodiode.artifact = art.photo.photo_artfctdef{i_block};
    %a.eye.artifact = art.eye.eyes_artfctdef{i_block};
    a.eog.artifact = art.eog.eog_artfctdef{i_block};
    cfg = []; % Put it into a cfg
    cfg.artfctdef = a;
    d = ft_rejectartifact(cfg, d);
    
    % Save the preproc data in a cell obj
    data_by_block{i_block} = d;
    clear d
end

% Combine all the data
data = ft_appenddata([], data_by_block{block_info.main});

% Save the data
save_dir = [exp_dir 'preproc/photodiode/' segment_type '/'];
[~,~,~] = mkdir(save_dir, fname);
save([save_dir '/' fname '/preproc'], ...
    '-v7.3', ... % For files over 2GB
    '-nocompression', ... % Takes >2x as long to load without this
    'data')
end