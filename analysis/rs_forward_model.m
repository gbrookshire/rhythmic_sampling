function rs_forward_model(i_subject)

%%%% Check the cfg options for all of this!

rs_setup

% Read in the MRI scan
mri = ft_read_mri(subject_info.T1(i_subject));
% Segmentation
cfg = [];
cfg.write = 'no';
segmentedmri = ft_volumesegment(cfg, mri);
% Prepare the head model
cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, segmentedmri);
% Prepare the lead field
cfg = [];
cfg.grad = FREQ_DATA.grad; % Fill in with output from ft_freqanalysis
cfg.headmodel = headmodel;
cfg.reducerank = 2;
cfg.channel = FREQ_DATA.label;
cfg.grid.resolution = 1; % use a 3-D grid with a 1 cm resolution
cfg.grid.unit = 'cm';
grid = ft_prepare_leadfield(cfg);
% Save the results
save([exp_dir 'forward_model/' fname '/fm'], 'headmodel', 'grid')
