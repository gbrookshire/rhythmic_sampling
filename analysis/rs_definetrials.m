function rs_definetrials(i_subject)

% Define the trials and save the output

rs_setup

fname = subject_info.meg{i_subject};
[~,~,~] = mkdir([exp_dir 'trialdef/'], fname);
trl_dir = [exp_dir 'trialdef/' fname '/'];

for i_block = block_info.all
    dataset = [exp_dir 'raw/' fname '/' num2str(i_block) '.fif'];
    if ~exist(dataset, 'file')
        warning('No MEG data file: %s\\%d.fif', fname, i_block)
        input('Press ENTER to continue')
        continue
    end

    trl = rs_trialfun(dataset);
    save([trl_dir num2str(i_block)], 'trl')
 end
