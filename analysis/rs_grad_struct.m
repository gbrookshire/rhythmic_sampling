function rs_grad_struct(i_subject)

% Save grad structures for each subject (only from one block of trials)

rs_setup

fname = subject_info.meg{i_subject};
dataset = [exp_dir 'raw/' fname '/3.fif'];
hdr = ft_read_header(dataset);
grad = hdr.grad;
[~,~,~] = mkdir([exp_dir 'grad/'], fname);
save([exp_dir 'grad/' fname '/grad'], 'grad')
