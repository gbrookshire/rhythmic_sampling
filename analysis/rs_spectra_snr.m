function rs_spectra_snr(i_subject)

% Compute the SNR over the spectra for each trial

rs_setup

fname = subject_info.meg{i_subject};
spec = load([exp_dir 'spectra/' fname '/spectra']);
spec = spec.freq_data;

f_tag = exp_params.tagged_freqs; % Tagged frequencies

snr = cell(size(f_tag));
for i_freq = 1:length(f_tag)
    snr{i_freq} = rs_snr(spec, f_tag(i_freq));
end

save([exp_dir 'spectra/' fname '/snr'], 'snr')
