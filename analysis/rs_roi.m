function roi = rs_roi(fname, roi_size)

% Return a selection of channels with the highest SNR for each freq
% fname: subject name
% roi_size: number of channels to include in the ROI

rs_base_path
exp_dir = [base_dir 'rhythmic_sampling_data/'];

roi = cell([1 2]);
for i_freq = 1:2
    snr = load([exp_dir 'spectra\' fname '\snr']);
    s = snr.snr{i_freq};

    % Get the best channels
    [~, inx] = sort(s.powspctrm, 'descend');
    roi{i_freq} = snr.snr{i_freq}.label(inx <= roi_size);
end