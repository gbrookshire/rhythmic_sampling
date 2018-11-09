function roi = rs_roi(subj)

% Return a selection of channels with the highest SNR

roi_size = 4; % Number of channels to return

% Take the SNR at 63 Hz, because it gives a clearer signal
snr = load([exp_dir 'spectra\' fname '\snr']);
s = snr.snr{1};

% Get the best channels
[~, inx] = sort(s.powspctrm, 'descend');

roi = snr.snr{1}.label(inx <= roi_size);
