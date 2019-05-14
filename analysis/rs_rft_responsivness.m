function rri = rs_rft_responsivness()

% Get an index of RFT strength

clear variables
rs_setup

win_size = 0.2;
win_str = sprintf('win_%.1fs', win_size);

% Load the data
overall_data = nan(height(subject_info), 46); % Subj * Freq
for i_subject = 1:height(subject_info)
    if subject_info.exclude(i_subject)
        continue
    end
    
    fname = subject_info.meg{i_subject};
    disp(fname)
    hf_fname = [exp_dir 'tfr/' win_str '/trial/' fname '/high.mat'];
    %finfo = dir(hf_fname);
    %disp(finfo.date)
    hf = load(hf_fname);
    freq_data = hf.high_freq_data;
    clear hf;
    
    % Average over trials and channels
    % Only keep time after the onset transent, before the end
    cfg = [];
    cfg.avgoverrpt = 'yes';
    cfg.nanmean = 'yes';
    cfg.avgoverchan = 'yes'; % L&R RESS channels
    cfg.latency = [1 3];
    cfg.avgovertime = 'yes';
    freq_data = ft_selectdata(cfg, freq_data);
    overall_data(i_subject,:) = freq_data.powspctrm;
end

% Make the RFT responsiveness index
approx_eq = @(x,y) abs(x - y) < 0.1;
inx_63 = approx_eq(freq_data.freq, 63);
inx_78 = approx_eq(freq_data.freq, 78);
% keep_subjs = ~isnan(overall_data(:,1));
x = overall_data(:,:);

rri = mean(log10([x(:,inx_63), x(:,inx_78)]), 2);
end