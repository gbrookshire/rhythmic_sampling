
% Compute a measure of the effect in the low frequencies
% L-Hit vs R-Hit at (the slightly different selection of) occipital chans
lf = load([exp_dir 'tfr/lf_0.3s/' 'agg']);
% Array for all data: Subj x Chan x TFRfreq x Time x TargSide x Hit 
% Compute modulation index for all data simultaneously
d_a = lf.agg_data_orig(:,:,:,:,:,2); % Get all hits
d_l = d_a(:,:,:,:,1,:); % Targets on the left
d_r = d_a(:,:,:,:,2,:); % Targets on the right
% Compare left-hits vs right-hits
d_x = (d_l - d_r) ./ (d_l + d_r);
% Average pre-stimulus activity
t_sel = (lf.d.time > -0.4) & (lf.d.time < -0.2);
d_x = mean(d_x(:,:,:,t_sel), 4);
% Snip out the alpha band
f_sel = (lf.d.freq >= 8) & (lf.d.freq <= 14);
d_x = mean(d_x(:,:,f_sel), 3);
d_lf = d_x;
clear d_a d_l d_r d_x t_sel f_sel

% Compute a measure of the effect at high frequencies
% L-Hit vs R-Hit at occipital channels
% t = 0
% Collapse across RFT frequencies
vers = 'raw_chans/all';
segment_type = 'target';
hf = load([exp_dir 'tfr/' vers '/' segment_type '/' 'agg']);
% x: Subject x Channel X TaggedFreq x Time
d_x = hf.x;
% Take activity at t = 0
d_x = d_x(:,:,:, hf.d.time == 0);
% Average over the two tagged frequencies
d_x = mean(d_x, 3);
d_hf = d_x;
clear d_x

% Look for correlations between these measures at each channel
r = nan(size(lf.d.label));
p = r;
for i_chan = 1:length(lf.d.label)
    [rho, pval] = corr(d_lf(:,i_chan), d_hf(:,i_chan), ...
        'type', 'Spearman', ...
        'rows', 'complete');
    r(i_chan) = rho;
    p(i_chan) = pval;
end

%% Plot it

i_param = 1;
params = {'r' 'p'};

for i_param = 1:2
    param = params{i_param};
    
    d = [];
    d.label = lf.d.label;
    d.time = 0;

    if strcmp(param, 'r')
        d.avg = r;
        z_lim = 'maxabs';
    elseif strcmp(param, 'p')
        d.avg = p;
        z_lim = [-2 0];
    end

    fn = [exp_dir 'grad/' subject_info.meg{1} '/grad'];
    grad = load(fn);
    d.grad = grad.grad;
    cfg = [];
    cfg.method = 'sum';
    cfg.updatesens = 'yes';
    d = ft_combineplanar(cfg, d);
    % This took the sum -- change to mean
    d.avg = d.avg / 2;
    if strcmp(param, 'p')
        % Take log10 of the p-value
        d.avg = log10(d.avg);
    end
    
    subplot(1, 2, i_param)
    cfg = [];
    cfg.layout = chan.grad_cmb.layout;
    cfg.zlim = z_lim;
    cfg.colorbar = 'no';
    cfg.style = 'straight';
    cfg.comment = 'no';
    cfg.shading = 'interp';
    cfg.markersymbol = '.';
    cfg.gridscale = 200;
    cfg.colorbar = 'yes';
    ft_topoplotER(cfg, d)
    title(param)
end

fn = [exp_dir 'plots/alpha_power/hf_correl'];
print('-dpng', fn)
