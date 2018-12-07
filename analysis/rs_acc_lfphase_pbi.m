function rs_acc_lfphase_pbi(i_subject)

% Compute PBI on low-frequency phase between hits and misses.
% Phase bifurcation index as in Busch, Dubois, VanRullen (2009)
% PBI_{t,f} = (C_{hits(t,f)} - C_{all(t,f)}) * (C_{misses(t,f)} - C_{all(t,f)})
%   where C is the inter-trial phase coherence for that condition

rs_setup

% Read in the data segmented around targets
fname = subject_info.meg{i_subject};
fn = [exp_dir 'tfr/target/' fname '/low'];
d = load(fn);
d = d.low_freq_data;

% Select only hits and misses
[hits, nans] = rs_resptype(i_subject);

% Exclude the trials with NaNs in the target trialdef
hits = hits(~nans);
hits_and_misses_inx = ismember(hits, [0 1]);

% Keep only the hits and misses (no FAs or late responses)
cfg = [];
cfg.trials = hits_and_misses_inx;
d = ft_selectdata(cfg, d);
hit = hits(hits_and_misses_inx);
clear hits nan

% Get ITPC for hits, misses, and everything together
cfg = [];
cfg.trials = hit == 1;
d_hit = ft_selectdata(cfg, d);
c_hit = itpc(d_hit);
c_hit = c_hit.itpc;

cfg = [];
cfg.trials = hit == 0;
d_miss = ft_selectdata(cfg, d);
c_miss = itpc(d_miss);
c_miss = c_miss.itpc;

c_all = itpc(d);
c_all = c_all.itpc;

% Compute the PBI
pbi = (c_hit - c_all) .* (c_miss - c_all);

label = d.label;
freq = d.freq;
save([exp_dir 'tfr/target/' fname '/low_acc_pbi'], ...
    'pbi', 'label', 'freq')

end


function itc = itpc(freq)

% Compute intertrial phase coherence
% From http://www.fieldtriptoolbox.org/faq/itc/ 0n 2018-12-07
% freq: output of ft_freqanalysis with cfg.output = 'fourier'

% make a new FieldTrip-style data structure containing the ITC
% copy the descriptive fields over from the frequency decomposition

itc = [];
itc.label     = freq.label;
itc.freq      = freq.freq;
itc.time      = freq.time;
itc.dimord    = 'chan_freq_time';

F = freq.fourierspctrm;   % copy the Fourier spectrum
N = size(F,1);           % number of trials

% compute inter-trial phase coherence (itpc) 
itc.itpc      = F./abs(F);         % divide by amplitude  
itc.itpc      = sum(itc.itpc,1);   % sum angles
itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension

% compute inter-trial linear coherence (itlc)
itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
itc.itlc      = abs(itc.itlc);     % take the absolute value, i.e. ignore phase
itc.itlc      = squeeze(itc.itlc); % remove the first singleton dimension
end