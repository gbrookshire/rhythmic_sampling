""" Preliminary analysis of the QUEST procedures.
"""

import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
# plt.style.use('classic')

exp_dir = 'C:/Users/brookshg/Documents/rhythmic_sampling/data/logfiles/'
n_subj = 8

subject_info = pd.read_csv(exp_dir + '../subject_info.csv')
subj = subject_info['behav'][n_subj-1]
fname_stem = '{subj}_{side}_{freq}.0.psydat'

line_type = {63: '-', 78: '--'}
line_col = {'left': 'green', 'right': 'purple'}
for freq in (63, 78):
    for side in ('left', 'right'):
        fname = fname_stem.format(subj=subj, side=side, freq=freq)
        with open(exp_dir + fname) as f:
            d = pickle.load(f)
            x = d.intensities
            # x = np.delete(x, np.nonzero(np.diff(x) == 0))
            plt.plot(range(len(x)),
                     x,
                     linestyle=line_type[freq],
                     color=line_col[side],
                     label='{} Hz, {}'.format(freq, side))
            plt.plot(len(x), d.quantile(0.5),
                     '<', color=line_col[side])
plt.legend()
plt.ylabel('Target opacity')
plt.xlabel('Trial number')
plt.ylim([0, 0.5])
plt.title('QUEST thresholding')

plt.savefig(exp_dir + '../plots/quest/' + subj + '.png')
plt.show()

#########
# Is this actually working?
# Why are there 40 obs. per cell when there should be 32?
