from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fs = 1000
freqs = [4, 63, 78]
t = np.arange(250) / fs

for freq in freqs:
    x = np.cos(t * 2 * np.pi * freq)
    plt.plot(t, x)

for t_start in np.arange(0, 0.25, 0.05):
    rect = mpatches.Rectangle((t_start, -1), width=0.033, height = 0.1,
                              edgecolor='none')
    plt.gca().add_patch(rect)

plt.show()
