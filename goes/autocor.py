"""
Calculate Kendall tau autocorrelations for blocks of consecutive observations
taken from the much longer X-ray flux time series.

The block size can be set below via the variable 'bs'.
"""

import numpy as np
from read import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as stats

pdf = PdfPages("autocor_py.pdf")

df = get_goes(2017)

# Block size, if bs=4000 the overall time of each block
# is around 2 hours and 15 minutes.
bs = 4000

# Make blocks of 'bs' consecutive time points with
# approximately 2-second spacing.
tix, flx = make_blocks(df, 2000, 0)

n, p = flx.shape

# Consider autocorrelation at these time lags
d = 25
dlags = np.arange(1, 200, 10)

# Convert lags to time in minutes
dtime = dlags * 2 / 60

# Calculate these quantiles across blocks of the autocorrelations.
pr = [0.25, 0.5, 0.75]

qa = []
for d in dlags:
    # Get the autocorrelation for each block
    r = np.zeros(n)
    for i in range(flx.shape[0]):
        r[i] = stats.kendalltau(flx[i, 0:p-d], flx[i, d:]).correlation
    qa.append(np.quantile(r[np.isfinite(r)], pr))
qa = np.asarray(qa)

plt.clf()
plt.axes([0.1, 0.1, 0.75, 0.8])
plt.grid(True)
for i, p in enumerate(pr):
    plt.plot(dtime, qa[:, i], label="%.2f" % p)
ha, lb = plt.gca().get_legend_handles_labels()
leg = plt.figlegend(ha, lb, "center right")
leg.draw_frame(False)
plt.xlabel("Time lag (minutes)", size=15)
plt.ylabel("Tau autocorrelation quantile", size=15)
pdf.savefig()

pdf.close()
