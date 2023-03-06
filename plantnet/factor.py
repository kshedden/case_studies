import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
import os

# Adjust the path below to find the data files.
pa = "/home/kshedden/mynfs/data/Teaching/plantnet"
df = pd.read_csv(os.path.join(pa, "plants_occurrences.csv.gz"))
dz = pd.read_csv(os.path.join(pa, "plants_locations.csv.gz"))

pdf = PdfPages("factor_py_plots.pdf")

# Pivot the data so that the species are in the columns
# and the dates are in the rows.
dx = df[["Date", "scientificName", "nobs"]]
dx = dx.set_index(["Date", "scientificName"]).unstack()
dx.columns = [x[1] for x in dx.columns]

# Variance stabilizing transformation
dx = np.sqrt(dx)

speciesmeans = dx.mean(0)
datemeans = dx.mean(1)

# Double center the data
dx -= dx.mean().mean()
dx -= dx.mean(0)
dx -= dx.mean(1)[:, None]

# Plot date means
plt.clf()
plt.axes([0.1, 0.2, 0.8, 0.7])
plt.grid(True)
plt.plot(dx.index, datemeans, "-")
plt.gca().xaxis.set_major_locator(mdates.YearLocator(5))
for x in plt.gca().xaxis.get_ticklabels():
    x.set_rotation(-90)
plt.xlabel("Date", size=15)
plt.ylabel("Date mean (occurrences/species)", size=15)
pdf.savefig()

# Plot order statistics of species means
plt.clf()
plt.axes([0.15, 0.2, 0.8, 0.7])
plt.grid(True)
plt.plot(np.arange(len(speciesmeans)), np.sort(speciesmeans), "-")
plt.xlabel("Position", size=15)
plt.ylabel("Species mean (occurrences/day)", size=15)
pdf.savefig()

# Make sure the count data and location data are aligned
# by species.
assert((dx.columns == dz.scientificName).all())

# Factor the matrix once, then sort so that the
# species scores are increasing for the first
# factor.
u, s, vt = np.linalg.svd(dx, 0)
v = vt.T
ii = np.argsort(v[:, 0])
dx = dx.iloc[:, ii]
dz = dz.iloc[ii, :]

# Factor the matrix again.
u, s, vt = np.linalg.svd(dx, 0)
v = vt.T

# Scree plot
plt.clf()
plt.grid(True)
plt.plot(s, "-")
plt.xlabel("SVD component", size=15)
plt.ylabel("Singular value", size=15)
pdf.savefig()

# Log/log scree plot
s1 = s[s > 1e-8]
plt.clf()
plt.grid(True)
plt.plot(np.log(np.arange(1, len(s1) + 1)), np.log(s1), "-")
plt.xlabel("SVD component", size=15)
plt.ylabel("Singular value", size=15)
pdf.savefig()

for j in range(10):

    # Plot the date factor scores (only for most recent dates)
    plt.clf()
    plt.axes([0.13, 0.2, 0.8, 0.7])
    plt.grid(True)
    plt.plot(dx.index[-2500:], u[-2500:, j], "-")
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
    for x in plt.gca().xaxis.get_ticklabels():
        x.set_rotation(-90)
    plt.xlabel("Date", size=15)
    plt.ylabel("Date %d factor score" % (j + 1), size=15)
    plt.title("Factor %d" % (j + 1))
    pdf.savefig()

    # Plot the species factor scores
    plt.clf()
    plt.grid(True)
    plt.title("Factor %d" % (j + 1))
    plt.xlabel("Species", size=15)
    plt.ylabel("Species factor %d" % (j + 1), size=15)
    plt.plot(v[:, j], "-")
    pdf.savefig()

    # Plot the species factor scores against the mean latitude for
    # the species.
    plt.clf()
    plt.grid(True)
    plt.plot(dz["decimalLatitude"], v[:, j], "o", mfc="none")
    plt.xlabel("Mean latitude", size=15)
    plt.ylabel("Factor %d score" % (j + 1), size=15)
    pdf.savefig()

    # Plot the species factor scores against the mean longitude for
    # the species.
    plt.clf()
    plt.grid(True)
    plt.plot(dz["decimalLongitude"], v[:, j], "o", mfc="none")
    plt.xlabel("Mean longitude", size=15)
    plt.ylabel("Factor %d score" % (j + 1), size=15)
    pdf.savefig()

    # Plot the species factor scores against the mean elevation for
    # the species.
    plt.clf()
    plt.grid(True)
    plt.plot(dz["elevation"], v[:, j], "o", mfc="none")
    plt.xlabel("Mean elevation", size=15)
    plt.ylabel("Factor %d score" % (j + 1), size=15)
    pdf.savefig()

pdf.close()
