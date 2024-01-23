## Analyzing demographic variation among US counties using biplots

# The data considered here are a single year of population counts for
# US counties.  The population within each county is partitioned into
# 2 x 2 x 4 x 19 = 304 demographic cells (sex x Hispanic ethnicity
# status x race x age).  See the prep.py script for more information.

import numpy as np
import pandas as pd
from prep import demog, births, pop, rucc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Get the FIPS codes that are used in the principal components
# regression (PCR). These are the larger counties.

da = pd.merge(births, pop, on="FIPS", how="left")
da = pd.merge(da, rucc, on="FIPS", how="left")
da["logPop"] = np.log(da["Population"])
da = da.dropna()
fips = da.FIPS.dropna().unique()

# Restrict the demographics data to the larger counties.

demogx = demog.loc[demog.index.isin(fips)]
demogx = np.asarray(demogx)

# The population counts have already been square root transformed.
# Square root transform again to better symmetrize the data.

demogx = np.sqrt(demogx)

# Double center the data.

demogx -= demogx.mean()
demogx -= demogx.mean(0)
demogx -= demogx.mean(1)[:, None]

# Factor the data matrix

u,s,vt = np.linalg.svd(demogx, 0)
v = vt.T

# For biplots, the singular values are partitioned between the left
# and right singular vectors. alpha = 1 gives a distance
# interpretation for rows (counties), alpha = 0 gives a distance
# interpretation for columns (demographic categories), alpha = 0.5
# does not have a strict distance interpretation.

alpha = 0.5
uu = np.dot(u, np.diag(s**alpha))
vv = np.dot(v, np.diag(s**(1-alpha)))

# Remove outlier counties to make the plot easier to read.

for j in range(2):
    qq = np.quantile(uu[:, j], [0.25, 0.5, 0.75])
    iqr = qq[2] - qq[1]
    ii = np.flatnonzero(np.abs(uu[:, j] - qq[1]) < 5*iqr)
    print("Dropping %d outliers in biplot" % (uu.shape[0] - len(ii)))
    uu = uu[ii, :]

# Specify some parameters for plotting.

colors = {"A": "purple", "B": "orange", "N": "lime", "W": "red"}
lt = {"F": "-", "M": ":"}
sym = {"H": "s", "N": "o"}
ages = range(0, 19)

pdf = PdfPages("biplots_py.pdf")

def generate_biplot(uu, vv, c, title):
    """
    Produce a biplot based on the row scores in 'uu' and the column
    scores in vv.  The column labels are in 'c' and the plot is given
    the title 'title'.
    """

    plt.clf()
    plt.figure(figsize=(10, 8))
    ax = plt.axes([0.1, 0.1, 0.76, 0.8])
    ax.grid(True)

    # Plot the counties as grey points
    plt.plot(uu[:, 0], uu[:, 1], 'o', color="grey", alpha=0.3)

    # Plot the demographic categories as colored points, joined
    # by lines connecting the age groups.
    for race in ["A", "B", "N", "W"]:
        for eth in ["N", "H"]:
            la = "%s_%s_%s" % (race, eth, sex)
            ii = [i for (i,x) in enumerate(c) if x.startswith(la)]
            ax.plot(vv[ii,0], vv[ii,1], "-o", color=colors[race], label=la, ms=3)
            ax.text(vv[ii[-1],0], vv[ii[-1],1], eth, ha="left", va="top", color=colors[race])

    ax.set_xlabel("Component 1", size=18)
    ax.set_ylabel("Component 2", size=18)

    ha, lb = ax.get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, loc="center right")
    leg.draw_frame(False)
    ax.set_title(title)

    pdf.savefig()
    plt.show()

# The demographic category labels

c = demog.columns.to_list()
cx = [x.split("_") for x in c]

# To reduce overplotting, produce separate biplots for females and for
# males.

for sex in ["F", "M"]:
    ii = [i for (i,x) in enumerate(cx) if x[2] == sex]
    ii = np.asarray(ii, dtype=int)
    generate_biplot(uu, vv[ii, :], [c[i] for i in ii], "Female" if sex == "F" else "Male")

pdf.close()
