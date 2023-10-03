import numpy as np
import pandas as pd
from prep import demog, births, pop, rucc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Get the FIPS codes that are used in the PCR
da = pd.merge(births, pop, on="FIPS", how="left")
da = pd.merge(da, rucc, on="FIPS", how="left")
da["logPop"] = np.log(da["Population"])
da = da.dropna()
fips = da.FIPS.dropna().unique()

# Restrict the demographics data to the larger counties
demogx = demog.loc[demog.index.isin(fips)]

demogx -= demogx.mean()
u,s,vt = np.linalg.svd(demogx, 0)
v = vt.T

#alpha = 1 # distance interpretation for rows (counties)
alpha = 0 # distance interpretation for columns (demographic categories)
uu = np.dot(u, np.diag(s**alpha))
vv = np.dot(v, np.diag(s**(1-alpha)))

for j in range(2):
    qq = np.quantile(uu[:, j], [0.25, 0.5, 0.75])
    iqr = qq[2] - qq[1]
    ii = np.flatnonzero(np.abs(uu[:, j] - qq[1]) < 5*iqr)
    print("Dropping %d outliers in biplot" % (uu.shape[0] - len(ii)))
    uu = uu[ii, :]

colors = {"A": "purple", "B": "orange", "N": "lime", "W": "red"}
lt = {"F": "-", "M": ":"}
sym = {"H": "s", "N": "o"}
ages = range(0, 19)

pdf = PdfPages("biplots_py.pdf")

plt.clf()
plt.figure(figsize=(10, 8))
ax = plt.axes([0.1, 0.1, 0.76, 0.8])
ax.grid(True)
plt.plot(uu[:, 0], uu[:, 1], 'o', color="grey", alpha=0.3)

c = demog.columns.to_list()
for race in ["A", "B", "N", "W"]:
    for eth in ["H", "N"]:
        for sex in ["F", "M"]:
            if sex == "M":
                continue
            la = "%s_%s_%s" % (race, eth, sex)
            ii = [i for (i,x) in enumerate(c) if x.startswith(la)]
            ax.plot(v[ii,0], v[ii,1], "-o", color=colors[race], label=la, ms=3)
            ax.text(v[ii[-1],0], v[ii[-1],1], eth, ha="left", va="top", color=colors[race])

ax.set_xlabel("Component 1", size=18)
ax.set_ylabel("Component 2", size=18)

ha, lb = ax.get_legend_handles_labels()
leg = plt.figlegend(ha, lb, loc="center right")
leg.draw_frame(False)

pdf.savefig()

pdf.close()

