import lmom
import pandas as pd
import numpy as np
from read import get_goes
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages("lmoments_py.pdf")

df = get_goes(2017)

# Calculate the raw L-moments
lm = []
for (k, dv) in df.groupby(["Year", "Month", "Day"]):
    v = np.sort(dv["Flux1"].values)
    row = [lmom.l1(v), lmom.l2(v), lmom.l3(v), lmom.l4(v)]
    lm.append(row)

lm = np.asarray(lm)
lm = pd.DataFrame(lm, columns=["l1", "l2", "l3", "l4"])

# Standardized L-moments.
lm["l3s"] = lm["l3"] / lm["l2"]
lm["l4s"] = lm["l4"] / lm["l2"]

v = ["l1", "l2", "l3s", "l4s"]
na = ["L-mean", "L-dispersion", "Standardized L-skew", "Standardized L-kurtosis"]
for j in range(4):
    for k in range(j):
        for dolog in [False, True]:

            # Don't log the skew
            if dolog and 3 in [j, k]:
                continue

            plt.figure(figsize=(6.5, 4.5))
            plt.clf()
            plt.grid(True)
            if dolog:
                plt.plot(np.log(lm[v[j]]), np.log(lm[v[k]]), "o", alpha=0.5, mfc="none")
                plt.xlabel("log %s" % na[j], size=15)
                plt.ylabel("log %s" % na[k], size=15)
            else:
                plt.plot(lm[v[j]], lm[v[k]], "o", alpha=0.5, mfc="none")
                plt.xlabel(na[j], size=15)
                plt.ylabel(na[k], size=15)
            pdf.savefig()

pdf.close()
