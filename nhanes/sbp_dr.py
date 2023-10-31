from statsmodels.regression.dimred import SIR
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from matplotlib.backends.backend_pdf import PdfPages
from read import *

pdf = PdfPages("sbp_dr_py.pdf")

vx = ["RIAGENDR", "RIDAGEYR", "BMXWT", "BMXHT", "BMXBMI", "BMXLEG",
      "BMXARML", "BMXARMC", "BMXWAIST", "BMXHIP"]
vn = ["BPXSY1"] + vx

dx = df.loc[:, vn].dropna()

dx["RIAGENDRx"] = dx.RIAGENDR.replace({"F": 1, "M": -1})

# Mean center all variables
for m in dx.columns:
    if dx[m].dtype == np.float64:
        dx[m] -= dx[m].mean()

y = np.asarray(dx["BPXSY1"])
vz = [x.replace("RIAGENDR", "RIAGENDRx") for x in vx]
X = np.asarray(dx[vz])
m = SIR(y, X)
r = m.fit()
scores = np.dot(X, r.params)

# Stratify on the j'th score, plot the mean of SBP with respect to the
# k'th score.
def plotstrat(j, k):

    dp = pd.DataFrame({"strat": scores[:, j], "x": scores[:, k], "y": dx.BPXSY1})
    dp["strat"] = pd.qcut(scores[:, j], 5)

    plt.clf()
    plt.figure(figsize=(7, 4))
    plt.axes([0.12, 0.12, 0.65, 0.8])
    plt.grid(True)
    for ky, dv in dp.groupby("strat"):
        xx = np.linspace(dv.x.min(), dv.x.max(), 100)
        m = lowess(dv.y, dv.x)
        f = interp1d(m[:, 0], m[:, 1])
        la = "%.2f-%.2f" % (ky.left, ky.right)
        plt.plot(xx, f(xx), "-", label=la)

    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, loc="center right")
    leg.draw_frame(False)
    leg.set_title("Score %d" % (j + 1))

    plt.xlabel("Score %d" % (k + 1), size=15)
    plt.ylabel("SBP (centered)", size=15)
    pdf.savefig()

plotstrat(1, 0)
plotstrat(0, 1)

cols = {"F": "orange", "M": "purple"}

# Plot the DV (SBP) against each score, or plot each score against every covariate.
for j in range(2):
    for x in dx.columns:
        if x == "RIAGENDR":
            continue
        plt.figure(figsize=(7, 5))
        plt.clf()
        plt.grid(True)
        if x != "BPXSY1":
            plt.xlabel(x, size=15)
            plt.ylabel("Score %d" % (j + 1), size=15)
            dp = pd.DataFrame({"x": dx[x], "y": scores[:,j], "sex": dx["RIAGENDR"]})
            for sex in "F", "M":
                ii = dp.sex == sex
                dz = dp.loc[ii, :]
                lw = lowess(dz["y"], dz["x"])
                plt.plot(dz["x"], dz["y"], "o", mfc="none", alpha=0.2, color=cols[sex],
                         label=sex, rasterized=True)
                plt.plot(lw[:, 0], lw[:, 1], "-", color=cols[sex])
        else:
            plt.ylabel(x, size=15)
            plt.xlabel("Score %d" % (j + 1), size=15)
            dp = pd.DataFrame({"y": dx[x], "x": scores[:,j], "sex": dx["RIAGENDR"]})
            for sex in "F", "M":
                ii = dp.sex == sex
                dz = dp.loc[ii, :]
                lw = lowess(dz["y"], dz["x"])
                plt.plot(dz["x"], dz["y"], "o", mfc="none", color=cols[sex],
                         alpha=0.2, label=sex, rasterized=True)
                plt.plot(lw[:, 0], lw[:, 1], "-", color=cols[sex])
        ha, lb = plt.gca().get_legend_handles_labels()
        leg = plt.figlegend(ha, lb, loc="center right")
        leg.draw_frame(False)
        pdf.savefig()

pdf.close()
