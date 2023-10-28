import numpy as np
import pandas as pd
from statsmodels.gam.api import GLMGam, BSplines
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from matplotlib.backends.backend_pdf import PdfPages
from read import *

pdf = PdfPages("sbp_gam_py.pdf")

vn = ["BPXSY1", "RIAGENDR", "RIDAGEYR", "BMXBMI"]
dx = df.loc[:, vn]
dx = dx.dropna()

# Build spline bases for two quantitative variables to
# include as smooth terms in a GAM predicting SBP.
x_spline = dx[["RIDAGEYR", "BMXBMI"]]
bs = BSplines(x_spline, df=[12, 10], degree=[3, 3])

# Fit a GAM including additive effects for age, BMI, and sex.
f0 = "BPXSY1 ~ RIDAGEYR + BMXBMI + RIAGENDR"
alpha = np.r_[20., 20.]
m0 = GLMGam.from_formula(f0, data=dx, smoother=bs, alpha=alpha)
r0 = m0.fit()

# Create a component plus residual plot of the GAM specied by
# the GAM `res`, for the variable with index `idx`.
def cpr_plot(res, idx):
    from statsmodels.graphics.utils import create_mpl_ax
    mod = res.model
    smoothers = mod.smoother.smoothers
    x = smoothers[idx].x
    y_est, se = res.partial_values(idx)
    ii = np.argsort(x)
    x = x[ii]
    y_est = y_est[ii]
    se = se[ii]
    xname = smoothers[idx].variable_name
    y_est += res.params[xname] * x
    fig, ax = create_mpl_ax(None)
    ax.grid(True)
    ax.plot(x, y_est + res.resid_working.iloc[ii], "o", color="grey", alpha=0.2)
    ax.plot(x, y_est, "-")
    ax.set_xlabel(xname)
    return fig

fig = cpr_plot(r0, 0)
fig.axes[0].set_ylabel("BPXSY1")
pdf.savefig()

fig = cpr_plot(r0, 1)
fig.axes[0].set_ylabel("BPXSY1")
fig.axes[0].set_xlim(10, 50)
pdf.savefig()

# This can be used to get a sense of good tuning parameters.
# It is very slow
#alpha = m0.select_penweight(niter=20)[0]

pdf.close()
