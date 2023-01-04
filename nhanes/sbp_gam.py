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

x_spline = dx[["RIDAGEYR", "BMXBMI"]]
bs = BSplines(x_spline, df=[12, 10], degree=[3, 3])

f0 = "BPXSY1 ~ RIDAGEYR + BMXBMI + RIAGENDR"
alpha = np.r_[20., 20.]
m0 = GLMGam.from_formula(f0, data=dx, smoother=bs, alpha=alpha)
r0 = m0.fit()

r0.plot_partial(0, cpr=True)
pdf.savefig()

r0.plot_partial(1, cpr=True)
pdf.savefig()

# Very slow
#alpha = m0.select_penweight(niter=20)[0]

pdf.close()
