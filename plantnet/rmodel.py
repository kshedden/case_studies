import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
from read import *

pdf = PdfPages("rmodel_py.pdf")

# Each of these variables can be the outcome, and the other
# two are predictors.
va = ["decimalLatitude", "decimalLongitude", "elevation"]

# Names for the variables in 'va' used for labeling plots.
titles = {"decimalLatitude": "Latitude", "decimalLongitude": "Longitude", "elevation": "Elevation"}

df = df.dropna()

# Statsmodels MixedLM is slow so use this during demonstration/debugging
# to reduce the sample size.
df = df.sample(10000)

def fitmodel(v):

    # Control for the variables except for the outcome
    te = [x for x in va if x != v]
    tx = ["bs(%s, 5)" % x for x in te]

    s = " + ".join(tx)
    fml = "%s ~ year + %s" % (v, s)

    # Fit the mixed model and print the model summary.
    m0 = sm.MixedLM.from_formula(fml, groups="scientificName", re_formula="1 + decade", data=df)
    r0 = m0.fit()
    print("n=%d observations" % r0.nobs)
    rr = r0.random_effects
    print("n=%d species" % len(rr))
    print(r0.summary())
    print("\n\n")

    return r0, rr

# Generate a plot of predicted trends for each species.
def make_plots(mm, rr, va):

    rx = [np.asarray(x) for x in rr.values()]
    rx = np.asarray(rx)

    # The fitted value at the central location, used as an intercept.
    xx = mm.model.exog
    xm = xx.mean(0)
    xm[1] = 0
    icept = np.dot(xm, mm.fe_params)

    # The year slope (in year units.)
    ys = mm.params[1]

    # Plot these years
    yr = np.linspace(2010, 2020, 20)

    # The plotted years as decades
    yx = (yr - yr.mean()) / 10

    plt.clf()
    plt.grid(True)
    for i in range(rx.shape[0]):
        plt.plot(yr, icept + rx[i, 0] + ys*yr + rx[i, 1]*yx, "-", color="grey", alpha=0.4)
    plt.xlabel("Year", size=15)
    plt.ylabel(titles[va], size=15)
    pdf.savefig()


for v in va:
    mm, rr = fitmodel(v)
    make_plots(mm, rr, v)

pdf.close()
