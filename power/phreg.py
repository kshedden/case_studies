"""
A researcher plans to study the time duration from when a startup is
founded until it fails due to bankruptcy.  There are two main
covariates of interest in the analysis: leverage (x) and size (z).
The researcher plans to use a proportional hazards regression approach
with the time origin being the date of founding and the event time
being the date of bankrupcty.  Firms that are currently operating are
censored on the date of data collection.

The firms belong to different industry sectors and it is plausible
that different sectors have different bankruptcy risks.  To account
for sector effects, the researcher will stratify the analysis on
sector.

See here for an actual study with similar aims:

https://www.tandfonline.com/doi/full/10.1080/00472778.2020.1750302
"""

import numpy as np
import statsmodels.api as sm
import pandas as pd
from scipy.stats.distributions import norm

# Generate n observations from a censored Weibull population
# with hazard proportional to exp(lhr' * [x, z]); x and z
# are standardized covariates such that cor(x, z) = r.
#
# The hazard function depends on shape; if shape = 0 the
# hazard function is constant, if shape > 0 the hazard
# function is increasing in time and if shape < 0 the hazard
# function is decreasing in time.
def gendat(n, r, shape, lhr):
    x = np.random.normal(size=n)
    z = np.random.normal(size=n)
    z = r*x + np.sqrt(1-r**2)*z

    lp = (lhr[0]*x + lhr[1]*z) / shape
    mn = np.exp(-lp)

    # Event time
    evttime = mn * np.random.weibull(shape, size=n)

    # Observation time
    obstime = np.random.weibull(mn.mean())

    status = 1*(evttime <= obstime)
    time = status*evttime + (1 - status)*obstime

    da = pd.DataFrame({"time": time, "status": status,
                       "obstime": obstime, "evttime": evttime,
                       "x": x, "z": z})
    return da

# Generate a sample of observations from a collection of Weibull
# distributions with different shape parameters, provided in the array
# shape.  For other parameters see the docstring for the gendat
# function.
def gendat_stratified(n, r, shape, lhr):
    da = []
    for (j, s) in enumerate(shape):
        dx = gendat(n, r, s, lhr)
        dx["group"] = j
        da.append(dx)
    da = pd.concat(da, axis=0)
    return da

# Use simulation to estimate the power for a given sample size (n)
# correlation between x and z (r), shape parameters, and log hazard
# ratios.
def runsim(n, r, shape, lhr):
    zs = np.zeros((nrep, 2))
    for i in range(nrep):
        da = gendat_stratified(n, r, shape, lhr)
        m = sm.PHReg.from_formula("time ~ x + z", status="status", strata="group", data=da)
        rr = m.fit()
        zs[i, :] = rr.tvalues
    mn = zs.mean(0)
    sd = zs.std(0)
    pw = norm(mn[0], sd[0]).cdf(-2) + 1 - norm(mn[0], sd[0]).cdf(2)
    return pw

# Estimate power for one or more values of n, r, shapes, lhr.
def runsims(n, r, shapes, lhr):
    rslt = []
    for n1 in n:
        for r1 in r:
            for shape1 in shapes:
                for lhr1 in lhr:
                    pw = runsim(n1, r1, shape1, lhr1)
                    row = [n1, r1, lhr1[0], lhr1[1], pw]
                    rslt.append(row)
    rslt = pd.DataFrame(rslt)
    rslt.columns = ["n", "r", "lhr[0]", "lhr[1]", "power"]
    return rslt

# nrep ~ 100 gives stable results but is slow
nrep = 50

shapes = [0.5, 1, 2, 4]
pw = runsims([150], [0.1, 0.2, 0.3, 0.4], [shapes], [np.r_[0.15, 0.4]])
