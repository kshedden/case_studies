"""
A researcher seeks to evaluate two risk factors for epilepsy:
sleep duration (x) and resting heart rate (z).  The number of
seizures in a year (y) was recorded and regression analysis will
be used to understand the conditional mean of y given the values
of x and z.  The data were collected at multiple clinics and
due to the possibility of systematic differences between
clinics that cannot be explained by covariates, GEE will be
used to fit the regression models.

This script illustrates how the statistical power can be assessed
for this type of study using simulation.
"""

import statsmodels.api as sm
import numpy as np
import pandas as pd
from scipy.stats.distributions import norm, gamma, poisson

# Generate a Gaussian n-vector with mean zero and exchangeable correlation icc.
def gen_exch(n, icc):
    return np.sqrt(icc)*np.random.normal() + np.sqrt(1 - icc)*np.random.normal(size=n)

# Generate clustered data with a Poisson or Gamma response variable (y),
# and two covariates (x, y).
#
# ngrp: number of groups
# m: average size of each group
# r: correlation between covariates x and z
# y_icc: exchangeable correlation of the outcome
# x_icc: exchangeable correlation of covariate x
# z_icc: exchangeable correlation of covariate z
# b: coefficients for x and z
# scale: the glm scale parameter
# fam: the family used to generate the data
def gendat(ngrp, m, r, y_icc, x_icc, z_icc, b, scale, gfam):

    x, z, g, e = [], [], [], []

    for i in range(ngrp):
        mm = 1 + np.random.poisson(m) # group size
        x1 = gen_exch(mm, x_icc)
        x.append(x1)
        z1 = gen_exch(mm, z_icc)
        z1 = r*x1 + np.sqrt(1 - r**2)*z1
        z.append(z1)
        e.append(gen_exch(mm, y_icc))
        g.append(i*np.ones(mm))

    x = np.concatenate(x)
    z = np.concatenate(z)
    e = np.concatenate(e)
    g = np.concatenate(g)

    # The marginal mean
    mn = np.exp(b[0]*x + b[1]*z)

    u = norm.cdf(e)

    if gfam == "gamma":
        # Mean is a*b, variance is a*b^2
        b = scale*mn
        a = 1 / scale
        y = gamma(a, scale=b).ppf(u)
    else:
        y = poisson(mn).ppf(u)

    da = pd.DataFrame({"x": x, "z": z, "y": y, "g": g})
    return da

# Estimate the power in one settting using simulation.
def run_power(ngrp, m, r, y_icc, x_icc, z_icc, b, scale, mfam, gfam):

    zs = np.zeros((nrep, 3))
    for i in range(nrep):
        da = gendat(ngrp, m, r, y_icc, x_icc, z_icc, b, scale, gfam)
        r0 = sm.GEE.from_formula("y ~ x + z", "g", family=mfam, data=da)
        m0 = r0.fit()
        zs[i, :] = m0.tvalues

    mn = zs[:, 1].mean(0)
    sd = zs[:, 1].std()
    pw = norm(mn, sd).cdf(-2) + 1 - norm(mn, sd).cdf(2)
    return pw

nrep = 100

# Generate gamma and Poisson data from the same mean structure, analyze it
# with either a gamma or Poisson GEE.
b = np.r_[0.15, 0.4]
for s in 0, 1:
    b1 = s*b

    print(["Null hypothesis is true:", "\nAlternative hypothesis is true:"][int(s)])

    mfam = sm.families.Gamma(link=sm.families.links.log())
    pw = run_power(200, 5, 0.3, 0.2, 0.2, 0.2, b1, 4, mfam, "gamma")
    print("%6.3f gamma data gamma model" % pw)
    pw = run_power(200, 5, 0.3, 0.2, 0.2, 0.2, b1, 4, mfam, "poisson")
    print("%6.3f Poisson data gamma model" % pw)

    mfam = sm.families.Poisson()
    pw = run_power(200, 5, 0.3, 0.2, 0.2, 0.2, b1, 4, mfam, "gamma")
    print("%6.3f gamma data Poisson model" % pw)
    pw = run_power(200, 5, 0.3, 0.2, 0.2, 0.2, b1, 4, mfam, "poisson")
    print("%6.3f Poisson data Poisson model" % pw)
