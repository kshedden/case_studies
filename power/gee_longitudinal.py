import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats.distributions import norm

"""
A researcher is interested in whether a specific behavioral factor (z)
associates with substance use (y).  She plans a longitudinal study
in which behavior and substance use are measured repeatedly over time
within subjects.  The focus is on the association between behavior
and substance use at the same time point, accounting for confounders.

The study will have a rolling recruitment in which m subjects are
recruited per quarter for 2 years (8 quarters).  The study will have a
total duration of 3 years (12 quarters).  The subjects will be
followed longitudinally until the end of the study (so subjects
recruited in quarter 1 will be followed for 12 quarters but subjects
in quarter 8 will only be followed for 5 quarters).

In each quarter, every subject has probability q of dropping out of
the study.  Once a subject drops out they do not return.

In each quarter, the subject reports their substance use on an
(approximately) continuous scale.  In addition, each subject is
assessed for a behavioral trait every second quarter.  A possible
confounding variable x1 is also measured every quarter.

We wish to assess how well the behavioral trait predicts substance
use.  Our plan is to fit a linear model using GEE with an independence
working correlation model.  GEE is used to account for correlations
among repeated measurements on the same subject.  We will use a "last
observation carried forward" approach with the behavioral measurent
(i.e. when the behavioral measurement is not taken we use the
measurement from the previous quarter).

Our goal is to consider the power for detecting an association between
the behavioral variable and the substance use variable, controlling
for the confounder, and accounting for within-subject correlations.
"""

# Generate a n x m matrix in which each row is an independent sample
# from a m-dimensional normal distribution with exchangeable
# covariance having correlation icc.
def gen_exch(n, m, icc):
    x = np.random.normal(size=(n, m))
    z = np.random.normal(size=n)
    x = np.sqrt(icc)*z[:, None] + np.sqrt(1 - icc)*x
    return x

# Generate data for assessing the power.
#
# es: effect size in SD units
# q: probability of dropping out in each quarter
# rx: the correlation between the measured confounder and behavior
# behav_icc: the exchangeable ICC of the behavior measurements
# subuse_icc: the exchangeable ICC of the substance use measurements
# confound_icc: the exchangeable ICC of the measured confounder
# m: the number of subject recruited per wave
def gendat(es, q, rx, behav_icc, subuse_icc, confound_icc, m):

    idx = 0

    ids = []
    behav = []
    wave = []
    subuse = []
    quarter = []
    xx = []

    # Loop over recruitment waves
    for j in range(8):

        # Number of longitudinal measurements for this wave.
        nmeas = 12 - j

        # Subject ids
        ids.append(np.repeat(np.arange(idx, idx+m), nmeas))
        idx += m

        # Quarter of observation
        quarter1 = np.tile(j + np.arange(nmeas), m)
        quarter.append(quarter1)

        # The behavioral data
        behav1 = gen_exch(m, nmeas, behav_icc)
        for j in range(1, behav1.shape[1], 2):
            # Last observation carried forward.
            behav1[:, j] = behav1[:, j-1]
        behav.append(behav1.flatten())

        # The confounder variable
        u = gen_exch(m, nmeas, confound_icc)
        xx1 = rx*behav1 + np.sqrt(1 - rx**2)*u
        xx.append(xx1.flatten())

        # The substance use variable
        u = gen_exch(m, nmeas, subuse_icc)
        subuse1 = es*behav1 + np.sqrt(1 - es**2)*u
        subuse.append(subuse1.flatten())

    da = pd.DataFrame({"id": np.concatenate(ids), "behav": np.concatenate(behav),
                       "quarter": np.concatenate(quarter), "subuse": np.concatenate(subuse),
                       "x1": np.concatenate(xx)})

    # Account for dropout
    dp = pd.DataFrame({"id": np.arange(8*m),
                       "lastquarter": np.ceil(np.log(np.random.uniform(size=8*m)) / np.log(1-q))})
    da = pd.merge(da, dp, on="id", how="left")
    da = da.loc[da["quarter"] <= da["lastquarter"], :]
    da = da.drop("lastquarter", axis=1)

    return da

nrep = 100

# Estimate the power in a number of scenarios.
def run_power(subuse_icc, behav_icc, rx, m, es):

    rslt = []
    for subuse_icc1 in subuse_icc:
        for behav_icc1 in behav_icc:
            for rx1 in rx:
                for m1 in m:
                    for es1 in es:
                        zs = np.zeros((nrep, 3))
                        for k in range(nrep):
                            da = gendat(es1, 0.02, rx1, behav_icc1, subuse_icc1, 0.1, m1)
                            m0 = sm.GEE.from_formula("subuse ~ behav + x1", groups="id", data=da)
                            r0 = m0.fit()
                            zs[k, :] = r0.params / np.sqrt(np.diag(r0.cov_params()))

                        z = zs[:, 1]
                        zm = z.mean()
                        pw = 1 - norm(zm, z.std()).cdf(2)
                        row = [m1, es1, rx1, behav_icc1, subuse_icc1, pw, zm]
                        rslt.append(row)
    rslt = pd.DataFrame(rslt)
    rslt.columns = ["m", "es", "rx", "behav_icc", "subuse_icc", "power", "zbar"]
    return rslt

# Compare different levels of confounding and different effect sizes
rslt = run_power([0.1], [0.1], [0, 0.5], [10], [0.1, 0.2])
print(rslt)
print("\n")

# Compare different ICC's in the behavioral data
rslt = run_power([0], [0, 0.25, 0.5, 0.75], [0.2], [10], [0.15])
print(rslt)
print("\n")

# Compare different ICC's in the substance use data
rslt = run_power([0, 0.25, 0.5, 0.75], [0], [0.2], [10], [0.15])
print(rslt)
print("\n")

# Compare different ICC's in the behavioral and substance use data
rslt = run_power([0, 0.5], [0, 0.5], [0.2], [10], [0.15])
print(rslt)
print("\n")
