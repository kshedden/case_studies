"""
Examine the lifespans of notable people using the BHHT data.

The main statistical tool here is a form of local regression
known as "LOESS", which is a form of local polynomial regression.
"""

import pandas as pd
import numpy as np
import plotille
import statsmodels.api as sm
from scipy.interpolate import interp1d
import os

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"

# Load the dataset.  Use the latin-1 encoding since there is some non-UTF
# data in the file.  Add "nrows=100000" when developing to reduce the run
# time (but use the complete data to get final results).
df = pd.read_csv(os.path.join(pa, "cross-verified-database.csv.gz"), encoding="latin-1")

# Create a lifespan variable (years of life).
df.loc[:, "lifespan"] = df.loc[:, "death"] - df.loc[:, "birth"]

# Examine lifespans of females and males in relation to year of birth.  To avoid
# censoring, exclude people born after 1920.  Also exclude people born before 1500.
dx = df.loc[(df.birth >= 1500) & (df.birth <= 1920), ["birth", "lifespan", "gender", "level1_main_occ"]]
dx = dx.dropna()

# There are a small number of people with missing or "Other" gender but it
# is too small of a sample to draw conclusions.
dx = dx.loc[dx.gender.isin(["Female", "Male"]), :]

# Drop uninformative occupation codes.
dx = dx.loc[~dx.level1_main_occ.isin(["Missing", "Other"]), :]

# plotille is a package for plotting in the terminal.  Feel free to use
# Matplotlib/PyPlot or Seaborn if you want bitmapped graphs.
fig = plotille.Figure()

# Plot the conditional dispersion for loess compared
# to a linear model.
dd = dx.sort_values(by="birth")
m0 = sm.OLS.from_formula("lifespan ~ birth", data=dd).fit()
ares0 = np.abs(m0.resid.values)
m1 = sm.nonparametric.lowess(dd.lifespan, dd.birth)
ares1 = np.abs(dd.lifespan - m1[:, 1])
ma0 = sm.nonparametric.lowess(ares0, dd.birth)
fig.plot(ma0[:, 0], ma0[:, 1], label="Linear")
ma1 = sm.nonparametric.lowess(ares1, dd.birth)
fig.plot(ma1[:, 0], ma1[:, 1], label="LOESS")
print(fig.show(legend=True))

# Estimate the conditional mean lifespan given year of birth for
# females and males.
fig = plotille.Figure()
for (la, dd) in dx.groupby("gender"):
    print("%s %d" % (la, dd.shape[0]))
    dd = dd.sort_values(by="birth")
    ll = sm.nonparametric.lowess(dd.lifespan, dd.birth)

    # It is sufficient to plot every 1000'th point.
    fig.plot(ll[::1000, 0], ll[::1000, 1], label=la)

print(fig.show(legend=True))

fig = plotille.Figure()

for (la, dd) in dx.groupby("level1_main_occ"):
    print("%s %d" % (la, dd.shape[0]))
    dd = dd.sort_values(by="birth")
    ll = sm.nonparametric.lowess(dd.lifespan, dd.birth)
    fig.plot(ll[::1000, 0], ll[::1000, 1], label=la)

print(fig.show(legend=True))

# Use bootstrapping to get a sense of the uncertainty
nboot = 10
birthyear = np.linspace(1520, 1900, 100)
sd = []
for (la, dd) in dx.groupby("gender"):
    dd = dd.sort_values(by="birth")
    n = dd.shape[0]
    y = np.zeros((len(birthyear), nboot))
    for b in range(nboot):
        ii = np.random.choice(n, n)
        ll = sm.nonparametric.lowess(dd.lifespan.iloc[ii], dd.birth.iloc[ii])
        f = interp1d(ll[:, 0], ll[:, 1])
        y[:, b] = f(birthyear)
    sd.append(y.std(1))
psd = np.sqrt(sd[0]**2 + sd[1]**2)

fig = plotille.Figure()
fig.plot(birthyear, psd)
print(fig.show())
