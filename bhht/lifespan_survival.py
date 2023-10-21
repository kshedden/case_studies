"""
Examine the lifespans of notable people using the BHHT data.

This analysis uses survival analysis methods, allowing us
to use information from still-living people.
"""

import pandas as pd
import numpy as np
import plotille
import statsmodels.api as sm
import os

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"

# Load the dataset.  Use the latin-1 encoding since there is some non-UTF
# data in the file.  Add "nrows=100000" when developing to reduce the run
# time (but use the complete data to get final results).
df = pd.read_csv(os.path.join(pa, "cross-verified-database.csv.gz"), encoding="latin-1", nrows=100000)

# Create a lifespan variable (years of life).  It will be missing for people who are currently living.
df["lifespan"] = df["death"] - df["birth"]

# Exclude people born before 1500, there is too little data to gain a meaningful
# understanding of the trends in lifespan prior to this year.
dx = df.loc[df["birth"] >= 1500, ["birth", "lifespan", "gender", "un_region", "level1_main_occ"]]

# There are a small number of people with missing or "Other" gender but it
# is too small of a sample to draw conclusions.
dx = dx.loc[dx["gender"].isin(["Female", "Male"]), :]

# Drop uninformative occupation codes.
dx = dx.loc[~dx["level1_main_occ"].isin(["Missing", "Other"]), :]

# Censor at 2022
censor_year = 2022
dx["clifespan"] = dx["lifespan"].fillna(censor_year - dx["birth"])
dx["died"] = 1 - 1*dx["lifespan"].isnull()

# Now we can drop all rows with missing data
dx = dx.drop("lifespan", axis=1)
dx = dx.dropna()

# A categorical variable indicating the century in which a person was born.
dx["era"] = np.floor((dx["birth"] - 1500) / 100)

# Plot the survival functions for people born in each century
fig = plotille.Figure()
fig.set_x_limits(0, 100)
fig.set_y_limits(0, 1)
fig.x_label = "Age"
fig.y_label = "Proportion alive"
fig.width = 60
fig.height = 20
for k,g in dx.groupby("era"):
    sf = sm.SurvfuncRight(g.clifespan, g.died)
    la = "%.0f" % (1500 + k*100)
    fig.plot(sf.surv_times, sf.surv_prob, label=la)
print(fig.show(legend=True))

# Fit a proportional hazards regression model
fml = "clifespan ~ bs(birth, 4) + gender + level1_main_occ + un_region"
m0 = sm.PHReg.from_formula(fml, dx, status="died")
r0 = m0.fit()

# Plot the partial effect of birth
dp = dx.iloc[0:100, :].copy()
dp["gender"] = "Female"
dp["level1_main_occ"] = "Leadership"
dp["un_region"] = "Asia"
dp["birth"] = np.linspace(dx["birth"].min(), dx["birth"].max(), 100)
lhr = r0.predict(exog=dp).predicted_values
plt = plotille.plot(dp["birth"], lhr, X_label="Birth year (transformed)",
                    Y_label="Contribtution to the log hazard", width=60, height=20)
print(plt)
print("\n\n")

# Plot the estimated baseline cumulative hazard function
ti, cumhaz, surv = r0.baseline_cumulative_hazard[0]
plt = plotille.plot(ti, cumhaz,
                    X_label="Age", Y_label="Cumulative hazard",
                    width=60, height=20, x_min=0, x_max=100)
print(plt)
print("\n\n")

# Plot the estimated baseline cumulative hazard function on the log scale
ti, cumhaz, surv = r0.baseline_cumulative_hazard[0]
plt = plotille.plot(ti, np.log(np.clip(cumhaz, 1e-4, np.inf)),
                    X_label="Age", Y_label="Log cumulative hazard",
                    width=60, height=20, x_min=0, x_max=100)
print(plt)
print("\n\n")

# Plot the estimated baseline hazard function using numerical differentiation
ti, chaz, surv = r0.baseline_cumulative_hazard[0]
haz = np.diff(chaz) / np.diff(ti)
shaz = sm.nonparametric.lowess(np.log(haz), ti[0:-1])
plt = plotille.plot(shaz[:, 0], shaz[:, 1], X_label="Age", Y_label="log hazard",
                    width=60, height=20, x_min=0, x_max=100)
print(plt)
print("\n\n")

# Fit a sex-stratified proportional hazards regression model
fml = "clifespan ~ bs(birth, 4) + level1_main_occ + un_region"
m1 = sm.PHReg.from_formula(fml, dx, status="died", strata="gender")
r1 = m1.fit()

# Plot the baseline hazard function for each sex
bh = r1.baseline_cumulative_hazard
fig = plotille.Figure()
fig.set_x_limits(0, 100)
fig.x_label = "Age"
fig.y_label = "Log hazard"
fig.width = 60
fig.height = 20
snames = m1.surv.stratum_names
for k in 0,1:
    ti = bh[k][0]
    chaz = bh[k][1]
    haz = np.diff(chaz) / np.diff(ti)
    shaz = sm.nonparametric.lowess(np.log(haz), ti[0:-1])
    fig.plot(shaz[:,0], shaz[:,1], label=snames[k])
print(fig.show(legend=True))
