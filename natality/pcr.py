# Examine factors associated with birth count variation among US
# counties using Principal Components Regression and Poisson GLM/GEE.

import pandas as pd
import numpy as np
from prep import births, demog, pop, na, age_groups, rucc
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Create a dataframe for modeling.  Merge the birth data with
# population and RUCC data.
da = pd.merge(births, pop, on="FIPS", how="left")
da = pd.merge(da, rucc, on="FIPS", how="left")

da["logPop"] = np.log(da["Population"])
da = da.dropna()
da = da.sort_values(["FIPS", "year"])

da["yearc"] = da["year"] - da["year"].mean()
da["logPopc"] = da["logPop"] - da["logPop"].mean()
da["RUCC_2013c"] = da["RUCC_2013"] - da["RUCC_2013"].mean()

pdf = PdfPages("pcr_py.pdf")

# Calculate the mean and variance within each county to
# assess the mean/variance relationship.
mv = births.groupby("FIPS")["Births"].agg([np.mean, np.var])
lmv = np.log(mv)

# Regress log variance on log mean
mr = sm.OLS.from_formula("var ~ mean", lmv).fit()
print(mr.summary())

# Plot the log variance against the log mean.  If variance = phi*mean,
# then log(variance) = log(phi) + log(mean), i.e. the slope is 1 and
# the intercept is log(phi).
plt.clf()
plt.grid(True)
plt.plot(lmv["mean"], lmv["var"], "o", alpha=0.2, rasterized=True)
plt.axline((8, 8), slope=1, color="grey")
plt.axline((lmv["mean"].mean(0), lmv["var"].mean(0)), slope=1, color="purple")
plt.axline((8, mr.params[0]+8*mr.params[1]), slope=mr.params[1], color="orange")
plt.xlabel("Log mean", size=16)
plt.ylabel("Log variance", size=16)
pdf.savefig()

# GLM, not appropriate since we have repeated measures on counties
fml = "Births ~ logPop + RUCC_2013"
m0 = sm.GLM.from_formula(fml, family=sm.families.Poisson(), data=da)
r0 = m0.fit() # Poisson
r0x = m0.fit(scale="X2") # Quasi-Poisson

# GEE accounts for the correlated data
m1 = sm.GEE.from_formula(fml, groups="FIPS", family=sm.families.Poisson(), data=da)
r1 = m1.fit() # Poisson and quasi-Poisson are the same for GEE
r1x = m1.fit(scale="X2")

# Use log population as an offset instead of a covariate
m2 = sm.GEE.from_formula("Births ~ RUCC_2013", groups="FIPS", offset="logPop",
                         family=sm.families.Poisson(), data=da)
r2 = m2.fit(scale="X2")

# A diagnostic plot for the variance structure that does not require
# there to be replicates.
plt.clf()
plt.grid(True)
lfv = np.log(r2.fittedvalues).values
apr = np.abs(r2.resid_pearson)
ii = np.argsort(lfv)
lfv = lfv[ii]
apr = apr[ii]
ff = sm.nonparametric.lowess(apr, lfv)
plt.plot(lfv, apr, "o", alpha=0.2, rasterized=True)
plt.plot(ff[:, 0], ff[:, 1], "-", color="orange")
plt.title("Poisson mean/variance model")
plt.xlabel("Log predicted mean", size=16)
plt.ylabel("Absolute Pearson residual", size=16)
pdf.savefig()

# Use Gamma family to better match the mean/variance relationship.
m3 = sm.GEE.from_formula("Births ~ RUCC_2013", groups="FIPS", offset="logPop",
                         family=sm.families.Gamma(link=sm.families.links.log()), data=da)
r3 = m3.fit(scale="X2")

# Diagnostic plot for mean/variance relationship with gamma model.
plt.clf()
plt.grid(True)
lfv = np.log(r3.fittedvalues).values
apr = np.abs(r3.resid_pearson)
ii = np.argsort(lfv)
lfv = lfv[ii]
apr = apr[ii]
ff = sm.nonparametric.lowess(apr, lfv)
plt.plot(lfv, apr, "o", alpha=0.2, rasterized=True)
plt.plot(ff[:, 0], ff[:, 1], "-", color="orange")
plt.title("Gamma mean/variance model")
plt.xlabel("Log predicted mean", size=16)
plt.ylabel("Absolute Pearson residual", size=16)
pdf.savefig()

# Use exchangeable correlation structure.  Since RUCC is constant within
# groups the parameter estimates and standard errors are the same
# as with the independence model.
m4 = sm.GEE.from_formula("Births ~ RUCC_2013", groups="FIPS", offset="logPop",
                         cov_struct=sm.cov_struct.Exchangeable(),
                         family=sm.families.Gamma(link=sm.families.links.log()), data=da)
r4 = m4.fit(scale="X2")

m5 = sm.GEE.from_formula("Births ~ RUCC_2013 + year", groups="FIPS", offset="logPop",
                         cov_struct=sm.cov_struct.Exchangeable(),
                         family=sm.families.Gamma(link=sm.families.links.log()), data=da)
r5 = m5.fit(scale="X2")

m6 = sm.GEE.from_formula("Births ~ RUCC_2013c * yearc", groups="FIPS", offset="logPop",
                         cov_struct=sm.cov_struct.Exchangeable(),
                         family=sm.families.Gamma(link=sm.families.links.log()), data=da)
r6 = m6.fit(scale="X2")

print("Score tests:")
print(r5.model.compare_score_test(r4))
print(r6.model.compare_score_test(r5))

# Get factors (principal components) from the demographic data
demog -= demog.mean(0)
u, s, vt = np.linalg.svd(demog)
v = vt.T

# Convert the coefficients back to the original coordinates
def convert_coef(c, npc):
    return np.dot(v[:, 0:npc], c/s[0:npc])

# The proportion of explained variance.
pve = s**2
pve /= sum(pve)

plt.clf()
plt.grid(True)
plt.plot(range(1, len(pve)+1), pve, "-")
plt.xlabel("Component", size=18)
plt.ylabel("Proportion of explained variance")
pdf.savefig()

plt.clf()
plt.grid(True)
plt.plot(np.log(range(1, len(pve)+1)), np.log(pve), "-")
plt.xlabel("Log component", size=18)
plt.ylabel("Log proportion of explained variance")
pdf.savefig()

plt.clf()
plt.grid(True)
plt.plot(range(1, len(pve)+1), np.cumsum(pve), "-")
plt.xlabel("Component", size=18)
plt.ylabel("Cumulative roportion of explained variance")
pdf.savefig()

# Put the demographic factors into a dataframe
m = {("pc%02d" % k) : u[:, k] for k in range(100)}
m["FIPS"] = demog.index
demog_f = pd.DataFrame(m)

# Merge demographic information into the births data
da = pd.merge(da, demog_f, on="FIPS", how="left")

# Include this number of factors in subsequent models
npc = 10

# A GLM, not appropriate since we have repeated measures on counties
fml = "Births ~ (logPopc + RUCC_2013c) * yearc + " + " + ".join(["pc%02d" % j for j in range(npc)])
m7 = sm.GLM.from_formula(fml, family=sm.families.Poisson(), data=da)
r7 = m7.fit(scale="X2")

# GEE accounts for the correlated data
m8 = sm.GEE.from_formula(fml, groups="FIPS",
         family=sm.families.Gamma(link=sm.families.links.log()), data=da)
r8 = m8.fit(scale="X2")

# Use log population as an offset instead of a covariate
fml = "Births ~ " + " + ".join(["pc%02d" % j for j in range(npc)])
m9 = sm.GEE.from_formula(fml, groups="FIPS", offset="logPop",
         family=sm.families.Gamma(link=sm.families.links.log()), data=da)
r9 = m9.fit(scale="X2")

# Restructure the coefficients so that the age bands are
# in the columns.
def restructure(c):
    ii = pd.MultiIndex.from_tuples(na)
    c = pd.Series(c, index=ii)
    c = c.unstack()
    return c

# This function fits a Poisson GLM to the data using 'npc' principal components
# as explanatory variables.
def fitmodel(npc):
    # A GEE using log population as an offset
    fml = "Births ~ 1" if npc == 0 else "Births ~ RUCC_2013*year + " + " + ".join(["pc%02d" % j for j in range(npc)])
    m = sm.GEE.from_formula(fml, groups="FIPS", family=sm.families.Gamma(link=sm.families.links.log()),
                            offset=da["logPop"], data=da)
    r = m.fit(scale="X2")

    # Convert the coefficients back to the original coordinates
    c = convert_coef(r.params[4:], npc)

    # Restructure the coefficients so that the age bands are
    # in the columns.
    c = restructure(c)

    return c, m, r

# Plot styling information
colors = {"A": "purple", "B": "orange", "N": "lime", "W": "red"}
lt = {"F": "-", "M": ":"}
sym = {"H": "s", "N": "o"}
ages = range(0, 19)

# Fit models with these numbers of PCs.
pcs = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

models = []
for npc in pcs:

    c, m, r = fitmodel(npc)
    models.append((m, r))

    plt.clf()
    plt.figure(figsize=(9, 7))
    ax = plt.axes([0.14, 0.18, 0.7, 0.75])
    ax.grid(True)
    for i in range(c.shape[0]):
        a = c.index[i]
        la = "/".join(a)
        ax.plot(ages, c.iloc[i, :], lt[a[2]] + sym[a[1]], color=colors[a[0]],
                label=la)

    # Setup the horizontal axis labels
    ax.set_xticks(ages)
    ax.set_xticklabels(age_groups)
    for x in plt.gca().get_xticklabels():
        x.set_rotation(-90)

    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, loc="center right")
    leg.draw_frame(False)

    plt.xlabel("Age group", size=17)
    plt.ylabel("Coefficient", size=17)
    plt.title("%d factors" % npc)
    pdf.savefig()

pdf.close()

# Use score tests to get a sense of the number of PC factors
# to include; also consider the PVEs calculated above.
print("Score tests for PC's:")
for k in range(10):
    st = models[k+1][0].compare_score_test(models[k][1])
    print("%d versus %d: p=%f" % (pcs[k+1], pcs[k], st["p-value"]))
