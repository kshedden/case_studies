import pandas as pd
import numpy as np
from scipy.stats.distributions import chi2
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

# The raw data file constructed by prep.jl should be available at this path.
pa = "/home/kshedden/myscratch/plantnet"

# Load the raw data
df = pd.read_csv(os.path.join(pa, "short.csv.gz"))
df["Date"] = pd.to_datetime(df.Date)

# The data are heavily skewed toward the more recent years, optionally
# restrict the analysis to these years.
firstyear = 2010
df = df.loc[df.Date >= pd.to_datetime("%4d-1-1" % firstyear), :]

# Generate some time variables
df["Date"] = pd.to_datetime(df["Date"])
df["year"] = [x.year for x in df.Date]
df["dayOfYear"] = [x.dayofyear for x in df.Date]

# There are very few records from southern hemisphere
# so remove them.
df = df.loc[df.decimalLatitude >= 0, :]

# Generate basis functions for latitude, longitude, elevation,
# and day within year (seasonality).  Seasonality and longitude
# are circular variables.
def setbasis(df):

    # Basis functions for season
    per = 365
    for k in range(4):
        df["sin_day_%d" % k] = np.sin(2*np.pi*df.dayOfYear/per)
        df["cos_day_%d" % k] = np.cos(2*np.pi*df.dayOfYear/per)
        per /= 2

    # Basis functions for longitude
    per = 180
    for k in range(4):
        df["sin_lon_%d" % k] = np.sin(2*np.pi*df.decimalLongitude/per)
        df["cos_lon_%d" % k] = np.cos(2*np.pi*df.decimalLongitude/per)
        per /= 2

    # Basis functions for latitude.
    x = (df["decimalLatitude"] - 45) / 100
    for k in range(1, 4):
        df["lat%d" % k] = x**k

    # Basis functions for elevation.
    x = (df["elevation"] - 100) / 1000
    for k in range(1, 4):
        df["elv%d" % k] = x**k

	# Basis functions for year.
    x = (df["year"] - 2010) / 100
    for k in range(1, 4):
        df["year%d" % k] = x**k

    return df

# These are the terms that will be used to build a regression model.
dayterms = ["sin_day_%d + cos_day_%d" % (j, j) for j in range(4)]
dayterms = " + ".join(dayterms)
lonterms = ["sin_lon_%d + cos_lon_%d" % (j, j) for j in range(4)]
lonterms = " + ".join(lonterms)
latterms = " + ".join(["lat%d" % k for k in range(1, 4)])
elvterms = " + ".join(["elv%d" % k for k in range(1, 4)])

pdf = PdfPages("bmodel_py.pdf")

# Adjust for everything except the response
def get_covariate_terms(response):
	terms = []
	v = ["day", "decimalLongitude", "decimalLatitude", "elevation"]
	for (j,te) in enumerate([dayterms, lonterms, latterms, elvterms]):
		if response != v[j]:
			terms.append(te)
	return " + ".join(terms)

# Fit a linear model for the response variable (should be one of
# latitude, longitude, or elevation) predicted by other control
# variables, with and without year.  If year is included, it
# is modeled linearly.
def fit_linear(response):
    pva = []
    terms = get_covariate_terms(response)
    for k, dv in df.groupby("scientificName"):
        dv = setbasis(dv)
        m0 = sm.OLS.from_formula("%s ~ %s" % (response, terms), data=dv)
        r0 = m0.fit()
        m1 = sm.OLS.from_formula("%s ~ year + %s" % (response, terms), data=dv)
        r1 = m1.fit()
        lrt = 2*(r1.llf - r0.llf)
        dof = r1.params.size - r0.params.size
        pv = 1 - chi2.cdf(lrt, dof)
        pva.append([pv, r1.params[1]])
    pva = np.asarray(pva)
    return pva

# Fit a nonlinear model for the response variable (should be one of
# latitude, longitude, or elevation) predicted by other control
# variables, with and without year.  If year is included, it
# is modeled nonlinearly with a 5-degree spline.
def fit_nonlin(response):
    pvb = []
    terms = get_covariate_terms(response)
    for k, dv in df.groupby("scientificName"):
        dv = setbasis(dv)
        m0 = sm.OLS.from_formula("%s ~ %s" % (response, terms), data=dv)
        r0 = m0.fit()
        m1 = sm.OLS.from_formula("%s ~ bs(year, 5) + %s" % (response, terms), data=dv)
        r1 = m1.fit()
        lrt = 2*(r1.llf - r0.llf)
        dof = r1.params.size - r0.params.size
        pv = 1 - chi2.cdf(lrt, dof)
        pvb.append([pv, 0])
    pvb = np.asarray(pvb)
    return pvb

def make_plots(pva, plot_slopes, title):
    plt.clf()
    plt.grid(True)
    pvas = np.sort(pva[:, 0])
    n = len(pvas)
    x = np.linspace(1/n, 1-1/n, n)
    plt.plot([0, 1], [0, 1], color="black")
    plt.plot(x, pvas, "o", mfc="none")
    plt.xlabel("Expected p-value", size=15)
    plt.ylabel("Observed p-value", size=15)
    plt.title(title)
    pdf.savefig()

    if plot_slopes:
        plt.clf()
        plt.grid(True)
        plt.plot(pva[:, 0], pva[:, 1], "o", mfc="none")
        plt.xlabel("p-value", size=15)
        plt.ylabel("Year slope", size=15)
        plt.title(title)
        pdf.savefig()

title = {"decimalLatitude": "Latitude", "decimalLongitude": "Longitude", "elevation": "Elevation"}
for v in ["decimalLatitude", "decimalLongitude", "elevation"]:
    pva = fit_linear(v)
    pvb = fit_nonlin(v)
    make_plots(pva, True, title[v])
    make_plots(pvb, False, title[v])

pdf.close()
