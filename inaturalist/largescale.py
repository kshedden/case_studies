import numpy as np
import pandas as pd
import statsmodels.api as sm
from pathlib import Path
from scipy.stats.distributions import norm, chi2
from scipy.stats.distributions import t as tdist
from statsmodels.stats.multitest import local_fdr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#pclass = "Pinopsida"
pclass = "Polypodiopsida"

pdf = PdfPages("%s_largescale.pdf" % pclass)

pa = Path("/home/kshedden/data/Teaching/inaturalist")

fn = pa / ("Plantae_%s.csv.gz" % pclass)

df = pd.read_csv(fn, parse_dates=["eventDate"])

df["day"] = (df["eventDate"] - pd.to_datetime("2010-01-01")).dt.days
df["day"] /= 1000

# Mean latitude per species
meanLat = df.groupby("species")["decimalLatitude"].aggregate(np.mean)
meanLat = meanLat.rename("meanLatitude")

# Treat longitude as a circular variable
df["lonrad"] = np.pi * df["decimalLongitude"] / 180
df["lonrad_sin"] = np.sin(df["lonrad"])
df["lonrad_cos"] = np.cos(df["lonrad"])

# Create a variable that cannot contain any information about the outcome.
df["fake"] = df["lonrad_cos"] + np.random.normal(size=df.shape[0])

# Group by species and fit a linear model predicting latitude with OLS for each species.
rr = []
for (sp,dx) in df.groupby("species"):

    if dx.shape[0] < 10:
        continue

    md = sm.OLS.from_formula("decimalLatitude ~ day + lonrad_sin + lonrad_cos + fake", data=dx)
    mr = md.fit()

    md2 = sm.OLS.from_formula("decimalLatitude ~ day * (lonrad_sin + lonrad_cos + fake)", data=dx)
    mr2 = md2.fit()

    lrt = 2 * (mr2.llf - mr.llf)
    dof = mr2.df_model - mr.df_model

    if dof == 3 and np.linalg.svd(md.exog,0)[1].min() > 1e-5:
        lrt_z = (lrt / dof)**(1/3)
        lrt_z -= 1 - 2/(9*dof)
        lrt_z /= np.sqrt(2/(9*dof))
    else:
        lrt_z = 0

    rr.append([sp, dx.shape[0], mr.params["day"], mr.bse["day"], mr.params["fake"], mr.bse["fake"], lrt_z])

# Construct Z-scores for parameters of interest
rr = pd.DataFrame(rr, columns=["species", "n", "day_slope", "day_slope_se", "fake_slope", "fake_slope_se", "lrt_z"])
rr["day_slope_z"] = rr["day_slope"] / rr["day_slope_se"]
rr["fake_slope_z"] = rr["fake_slope"] / rr["fake_slope_se"]

rr = pd.merge(rr, meanLat, left_on="species", right_index=True)

# Account for finite group sizes
def t_adjust(rr, vn):
    x = tdist.cdf(rr["%s_z" % vn], rr["n"] - 5)
    x = np.clip(x, 1e-12, 1-1e-12)
    return norm.ppf(x)

rr["day_slope_z"] = t_adjust(rr, "day_slope")
rr["fake_slope_z"] = t_adjust(rr, "fake_slope")

for vn in ["day_slope", "fake_slope", "lrt"]:
    plt.clf()
    plt.grid(True)
    plt.hist(rr["%s_z" % vn], bins=40, density=True)
    xx = np.linspace(-3, 3, 100)
    yy = np.exp(-xx**2/2) / np.sqrt(2*np.pi)
    plt.plot(xx, yy, "-")
    plt.xlabel("%s slope" % vn.title(), size=15)
    plt.ylabel("Frequency", size=15)
    pdf.savefig()

    n = rr.shape[0]
    zs = np.sort(rr["%s_z" % vn])
    xx = np.linspace(1/n, 1-1/n, n)
    yy = norm.ppf(xx)
    plt.clf()
    plt.grid(True)
    plt.plot(zs, yy, "-")
    plt.axline((0, 0), slope=1, color="grey")
    plt.xlabel("Observed %s quantiles" % vn, size=15)
    plt.ylabel("Normal quantiles", size=15)
    pdf.savefig()

# To control family-wise error rates using the Bonferroni approach,
# the Z-scores must exceed this value in magnitude.
bonf_z = norm.ppf(1 - 0.025 / rr.shape[0])

# Local False Discovery Rate (FDR)
rr["locfdr"] = local_fdr(rr["day_slope_z"])

# Plot the day slope Z-score against the mean latitude, to assess whether
# there are systematic trends relative to the equator.
plt.clf()
plt.grid(True)
plt.plot(rr["meanLatitude"], rr["day_slope_z"], "o", alpha=0.5)
plt.xlabel("Mean latitude", size=15)
plt.ylabel("Day slope (Z)", size=15)
pdf.savefig()

# Plot the local FDR against the day slope Z-score.
plt.clf()
plt.grid(True)
plt.plot(rr["day_slope_z"], rr["locfdr"], "o", alpha=0.5)
plt.xlabel("Day slope (Z)", size=15)
plt.ylabel("Local FDR", size=15)
pdf.savefig()

# Plot the day slope Z-score against the sample size.  If we are
# mainly limited by power then the larger Z-scores will be concentrated
# where the sample size is larger.
plt.clf()
plt.grid(True)
plt.plot(np.log(rr["n"]), rr["day_slope_z"], "o", alpha=0.5)
plt.xlabel("Log n", size=15)
plt.ylabel("Day slope (Z)", size=15)
pdf.savefig()

pdf.close()
