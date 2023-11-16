import numpy as np
import pandas as pd
import statsmodels.api as sm
from pathlib import Path
from scipy.stats.distributions import norm
from scipy.stats.distributions import t as tdist
from statsmodels.stats.multitest import local_fdr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pclass = "Pinopsida"
#pclass = "Polypodiopsida"

pdf = PdfPages("%s_largescale.pdf" % pclass)

pa = Path("/home/kshedden/data/Teaching/inaturalist")

fn = pa / ("Plantae_%s.csv.gz" % pclass)

df = pd.read_csv(fn, parse_dates=["eventDate"])

df["day"] = (df["eventDate"] - pd.to_datetime("2010-01-01")).dt.days
df["day"] /= 1000

# Treat longitude as a circular variable
df["lonrad"] = np.pi * df["decimalLongitude"] / 180
df["lonrad_sin"] = np.sin(df["lonrad"])
df["lonrad_cos"] = np.cos(df["lonrad"])

# Group by species and fit a linear model with OLS for each species.
rr = []
for (sp,dx) in df.groupby("species"):

    if dx.shape[0] < 10:
        continue

    md = sm.OLS.from_formula("decimalLatitude ~ day + lonrad_sin + lonrad_cos", data=dx)
    mr = md.fit()
    rr.append([sp, dx.shape[0], mr.params["day"], mr.bse["day"]])

rr = pd.DataFrame(rr, columns=["species", "n", "day_slope", "day_slope_se"])
rr["day_slope_z"] = rr["day_slope"] / rr["day_slope_se"]

# Account for sample size (degrees of freedom in t-distribution)
rr["day_slope_z"] = tdist.cdf(rr["day_slope_z"], rr["n"]-3)
rr["day_slope_z"] = rr["day_slope_z"].clip(1e-8, 1-1e-8)
rr["day_slope_z"] = norm.ppf(rr["day_slope_z"])

# To control family-wise error rates using the Bonferroni approach,
# the Z-scores must exceed this value in magnitude.
bonf_z = norm.ppf(1 - 0.025 / rr.shape[0])

# Local False Discovery Rate (FDR)
rr["locfdr"] = local_fdr(rr["day_slope_z"])

plt.clf()
plt.grid(True)
plt.plot(rr["day_slope"], rr["day_slope_z"], "o", alpha=0.5)
plt.xlabel("Day slope", size=15)
plt.ylabel("Day slope (Z)", size=15)
pdf.savefig()

plt.clf()
plt.grid(True)
plt.plot(rr["day_slope_z"], rr["locfdr"], "o", alpha=0.5)
plt.xlabel("Day slope (Z)", size=15)
plt.ylabel("Local FDR", size=15)
pdf.savefig()

plt.clf()
plt.grid(True)
plt.plot(np.log(rr["n"]), rr["day_slope_z"], "o", alpha=0.5)
plt.xlabel("n", size=15)
plt.ylabel("Day slope (Z)", size=15)
pdf.savefig()

pdf.close()
