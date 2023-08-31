import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import genpareto
import os

# A paper using extreme value techniques to study rainfall in Brazil:
# https://link.springer.com/article/10.1007/s42452-020-03199-8

# Estimate the parameters of a generalized Pareto distribution
# using the empirical Bayes method of Zhang and Stephens.
# https://www.jstor.org/stable/pdf/40586625.pdf
def gp_estimate(z):

    z = np.sort(z)
    n = len(z)
    xstar = z[int(np.round(n/4 + 0.5))]
    m = np.ceil(20 + np.sqrt(n))
    xmax = z.max()

    jj = np.arange(1, m+1)
    tgrid = 1/xmax + (1 - np.sqrt(m/(jj-0.5))) / (3 * xstar)

    def profile(theta):
        k = -np.log(1 - theta*z).mean()
        return n*(np.log(theta/k) + k - 1)

    ltg = np.asarray([profile(t) for t in tgrid])
    ltg -= ltg.max()
    Ltg = np.exp(ltg)
    Ltg /= Ltg.sum()
    theta_hat = np.dot(Ltg, tgrid)
    k_hat = -np.log(1 - theta_hat*z).mean()
    sigma_hat = k_hat / theta_hat

    return -k_hat, sigma_hat

def hill(z, k=200):
    z = np.sort(z)
    z = np.log(z)
    n = len(z)
    return (z[-k+1:] - z[-k]).mean()

def check_gp_estimate(shape, scale, thresh):
    z = genpareto.rvs(shape, scale=scale, size=100000)
    z = z[z > thresh] - thresh
    shape_hill = hill(z)
    shape_eb, scale_eb = gp_estimate(z)
    return shape_eb, scale_eb, shape_hill

1/0



target_dir = "/home/kshedden/data/Teaching/precip"

fname = "USW00094847.csv" # Detroit
#fname = "USW00012839.csv" # Miami

df = pd.read_csv(os.path.join(target_dir, fname + ".gz"), parse_dates=["DATE"])

df = df[["DATE", "PRCP"]].dropna()

# Convert precipitation to millimeters
df["PRCP"] /= 10

# Use this threshold for calculating exceedances
thresh = 5.0

# Median annual maximum
df["year"] = df["DATE"].dt.year
annmax = np.median(df.groupby("year")["PRCP"].agg(np.max))
