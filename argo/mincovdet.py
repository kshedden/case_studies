"""
Estimate the covariance matrix of the temperature and salinity
profiles in Argo using the "minimum determinant covariance"
appproach (MCD).  This estimator is robust to outliers and
can be used to identify unusual observations in large collections
of high-dimensional data.
"""

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib
from sklearn.covariance import MinCovDet
from statsmodels.nonparametric.smoothers_lowess import lowess
#from mpl_toolkits.basemap import Basemap
from read import *

pdf = PdfPages("mincovdet_py.pdf")

def run_mcd(X, label):

    # Run the fast MCD algorithm to estimate the covariance matrix.
    cc = MinCovDet().fit(X)

    # Compare the estimated means
    plt.clf()
    plt.grid(True)
    plt.plot(pressure, cc.location_, label="Robust mean")
    plt.plot(pressure, cc.raw_location_, label="Sample mean")
    plt.xlabel("Pressure", size=15)
    plt.ylabel("Mean %s" % label.lower(), size=15)
    ha, lb =plt.gca().get_legend_handles_labels()
    leg = plt.legend(ha, lb, loc="upper right")
    pdf.savefig()

    # Get the MCD covariance matrix
    cv = cc.covariance_

    # Get the conventional covariance matrix
    c0 = cc.raw_covariance_

    # Compare the spectra of the two covariance matrix estimates
    a,_ = np.linalg.eigh(cv)
    b,_ = np.linalg.eigh(c0)
    plt.clf()
    plt.grid(True)
    plt.plot(np.log(a), np.log(b), "-o")
    plt.xlabel("MCD covariance", size=15)
    plt.ylabel("Conventional covariance", size=15)
    plt.title("Comparison of spectra for %s" % label.lower())
    pdf.savefig()

    # Get the squared distance of each observation to the center of the
    # distribution using the MCD covariance.
    dmcd = cc.mahalanobis(X)

    # Get the squared distance of each observation to the center
    # of the distribution using the conventional covariance.
    Xc = X - cc.location_
    dcov = (Xc * np.linalg.solve(c0, Xc.T).T).sum(1)

    plt.clf()
    plt.grid(True)
    plt.title("Log distances to center for %s" % label.lower())
    plt.plot(np.log(dmcd), np.log(dcov), "o")
    plt.xlabel("Log MCD distance to center", size=15)
    plt.ylabel("Log conventional distance to center", size=15)
    pdf.savefig()

    # Plot "outliers"
    for di in [dmcd, dcov]:
        ii = np.argsort(-di)

        # Plot the outlier curves
        plt.clf()
        plt.grid(True)
        plt.title("%s outliers" % ("MCD" if di is dmcd else "Sample covariance"))
        plt.plot(pressure, cc.location_, label="Robust mean")
        for i in ii[0:15]:
            plt.plot(pressure, X[i, :], color="grey", alpha=0.8)
        ha, lb = plt.gca().get_legend_handles_labels()
        leg = plt.legend(ha, lb, loc="upper right")
        plt.xlabel("Pressure", size=15)
        plt.ylabel(label, size=15)
        pdf.savefig()

        # Plot the locations of the outliers
        plt.clf()
        plt.grid(True)
        plt.plot(lon, lat, "o", color="grey", alpha=0.3, rasterized=True)
        plt.plot(lon[ii[0:15]], lat[ii[0:15]], "o", color="red")
        plt.xlabel("Longitude", size=15)
        plt.ylabel("Latitude", size=15)
        pdf.savefig()

run_mcd(temp.T, "Temperature")
run_mcd(psal.T, "Salinity")

pdf.close()
