import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib
from statsmodels.nonparametric.smoothers_lowess import lowess
from mpl_toolkits.basemap import Basemap
from read import *

# See this reference for information about the support points
# algorithm: https://arxiv.org/pdf/1609.01811.pdf

# Equation 22 in Mak et al.
def update_support(X, Y):
    N, p = Y.shape
    n, _ = X.shape
    XX = np.zeros((n, p))

    for i in range(n):
        Dx = X[i, :] - X
        DxN = np.linalg.norm(Dx, axis=1)
        DxN[i] = np.inf
        Dy = X[i, :] - Y
        DyN = np.linalg.norm(Dy, axis=1)
        q = (1/DyN).sum()
        XX[i, :] = np.dot(1/DxN, Dx) * (N / n)
        XX[i, :] += np.dot(1/DyN, Y)
        XX[i, :] /= q

    return XX

# Calculate N support points for the data in Y.  The points
# are stored in the rows of Y.
def support(Y, N, maxiter=1000):

    n, p = Y.shape
    X = np.random.normal(size=(N, p))

    for i in range(maxiter):
        X1 = update_support(X, Y)
        ee = np.linalg.norm(X1 - X)
        X = X1
        if ee < 1e-8:
            break

    return X

pdf = PdfPages("support_py.pdf")

def plot_support_map(ii, title):
    plt.clf()
    plt.figure(figsize=(7, 8))
    plt.axes([0.1, 0.1, 0.8, 0.8])
    m = Basemap(llcrnrlon=-100.,llcrnrlat=-65.,urcrnrlon=30.,urcrnrlat=80.,
                resolution='l',projection='merc', lat_0=0.,lon_0=0.)

    for j in range(ii.max() + 1):
        jj = np.flatnonzero(ii == j)
        x, y = m(lon[jj], lat[jj])
        plt.scatter(x, y, s=8, label=str(1+j))

    ha,lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, "center right")
    leg.draw_frame(False)

    m.drawcoastlines()
    m.drawmapboundary()
    plt.title(title)
    pdf.savefig()

# Find the position of the closest support point
# in the rows of S to the vectors in the rows
# of X.
def support_neighbor(X, S):
    ii = np.zeros(X.shape[1]).astype(int)
    for i in range(X.shape[1]):
        d = ((X[:, i] - S)**2).sum(1)
        ii[i] = np.argmin(d)
    return ii

# Plot support points for temperature and salinity separately.
for (j,x) in enumerate([temp, psal]):

    # Make plots with different numbers of support points.
    for npt in 5, 10, 20:
        print("npt=", npt)
        X = support(x.T, npt, maxiter=200)
        plt.clf()
        plt.figure(figsize=(6.4,4.8))
        plt.grid(True)
        plt.title("%d support points" % npt)
        for i in range(npt):
            plt.plot(pressure, X[i, :], "-", color="grey")
        plt.xlabel("Pressure", size=15)
        plt.ylabel(["Temperature", "Salinity"][j], size=15)
        pdf.savefig()

# Plot support points for the combined temperature and salinity
# trajectories.  Normalize the ranges of temperature and salinity
# so that the support points are more equally based on the two
# variables.
cm = matplotlib.cm.get_cmap("tab10")
tempz = (temp - temp.mean()) / temp.std()
psalz = (psal - psal.mean()) / psal.std()
pt = np.vstack([tempz, psalz])
for npt in 5,:
    print("npt=", npt)
    X = support(pt.T, npt, maxiter=200)
    plt.clf()
    plt.figure(figsize=(6.4,4.8))
    plt.axes([0.1, 0.1, 0.78, 0.8])
    plt.grid(True)
    plt.title("%d support points" % npt)
    ax1 = plt.gca()
    for i in range(npt):
        ax1.plot(pressure, X[i, 0:100], "-", color=cm(i/10))
    ax1.set_ylabel("Temperature (solid lines)", size=15)
    ax2 = ax1.twinx()
    for i in range(npt):
        ax2.plot(pressure, X[i, 100:200], ":", color=cm(i/10))
    ax2.set_ylabel("Salinity (broken lines)", size=15)
    ax1.set_xlabel("Pressure", size=15)
    pdf.savefig()

# Make maps showing the distribution of points falling closest
# to each support point.
npt = 10
for j,x in enumerate([temp, psal]):
    S = support(x.T, npt, maxiter=200)
    ii = support_neighbor(x, S)
    plot_support_map(ii, ["Temperature", "Salinity"][j])

pdf.close()
