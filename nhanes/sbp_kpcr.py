import numpy as np
import statsmodels.api as sm
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.distance import cdist
from scipy.sparse.linalg import eigsh
from read import df

pdf = PdfPages("sbp_kpcr.pdf")

vx = ["RIAGENDR", "RIDAGEYR", "BMXWT", "BMXHT", "BMXBMI", "BMXLEG",
      "BMXARML", "BMXARMC", "BMXWAIST", "BMXHIP"]
vn = ["BPXSY1"] + vx

dx = df.loc[:, vn].dropna()

# Recode everything as numeric, convert tpo numpy arrays.
dx["RIAGENDRx"] = dx.RIAGENDR.replace({"F": 1, "M": -1})
vxx = ["RIAGENDRx" if x == "RIAGENDR" else x for x in vx]
x = dx[vxx].values
y = dx["BPXSY1"].values

# Standardize all variables
x -= x.mean(0)
x /= x.std(0)
y -= y.mean(0)

# Fit a linear model using OLS for comparison.
u, s, vt = np.linalg.svd(x, 0)
yh = np.dot(u, np.dot(u.T, y))
rmse_ols = np.sqrt(np.mean((y - yh)**2))

# Use up to this number of basis functions.
nev = 20

# Squared exponential kernel with a given scale parameter.
def ker_sqexp(K, scale):
    return np.exp(-K**2/(2*scale**2))

# Polynomial kernel with exponent m.
def ker_poly(K, m):
    return (1 + np.dot(x, x.T))**m

# Get the basis functions from the data matrix x using the given kernel.
def get_basis(x, ker, par, nev):
    K = cdist(x, x)
    K = ker(K, par)
    ei,ev = eigsh(K, k=nev)

    # Sort by decreasing eigenvalue.
    ii = np.argsort(-ei)
    ei = ei[ii]
    ev = ev[:, ii]
    return ei, ev

def run_kernels(nev, ker, par_range, ti):
    plt.clf()
    plt.axes([0.12, 0.12, 0.75, 0.8])
    plt.grid(True)

    # Sweep out a range of kernel parameters.
    for par in par_range:

        ei, ev = get_basis(x, ker, par, nev)

        rmse = []
        for k in range(1, nev):
            b = ev[:, 0:k]
            f = np.dot(b, np.dot(b.T, y))
            e = np.sqrt(np.mean((y - f)**2))
            rmse.append(e)
        plt.plot(range(1,nev), rmse, "-", label=str(par))

    plt.title(ti)
    plt.axhline(rmse_ols, label="OLS", color="black")
    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, "center right")
    leg.draw_frame(False)
    plt.xlabel("Number of components", size=15)
    plt.ylabel("RMSE", size=15)
    pdf.savefig()

run_kernels(nev, ker_sqexp, [1, 2, 3, 4, 5, 10], "Squared exponential kernel")
run_kernels(nev, ker_poly, [1, 2, 3], "Polynomial kernel")

pdf.close()
