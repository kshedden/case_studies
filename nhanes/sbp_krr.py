import numpy as np
import statsmodels.api as sm
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.distance import cdist
from scipy.sparse.linalg import eigsh
from read import df

pdf = PdfPages("sbp_krr.pdf")

vx = ["RIAGENDR", "RIDAGEYR", "BMXWT", "BMXHT", "BMXBMI", "BMXLEG",
      "BMXARML", "BMXARMC", "BMXWAIST", "BMXHIP"]
vn = ["BPXSY1"] + vx

dx = df.loc[:, vn].dropna()

# Recode everything as numeric, convert to numpy arrays.
dx["RIAGENDRx"] = dx.RIAGENDR.replace({"F": 1, "M": -1})
vxx = ["RIAGENDRx" if x == "RIAGENDR" else x for x in vx]
X = dx[vxx].values
y = dx["BPXSY1"].values

# Standardize all variables
xmn = X.mean(0)
X -= xmn
xsd = X.std(0)
X /= xsd
ymn = y.mean()
y -= ymn
ysd = y.std()
y /= ysd

# Fit a linear model using OLS for comparison.
u, s, vt = np.linalg.svd(X, 0)
yh = np.dot(u, np.dot(u.T, y))
rmse_ols = np.sqrt(np.mean((y - yh)**2))

# Squared exponential kernel with a given scale parameter.
def ker_sqexp(X, scale, Y=None):
    if Y is None:
        Y = X
    D = cdist(X, Y)
    return np.exp(-D**2 / (2*scale**2))

# Plot the RMSE for various regularization parameters and
# kernel scale parameters.
def goodness_of_fit(ker, lam_range, scale_range, ti):
    plt.clf()
    plt.axes([0.12, 0.12, 0.7, 0.8])
    plt.grid(True)

    # Sweep out a range of kernel parameters.
    for scale in scale_range:

        K = ker(X, scale)
        S, B = np.linalg.eigh(K)
        Bty = np.dot(B.T, y)

        # Calculate the RMSE for models with different
        # degrees of regularization.
        rmse = []
        for lam in lam_range:
            alpha_hat = np.dot(np.dot(B, np.diag(S/(S**2 + lam))), Bty)
            yhat = np.dot(K, alpha_hat)
            e = np.sqrt(np.mean((y - yhat)**2))
            rmse.append(e)
        plt.plot(lam_range, rmse, "-", label="%.1f" % scale)

    plt.title(ti)
    plt.axhline(rmse_ols, label="OLS", color="black")
    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, loc="center right")
    leg.draw_frame(False)
    leg.set_title("Kernel scale")
    plt.xlabel("Regularization", size=15)
    plt.ylabel("RMSE", size=15)
    pdf.savefig()

# Plot the fitted values for females and males, using the given
# kernel, scale, and regularization parameters.
def plot_fit(ker, scale, lam):

    K = ker(X, scale)
    S, B = np.linalg.eigh(K)
    alpha_hat = np.dot(np.dot(B, np.diag(S/(S**2 + lam))), np.dot(B.T, y))

    females = np.flatnonzero(X[:, 0] > 0)
    males = np.flatnonzero(X[:, 0] < 0)

    ages = np.linspace(18, 80, 50)
    ages_std = (ages - xmn[1]) / xsd[1]
    Xp = np.zeros_like(X[0:100, :])
    Xp[0:50, 0] = X[females[0], 0]
    Xp[50:, 0] = X[males[0], 0]
    Xp[0:50, 1] = ages_std
    Xp[50:, 1] = ages_std

    Kp = ker(X, scale, Xp)
    yhat = np.dot(Kp.T, alpha_hat)
    yhat = ymn + ysd * yhat

    plt.clf()
    plt.axes([0.12, 0.12, 0.7, 0.8])
    plt.grid(True)
    plt.plot(ages, yhat[0:50], "-", color="purple", label="Female")
    plt.plot(ages, yhat[50:100], "-", color="orange", label="Male")
    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, loc="center right")
    leg.draw_frame(False)
    plt.title("scale=%.2f, regularization=%.2f" % (scale, lam))
    plt.xlabel("Age", size=15)
    plt.ylabel("SBP", size=15)
    pdf.savefig()

scale_range = [0.5, 1]
lam_range = [0.5, 1, 2, 3]

goodness_of_fit(ker_sqexp, lam_range, scale_range, "Squared exponential kernel")

for scale in scale_range:
    for lam in lam_range:
        plot_fit(ker_sqexp, scale, lam)

pdf.close()
