import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from read import get_goes

df = get_goes(2017)

# Estimate the tail index.  If P(X > t) ~ 1/x^a, then
# the Hill estimator estimates a.
def hill(x, p=0.001):
    x = np.sort(x)
    n = len(x)
    k = int(p*n)
    m = n - k
    lx = np.log(x[m:])
    lx -= lx[0]
    alpha = 1 / np.mean(lx)
    return alpha

# Make a Pareto plot using the upper p proportion
# of the data in x.
def pareto_plot(x, title, p=0.01):
    n = len(x[0])
    k = int((1 - p) * n)
    q = np.arange(1, n+1) / (n+2)
    q = q[k:]

    plt.clf()
    plt.title(title)
    plt.grid(True)
    for y in x:
        y = np.sort(y)
        y = y[k:]
        plt.plot(np.log(y), np.log(1 - q), color="grey", alpha=0.8)
    plt.xlabel("Observed log quantile", size=15)
    plt.ylabel("Log complementary probability", size=15)
    pdf.savefig()

# A more ad-hoc measure of tail thickness.
def tailratio(x, p=0.99, q=0.75, r=0.5, ref=stats.norm):
    qp = np.quantile(x, p)
    qq = np.quantile(x, q)
    qr = np.quantile(x, r)
    numer = (qp - qr) / (qq - qr)
    denom = (ref.ppf(p) - ref.ppf(r)) / (ref.ppf(q) - ref.ppf(r))
    return numer / denom

pdf = PdfPages("tail_py.pdf")

# Check the Hill estimator using Pareto data
print("Pareto data:")
n = int(1e6)
for b in [1, 2, 3, 4]:
    f1 = stats.pareto.rvs(b, size=n)
    pareto_plot([stats.pareto.rvs(b, size=n) for _ in range(10)], "Pareto data with b=%.1f" % b)
    for p in [1e-5, 1e-4, 1e-3, 1e-2]:
        alpha = hill(f1, p)
        print("%5d    %8.6f %8d %12.2f" % (b, p, int(p*len(f1)), alpha))

# Check the Hill estimator using non-Pareto data with a Pareto tail
print("\nNon-Pareto data with Pareto tail:")
n = int(1e6)
for b in [1, 2, 3, 4]:
    x = []
    for _ in range(10):
        f1 = stats.pareto.rvs(b, size=n)
        ii = np.random.choice(range(n), n//2)
        f1[ii] = -np.log(np.random.uniform(len(ii)))
        x.append(f1)
    pareto_plot(x, "Pareto/exponential mixture with b=%.1f" % b)
    f1 = x[0]
    for p in [1e-5, 1e-4, 1e-3, 1e-2]:
        alpha = hill(f1, p)
        print("%5d    %8.6f %8d %12.2f" % (b, p, int(p*len(f1)), alpha))

# What does the Hill estimator do when the tail is not heavy?
n = int(1e6)
x = [np.random.normal(size=n) for _ in range(10)]
pareto_plot(x, "Gaussian data")
f1 = x[0]
print("\nGaussian (light-tailed) data:")
for p in [1e-5, 1e-4, 1e-3, 1e-2]:
    alpha = hill(f1, p)
    print("%8.6f %8d %12.2f" % (p, int(p*len(f1)), alpha))

# Make Pareto plots of the GOES-flux data and first differences.
pareto_plot([df["Flux1"].values], "GOES Flux-1 data")
pareto_plot([df["Flux2"].values], "GOES Flux-2 data")
pareto_plot([np.diff(df["Flux1"].values)], "GOES Flux-1 data (differenced)")
pareto_plot([np.diff(df["Flux2"].values)], "GOES Flux-2 data (differenced)")

# Estimate tail parameters for the GOES-flux data and first differences.
for d in [False, True]:
    f1 = df["Flux1"].values
    if d:
        print("\nX-ray flux data (differenced):")
        f1 = np.diff(f1)
    else:
        print("\nX-ray flux data:")
    for p in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]:
        alpha = hill(f1, p)
        print("%8.6f %8d %12.2f" % (p, int(p*len(f1)), alpha))

# Make plots of the tail ratios
def get_tailratios(x):
    zr = []
    for r in [0.5, 0.75, 0.9]:
        for f in [0.5, 0.9]:
            q = r + (1 - r)*f
            for g in [0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999]:
                p = q + (1 - q)*g
                tw = tailratio(x, p=p, q=q, r=r, ref=stats.expon)
                zr.append([p, q, r, tw])
    zr = pd.DataFrame(np.asarray(zr), columns=["p", "q", "r", "tw"])
    return zr

def plot_tailratios(zr, title):
    plt.clf()
    plt.axes([0.12, 0.12, 0.68, 0.78])
    plt.grid(True)
    for (k, dv) in zr.groupby(["q", "r"]):
        plt.plot(-np.log(1-dv["p"]), dv["tw"], "-", label="%.2f/%.2f" % tuple(k))
    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, "center right")
    leg.draw_frame(False)
    leg.set_title("q/r")
    plt.xlabel("-log(1-p)", size=15)
    plt.ylabel("Tail ratio", size=15)
    plt.title(title)
    pdf.savefig()

fl1 = df["Flux1"].values
zr = get_tailratios(fl1)
plot_tailratios(zr, "Flux-1 tail ratios")

fl2 = df["Flux2"].values
zr = get_tailratios(fl2)
plot_tailratios(zr, "Flux-2 tail ratios")

fl1d = np.diff(fl1)
zr = get_tailratios(fl1d)
plot_tailratios(zr, "Differenced flux-1 tail ratios")

fl2d = np.diff(fl2)
zr = get_tailratios(fl2d)
plot_tailratios(zr, "Differenced flux-2 tail ratios")

pdf.close()
