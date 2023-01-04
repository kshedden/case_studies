"""
Use ridge regression to predict X-ray flux at a given time
distance into the future, using a block of consecutive values.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from read import *

pdf = PdfPages("autoreg_py.pdf")

df = get_goes(2017)

# Use blocks of size m, and use the first q observations to predict
# the final observation.
m = 1000
q = 200

# The time points of the predictor information relative to the time
# being predicted.
tax = np.arange(-2*m, -2*(m-q))[0:q]
tax = tax / 60

# Make blocks of 'm' consecutive time points with
# approximately 2-second spacing.
tix, flx = make_blocks(df, m, 0)

flx = np.log(1e-8 + flx)

# Create a design matrix and response vector
x = flx[:, 0:q]
y = flx[:, -1]

# Center the data
y -= y.mean()
x -= x.mean(0)

# Regress y on x using ridge regression, with penalty parameter f.
def ridge(x, y, f):
    u, s, vt = np.linalg.svd(x, 0)
    v = vt.T
    g = s / (s**2 + f)
    b = np.dot(v, np.dot(u.T, y) * g)
    return b

# Consider how the regression coefficients look for various values
# of the penalty parameter.
for f in [1, 10, 100, 1000, 10000]:
    b = ridge(x, y, f)
    plt.clf()
    plt.grid(True)
    plt.plot(tax, b, "-")
    plt.ylabel("Coefficient", size=15)
    plt.xlabel("Minutes before current time", size=15)
    plt.title("f=%d" % f)
    pdf.savefig()

pdf.close()
