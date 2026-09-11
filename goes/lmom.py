import numpy as np
import warnings

def pwm_coef(n, r):
    """
    Returns the coefficients for estimating the r^th probability weighted
    moment of a sample of size n.  The probability weighted moment estimate
    is a linear function of the order statistics with these coefficients.
    """
    J = np.arange(1, n+1)
    lognumer = np.zeros(n)
    for k in range(r):
        lognumer[r:] += np.log(J[r:] - k - 1)
    logcoef = lognumer - np.sum(np.log(np.arange(n-r, n+1)))
    coef = np.exp(logcoef)
    coef[0:r] = 0
    return coef

def lmom_coef(n, k):
    """
    Returns the coefficients for estimating the k^th estimated L-moment of
    a sample of size n.  The L-moment estimate is a linear combination of
    the order statistics with these coefficients.
    """
    b = [pwm_coef(n, r) for r in range(k)]
    if k == 1:
        return b[0]
    elif k == 2:
        return 2*b[1] - b[0]
    elif k == 3:
        return 6*b[2] - 6*b[1] + b[0]
    elif k == 4:
        return 20*b[3] - 30*b[2] + 12*b[1] - b[0]
    else:
        raise("Unimplemented L-moment order")

def lmom(x, k, standardize=False):
    """
    Returns the estimated L-moment of order k.
    """
    xs = np.sort(x)
    lm = np.dot(lmom_coef(len(x), k), xs)
    if standardize:
        if k <= 2:
            warnings.warn("Standardizing L-moment of order <3.")
        lm /= np.dot(lmom_coef(len(x), 2), xs)
    return lm

def l1(x):
    """
    Returns the estimated L-moment of order 1.
    """
    return lmom(x, 1)

def l2(x):
    """
    Returns the estimated L-moment of order 2.
    """
    return lmom(x, 2)

def l3(x):
    """
    Returns the estimated L-moment of order 3.
    """
    return lmom(x, 3)

def l4(x):
    """
    Returns the estimated L-moment of order 4.
    """
    return lmom(x, 4)

def lmom_del(x, k):
    """
    Returns the k^th estimated L-moment of x, deleting each observation in turn.
    """
    ii = np.argsort(x)
    ir = np.argsort(ii)
    z = x[ii] # order statistics
    n = len(x)
    c = lmom_coef(n-1, k)
    f = np.zeros(n)
    f[1:] = np.cumsum(c * z[0:-1])
    b = np.zeros(n)
    b[0:-1] = np.cumsum(c[::-1] * z[1:][::-1])[::-1]
    return (f + b)[ir]

def lmom_jack(x, k):
    """
    Returns the jack-knife pseudo-observations for the k^th L-moments of x.
    """
    n = len(x)
    lm = lmom(x, k)
    lmd = lmom_del(x, k)
    return n*lm - (n - 1)*lmd

def lcomom(x, y, k, standardize=False):
    """
    Returns the estimated k^th L-comoment of x and y.
    """
    ii = np.argsort(y)
    z = x[ii]
    c = lmom_coef(len(x), k)
    cm = np.dot(z, c)
    if standardize:
        cm /= lmom(x, k)
    return cm
