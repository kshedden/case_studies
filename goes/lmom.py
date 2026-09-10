import numpy as np

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

def l1(x):
    """
    Returns the estimated L-moment of order 1.
    """
    return np.dot(lmom_coef(len(x), 1), np.sort(x))

def l2(x):
    """
    Returns the estimated L-moment of order 2.
    """
    return np.dot(lmom_coef(len(x), 2), np.sort(x))

def l3(x):
    """
    Returns the estimated L-moment of order 3.
    """
    return np.dot(lmom_coef(len(x), 3), np.sort(x))

def l4(x):
    """
    Returns the estimated L-moment of order 4.
    """
    return np.dot(lmom_coef(len(x), 4), np.sort(x))
