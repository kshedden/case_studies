import numpy as np
from scipy.special import comb

def l1(x):
    return x.mean()

def l2(x):
    x = np.sort(x)
    n = len(x)
    ix = np.arange(n)
    m = np.dot(x, comb(ix, 1) - comb(n-ix-1, 1))
    m /= (2 * comb(n, 2))
    return m

def l3(x):
    x = np.sort(x)
    n = len(x)
    ix = np.arange(n)
    m = np.dot(x, comb(ix, 2) - 2*comb(ix, 1)*comb(n-ix-1, 1) + comb(n-ix-1, 2))
    m /= (3 * comb(n, 3))
    return m

def l4(x):
    x = np.sort(x)
    n = len(x)
    ix = np.arange(n)
    m = np.dot(x, comb(ix, 3) - 3*comb(ix, 2)*comb(n-ix-1, 1) + 3*comb(ix, 1)*comb(n-ix-1, 2) - comb(n-ix-1, 3))
    m /= (4 * comb(n, 4))
    return m
