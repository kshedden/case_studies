import numpy as np
from read import *

df = get_goes(2017)

nn = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]

# Estimate the Hurst parameter.
def hurst(df, nn, d):

    r = np.zeros((len(nn), 2))
    for j, m in enumerate(nn):

        # Generate a matrix of non-overlapping blocks of
        # size m.
        _, flx = make_blocks(df, m, d)

        # Calculate the sample mean of each block.
        bm = flx.mean(1)

        # Take the sample variance of the block means.
        r[j, :] = [m, bm.var()]

    # Estimate the Hurst exponent from the variances of
    # the sample means.
    rl = np.log(r)
    b = np.cov(rl[:, 0], rl[:, 1])[0, 1] / np.var(rl[:, 0])
    a = rl[:, 1].mean() - b*rl[:, 0].mean()

    return 1 + b/2

# As a check, estimate the Hurst parameter for 
# IID normal data (the true value is 1/2). 
dx = df.iloc[0:100000, :].copy()
dx["Flux1"] = np.random.normal(size=100000)
h0 = hurst(dx, nn, 0)
print("Estimated Hurst parameter for IID standard normal data:")
print(h0)

# As another check, simulate correlated data
# with short-range dependence (the true value
# is 1/2).
fx = np.random.normal(size=dx.shape[0])
r = 0.5
for i in range(1, len(fx)):
    fx[i] = r*fx[i-1] + np.sqrt(1 - r**2)*fx[i]
dx["Flux1"] = fx
h0 = hurst(dx, nn, 0)
print("\nEstimated Hurst parameter for short-range dependent normal data:")
print(h0)

# Estimate the Hurst Parameter for the GOES data.
print("\nEstimated Hurst parameter for blocks of Flux-1 data:")
for dx in np.array_split(df, 20):
    h0 = hurst(dx, nn, 0)
    h1 = hurst(dx, nn, 1)
    print([h0, h1])
