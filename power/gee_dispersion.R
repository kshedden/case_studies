# See gee_dispersion.py for background

library(geepack)
library(dplyr)

# Generate a Gaussian n-vector with mean zero and exchangeable correlation icc.
gen_exch = function(n, icc) {
    return(sqrt(icc)*rnorm(1) + sqrt(1 - icc)*rnorm(n))
}
# Generate clustered data with a Poisson or Gamma response variable (y),
# and two covariates (x, y).
#
# ngrp: number of groups
# m: average size of each group
# r: correlation between covariates x and z
# y_icc: exchangeable correlation of the outcome
# x_icc: exchangeable correlation of covariate x
# z_icc: exchangeable correlation of covariate z
# b: coefficients for x and z
# scale: the glm scale parameter
# fam: the family used to generate the data
gendat = function(ngrp, m, r, y_icc, x_icc, z_icc, b, scale, fam) {

    x = list()
    z = list()
    g = list()
    e = list()

    for (i in 1:ngrp) {
        mm = 1 + rpois(1, m) # group size
        x1 = gen_exch(mm, x_icc)
        x[[i]] = x1
        u = gen_exch(mm, z_icc)
        z[[i]] = r*x[[i]] + sqrt(1 - r^2)*u
        e[[i]] = gen_exch(mm, y_icc)
        g[[i]] = i*array(1, mm)
    }

    x = unlist(x)
    z = unlist(z)
    e = unlist(e)
    g = unlist(g)

    # The marginal mean
    mn = exp(b[1]*x + b[2]*z)

    u = pnorm(e)

    if (fam[[1]] == "gamma") {
        # Mean is a*b, variance is a*b^2
        b = scale*mn
        a = 1 / scale
        y = qgamma(u, a, scale=b)
    } else {
        y = qpois(u, mn)
    }

    da = data.frame(x=x, z=z, y=y, g=g)
    return(da)
}

# Estimate the power in one settting using simulation.
run_power = function(ngrp, m, r, y_icc, x_icc, z_icc, b, scale, mfam, gfam) {

    zs = array(0, c(nrep, 3))
    for (i in 1:nrep) {
        da = gendat(ngrp, m, r, y_icc, x_icc, z_icc, b, scale, gfam)
        if (mfam[[1]] == "Gamma") {
            da = da %>% mutate(y=y+0.1)
        } else {
            da = da %>% mutate(y=round(y))
        }
        r0 = geeglm(y ~ x + z, id=da$g, family=mfam, data=da)
        zs[i,] = coef(r0) / sqrt(diag(vcov(r0)))
    }
    mn = mean(zs[, 2])
    sd = sd(zs[, 2])
    pw = pnorm(-2, mean=mn, sd=sd) + 1 - pnorm(2, mean=mn, sd=sd)
    return(pw)
}

nrep = 100

# Generate gamma and Poisson data from the same mean structure, analyze it
# with either a gamma or Poisson GEE.
b = c(0.15, 0.4)
for (s in c(0, 1)) {
    b1 = s*b

    cat(c("Null hypothesis is true:\n", "\nAlternative hypothesis is true:\n")[1 + s])

    mfam = Gamma(link="log")
    pw = run_power(200, 5, 0.3, 0.2, 0.2, 0.2, b1, 4, mfam, "gamma")
    cat(sprintf("%6.3f gamma data gamma model\n", pw))
    pw = run_power(200, 5, 0.3, 0.2, 0.2, 0.2, b1, 4, mfam, "poisson")
    cat(sprintf("%6.3f Poisson data gamma model\n", pw))

    mfam = poisson()
    pw = run_power(200, 5, 0.3, 0.2, 0.2, 0.2, b1, 4, mfam, "gamma")
    cat(sprintf("%6.3f gamma data Poisson model\n", pw))
    pw = run_power(200, 5, 0.3, 0.2, 0.2, 0.2, b1, 4, mfam, "poisson")
    cat(sprintf("%6.3f Poisson data Poisson model\n", pw))
}
