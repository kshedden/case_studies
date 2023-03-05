# See phreg.py for background information.

library(survival)

# Generate n observations from a censored Weibull population
# with hazard proportional to exp(lhr' * [x, z]); x and z
# are standardized covariates such that cor(x, z) = r.
#
# The hazard function depends on shape; if shape = 0 the
# hazard function is constant, if shape > 0 the hazard
# function is increasing in time and if shape < 0 the hazard
# function is decreasing in time.
gendat = function(n, r, shape, lhr) {
    x = rnorm(n)
    z = rnorm(n)
    z = r*x + sqrt(1-r^2)*z

    lp = (lhr[1]*x + lhr[2]*z) / shape
    mn = exp(-lp)

    # Event time
    evttime = mn * rweibull(n, shape=shape)

    # Observation time
    obstime = mean(mn) * rweibull(n, 2)

    status = 1*(evttime <= obstime)
    time = status*evttime + (1 - status)*obstime

    da = data.frame(time=time, status=status, obstime=obstime,
                    evttime=evttime, x=x, z=z)
    return(da)
}

# Generate a sample of observations from a collection of Weibull
# distributions with different shape parameters, provided in the array
# shape.  For other parameters see the docstring for the gendat
# function.
gendat_stratified = function(n, r, shape, lhr) {
    da = NULL
    for (j in 1:length(shape)) {
        dx = gendat(n, r, shape[j], lhr)
        dx["group"] = j
        da = rbind(da, dx)
    }
    return(da)
}

# Use simulation to estimate the power for a given sample size (n)
# correlation between x and z (r), shape parameters, and log hazard
# ratios.
runsim = function(n, r, shape, lhr) {
    zs = array(0, c(nrep, 2))
    for (i in 1:nrep) {
        da = gendat_stratified(n, r, shape, lhr)
        m = coxph(Surv(time, status) ~ x + z + strata(group), data=da)
        zs[i,] = coef(m) / sqrt(diag(vcov(m)))
    }
    mn = apply(zs, 2, mean)
    sd = apply(zs, 2, sd)
    pw = pnorm(-2, mean=mn[1], sd=sd[1]) + 1 - pnorm(2, mean=mn[1], sd=sd[1])
    return(pw)
}

# Estimate power for one or more values of n, r, shapes, lhr.
runsims = function(n, r, shapes, lhr) {
    rslt = data.frame()
    for (n1 in n) {
        for (r1 in r) {
            for (shape1 in shapes) {
                for (lhr1 in lhr) {
                    pw = runsim(n1, r1, shape1, lhr1)
                    row = data.frame(n=n1, r=r1, lhr1=lhr1[1], lhr2=lhr1[2], power=pw)
                    rslt = rbind(rslt, row)
                }
            }
        }
    }
    return(rslt)
}

# nrep ~ 100 gives stable results but is slow
nrep = 50

shapes = c(0.5, 1, 2, 4)
pw = runsims(c(150), c(0.1, 0.2, 0.3, 0.4), list(c(shapes)), list(c(0.15, 0.4)))
