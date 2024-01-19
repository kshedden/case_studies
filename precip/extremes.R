library(readr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(ROOPSD)
library(dplyr)

# Estimate the parameters of a generalized Pareto distribution
# using the empirical Bayes method of Zhang and Stephens.
#
# https://www.jstor.org/stable/pdf/40586625.pdf
gp_estimate = function(z) {

    n = length(z)
    xstar = quantile(z, 0.25)
    m = ceiling(20 + sqrt(n))
    xmax = max(z)

    tgrid = 1/xmax + (1 - sqrt(m / (seq(1, m) - 0.5))) / (3 * xstar)

    profile = function(theta) {
        k = -mean(log(1 - theta*z))
        return(n*(log(theta/k) + k - 1))
    }

    ltg = sapply(tgrid, profile)
    ltg = ltg - max(ltg)
    Ltg = exp(ltg)
    Ltg = Ltg / sum(Ltg)
    theta_hat = sum(Ltg * tgrid)
    k_hat = -mean(log(1 - theta_hat*z))
    sigma_hat = k_hat / theta_hat

    return(list(scale=sigma_hat, shape=-k_hat))
}

# Use parameteric bootstrap with a Gaussian copula
# to assess the sampling performance of the empirical
# Bayes estimates of the GPD parameters.
gp_simstudy_copula = function(n, r, eb0, nrep) {
    scales = NULL
    shapes = NULL
    for (k in 1:nrep) {
        zz = rnorm(n)
        zz[2:n] = r*zz[1:n-1] + sqrt(1-r^2)*zz[2:n]
        uu = pnorm(zz)
        zz = qgpd(uu, scale=eb0$scale, shape=eb0$shape)
        eb = gp_estimate(zz)
        scales[length(scales)+1] = eb$scale
        shapes[length(shapes)+1] = eb$shape
    }

    dp = data.frame(scales=scales, shapes=shapes)
    return(dp)
}

# Assess the empirical Bayes estimate of the GPD parameters.
# Use copula to induce serial dependence.
gp_simstudy = function(z, nrep=500) {

    # Exceedances
    z = z[z > thresh] - thresh
    n = length(z)

    # Empirical Bayes estimate of Zhang and Stephens.
    eb0 = gp_estimate(z)

    dpl = list()
    rl = c(0.0, 0.5)

    cat("Simulation study for empirical Bayes estimates of GPD parameters:\n")
    for (r in rl) {
        dp = gp_simstudy_copula(1000, r, eb0, nrep)
        cat(sprintf("r=%.2f\n", r))
        cat(sprintf("Bias(shape)=%.2f\n", mean(dp$shapes) - eb0$shape))
        cat(sprintf("Bias(scales)=%.2f\n", mean(dp$scales) - eb0$scale))
        cat(sprintf("SD(shape)=%.3f\n", sd(dp$shapes)))
        cat(sprintf("SD(scale)=%.3f\n", sd(dp$scales)))
        dpl[[length(dpl)+1]] = dp
    }

    dp = dpl[[1]] # make plots for uncorrelated case
    p = ggplot(aes(x=scales), data=dp) + geom_histogram()
    p = p + xlab("Scale parameter") + ylab("Frequency")
    p = p + ggtitle("Sampling distribution of scale parameter")
    print(p)

    p = ggplot(aes(x=shapes), data=dp) + geom_histogram()
    p = p + xlab("Shape parameter") + ylab("Frequency")
    p = p + ggtitle("Sampling distribution of shape parameter")
    print(p)
}

# Fit a generalized extreme value distribution using maximum likelihood
# estimation.  Probability weighted moments are used to obtain starting
# values:
#
# https://www.stat.cmu.edu/technometrics/80-89/VOL-27-03/v2703251.pdf
fit_gev = function(x) {

    x = sort(x)
    n = length(x)

    # The first three probability weighted moments
    b = array(0, 3)
    jj = seq(1, n) / (n + 1)
    for (r in 1:3) {
        b[r] = sum(jj^(r-1) * x) / n
    }

    # The PWM estimator of Hoskins et al.
    c = (2*b[2] - b[1]) / (3*b[3] - b[1])  - log(2) / log(3)
    shape = 7.8590*c + 2.9554*c^2
    scale = (2*b[2] - b[1]) * shape / (gamma(1 + shape) * (1 - 1/2^shape))
    loc = b[1] + scale*(gamma(1 + shape) - 1) / shape

    # Get the MLE
    fll = function(par) {
        loc = par[1]
        scale = par[2]
        shape = par[3]
        ll = dgev(x, loc=loc, scale=scale, shape=shape, log=TRUE)
        return(-sum(ll))
    }
    x0 = c(loc, scale, -shape)
    rr = optim(x0, fll)
    pa = rr$par

    return(list(loc=pa[1], scale=pa[2], shape=pa[3]))
}

# Returns values x, p such that the slope of p on x
# estimates the shape parameter (tail index) of
# a distribution with power-law tails, or the rate parameter
# of a distribution with exponential tails.
# The upper p0 fraction of the data in z are used
# to produce the returned values in (x, p).  The values in x are the order
# statistics of z in the Exponential case, and the log
# order statistics of z in the Pareto case. The values in
# p are derived from probability points.
tail_shape = function(z, p0=0.1, family="powerlaw") {
    if (!(family %in% c("exponential", "powerlaw"))) {
        error(sprintf("Unknown family %s", family))
    }
    p = 1 - p0
    z = sort(z)
    n = length(z)
    m = as.integer(round(p*n))
    x = z[m:length(z)]
    if (family == "powerlaw") {
        x = log(x)
    }
    p = log(1 - seq(m, n) / (n + 1))
    return(list(x=x, p=p))
}

# Use least squares regression in the tail of a quantile plot
# to estimate the shape parameter.
fit_tail_reg = function(x, p0=0.99, family="powerlaw") {

    ts = tail_shape(x, p0=p0, family=family)

    # Estimate the tail index using a least squares fit to the order
    # statistics.
    alpha_hat = -cov(ts$p, ts$x) / var(ts$x)
    icept = mean(ts$p) + alpha_hat*mean(ts$x)

    return(list(icept=icept, alpha_hat=alpha_hat, x=ts$x, p=ts$p))
}

# Make a probability plot of the tails toassess goodness of fit toa
# power law or exponential tailed distribution.
plot_tails = function(z, p0, thresh, family) {

    n = length(z)

    # Exceedances
    z = z[z >= thresh]
    z = z - thresh

    # The number of selected observations
    m = as.integer(round(p0*n))

    xlabel = ifelse(family == "powerlaw", "log Q(p)", "Q(p)")
    ylabel = "log(1-p)"

    tr = fit_tail_reg(z, p0=p0, family=family)
    title = sprintf("%s model, threshold=%.1f, top %.1f%% (n=%d), alpha=%.3f",
                    str_to_title(family), thresh, 100*p0, m, tr$alpha)

    xx = tr$x
    pp = tr$p
    pf = tr$icept - tr$alpha_hat*xx
    dp = data.frame(x=xx, p=pp, pfit=pf)
    p = ggplot(aes(x=x, y=p), data=dp) + geom_point() + geom_line(aes(x=x, y=pfit))
    p = p + ggtitle(title) + xlab(xlabel) + ylab(ylabel)
    print(p)
}

# Calculate the maximum precipitation value for each complete year,
# and fit a generalized extreme value (GEV) distribution to the
# data.  Then use the fitted model to calculate returns for a sequence
# of time horizons, and create a QQ plot to assess goodness-of-fit.
block_max = function(annmax) {

    # Fit a generalized extreme value distribution to the block maxima.
    gev = fit_gev(annmax$PRCP)

    # m-observation returns
    mr = data.frame(Years=c(10, 100, 500, 1000))
    mr$Return = qgev(1 - 1/mr$Years, loc=gev$loc, scale=gev$scale, shape=gev$shape)
    cat("\nReturns based on GEV:\n")
    print(mr)

    # Make a QQ plot to assess goodness of fit
    z = sort(annmax$PRCP)
    n = dim(annmax)[1]
    pp = seq(1, n) / (n+1)
    qq = qgev(pp, loc=gev$loc, scale=gev$scale, shape=gev$shape)

    dp = data.frame(pp=pp, qq=qq, z=z)
    p = ggplot(aes(x=qq, y=z), data=dp) + geom_line()
    p = p + ylab("Order statistics") + xlab("GEV quantiles")
    p = p + ggtitle("GEV fit to annual maxima")
    print(p)

    return(gev)
}

# Calculate the m-observation returns for the data in z, using either
# an exponential or generalized Pareto model.
mobs_return = function(z, mr, thresh=thresh, family="exponential", gp=NULL) {

    n = length(z)

    # Select only extreme values and translate back to the origin
    ix = z >= thresh
    q = sum(ix) / n # proportion of values exceeding the threshold
    z = z[ix]
    z = z - thresh

    pr = 1 - 1 / (q*mr)

    if (family == "exponential") {
        mn = mean(z)
        m0 = thresh - mn*log(1 - pr)
    } else if (family == "generalizedpareto") {
        m0 = thresh + qgpd(pr, scale=gp$scale, shape=gp$shape)
    } else {
        stop("!!")
    }

    return(m0)
}

# Change this to point to the location of the data, matching the path
# name in get_data.R
target_dir = "/home/kshedden/data/Teaching/precip"

# Choose a specific location to analyze.
fname = "USW00094847.csv" # Detroit
#fname = "USW00012839.csv" # Miami

pdf_out = pdf("extremes_R.pdf")

df = read_csv(file.path(target_dir, sprintf("%s.gz", fname)))

df = df %>% select(DATE, PRCP) %>% drop_na

# Convert precipitation to millimeters
df = df %>% mutate(PRCP = PRCP / 10)

# Use this threshold for calculating exceedances
thresh = 5.0

# Annual maxima
df = df %>% mutate(year = year(DATE))
annmax = df %>% group_by(year) %>% summarize(PRCP = max(PRCP), n=n())
annmax = annmax %>% filter(n > 350)

# Quantile plots of the tail of the distribution of 24 hour rainfall totals
for (family in c("powerlaw", "exponential")) {
    for (p0 in c(0.5, 0.1, 0.05, 0.01)) {
        plot_tails(df$PRCP, p0, thresh, family)
    }
}

# Time series plot of precipitation data
p = ggplot(aes(x=DATE, y=PRCP), data=df) + geom_line()
p = p + xlab("Date") + ylab("Precipitation (mm)")
print(p)

# Histogram of precipitation data
p = ggplot(aes(y=PRCP), data=df) + geom_histogram()
p = p + ylab("Frequency") + xlab("Precipitation (mm)")
print(p)

# Plot the empirical cumulative distribution function (eCDF) of precipitation data
n = dim(df)[1]
dp = data.frame(y=(seq(n)-1/2)/n, x=sort(df$PRCP))
p = ggplot(aes(x=x, y=y), data=dp) + geom_step()
p = p + ylab("Cumulative probability") + xlab("Precipitation (mm)")
print(p)

# Plot the complementary cumulative distribution function (eCCDF) of precipitation data
n = dim(df)[1]
dp = data.frame(y=1-(seq(n)-1/2)/n, x=sort(df$PRCP))
p = ggplot(aes(x=x, y=y), data=dp) + geom_step()
p = p + ylab("Complementary cumulative probability") + xlab("Precipitation (mm)")
print(p)

# Construct exceedances
z = df$PRCP - thresh
z = z[z > 0]

eb = gp_estimate(z)
gp_simstudy(df$PRCP)
gev = block_max(annmax)

# Calculate m-observation returns based on various models fit
# to the 24 hour rainfall totals.
yr = c(1, 10, 100, 500, 1000)
cfg = list(list(f="exponential", gp=NULL), list(f="generalizedpareto", gp=eb))
for (rr in cfg) {
    f = rr$f
    cat(sprintf("\nM-observation returns for %s model:\n", f))
    print(f)
    mr = mobs_return(df$PRCP, 365 * yr, family=rr$f, thresh=thresh, gp=rr$gp)
    rr = data.frame(Years=yr, MR=mr)
    print(rr)
}

dev.off()
