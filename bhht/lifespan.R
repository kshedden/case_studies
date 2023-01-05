library(dplyr)
library(readr)
library(ggplot2)
library(locpol)

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"

da = read_csv(file.path(pa, "cross-verified-database.csv.gz"))

da$lifespan = da$death - da$birth

da = da[, c("birth", "lifespan", "gender", "level1_main_occ", "un_region")]

da = rename(da, "occ"="level1_main_occ")
da = rename(da, "reg"="un_region")
da = rename(da, "sex"="gender")

# There are too few people with "Other" gender to estimate
# the conditional mean lifespans.
da = filter(da, sex %in% c("Female", "Male"))

# Focus on the last 500 years.
da = filter(da, birth >= 1500)

# People born after 1920 may still be alive, which leads to censoring.
da = filter(da, birth <= 1920)

dx = da[, c("birth", "occ", "sex", "reg", "lifespan")]
dx = dx[complete.cases(dx),]
dx = dx[order(dx$birth),]

# Plot the conditional dispersion for local polynomial regression compared
# to a linear model.
m0 = lm(lifespan ~ birth, dx)
ares0 = abs(resid(m0))
bw = thumbBw(dx$birth, dx$lifespan, 1, EpaK)
m1 = locLinSmootherC(dx$birth, dx$lifespan, xeval=seq(1500, 1920, length.out=500), bw=bw, kernel=EpaK)
f = approxfun(m1$x, m1$beta0)
ares1 = abs(f(dx$birth) - dx$lifespan)
ma0 = lowess(dx$birth, ares0, f=0.2)
ma1 = lowess(dx$birth, ares1, f=0.2)
dp = data.frame(birth=c(ma0$x, ma1$x), resid=c(ma0$y, ma1$y))
dp$method = "linear"
n = length(ma0$x)
dp$method[(n+1):(2*n)] = "locpoly"
plt0 = ggplot(data=dp, aes(x=birth, y=resid, group=method)) + geom_line(aes(color=method))

# Estimate the conditional mean of the variable 'va' given the value of
# 'birth' on a grid of 100 years spanning from 1500 to 1920.
cmest = function(va, dx) {
    bb = seq(1500, 1920, length.out=100)
    rr = dx %>% group_by(!!sym(va)) %>% group_modify(~ {
        m = lowess(.x$birth, .x$lifespan)
        f = approxfun(m)
        data.frame(birth=bb, lifespan=f(bb))
    })
    return(rr)
}

# Use bootstrapping to estimate the standard errors for the
# quantities calculated by the function 'cmest'.
cmest_boot = function(va, dx, nboot=10) {
    br = NULL
    for (i in 1:nboot) {
        dxb = dx %>% slice_sample(n=dim(dx)[1], replace=T)
        rr = cmest(va, dxb)
        if (i == 1) {
            br = array(0, c(dim(rr)[1], nboot))
        }
        br[,i] = rr$lifespan
    }
    return(apply(br, 1, sd))
}

# Generate a plot of the estimated conditional mean lifespan given
# year of birth, for data stratified according to the variable in 'va'.
# If se=T plot a pointwise confidence band around the conditional mean
# estimates.
cmplot = function(va, dx, se=F) {
    rr = cmest(va, dx)
    if (se) {
        s = cmest_boot(va, dx)
        rr$lcb = rr$lifespan - 2*s
        rr$ucb = rr$lifespan + 2*s
        plt = ggplot(data=rr, aes(x=birth, y=lifespan, group=!!sym(va), color=!!sym(va), fill=!!sym(va)))
        plt = plt + geom_line()
        plt = plt + geom_ribbon(aes(ymin=lcb, ymax=ucb), alpha=0.3, linetype=0)
    } else {
        plt = ggplot() + geom_line(data=rr, aes(x=birth, y=lifespan, color=!!sym(va)))
    }
    return(plt)
}

# Estimate the difference of conditional means between sexes (female - male)
# for lifespan given year of birth, using the data in 'dx'.
cmdiff = function(dx) {
    rr = cmest("sex", dx)
    s = cmest_boot("sex", dx)
    n = dim(rr)[1]

    # This function only works for dataframes with m rows of
    # female data and m rows of male data.
    stopifnot(rr$sex[1] == "Female")
    stopifnot(rr$sex[n/2+1] == "Male")

    # Female - male mean lifespan
    d = rr$lifespan[1:(n/2)] - rr$lifespan[(n/2+1):n]

    # Standard error of d
    s = sqrt(s[1:(n/2)]^2 + s[(n/2+1):n]^2)

    # Confidence band
    lcb = d - 2*s
    ucb = d + 2*s

    r2 = data.frame(birth=rr$birth[1:(n/2)], d=d, lcb=lcb, ucb=ucb)
    return(r2)
}

# Generate a plot displaying the difference of the conditional mean lifespan
# for females relative to males (females - males), using the data in 'dx'.
cmdiff_plot = function(dx) {
    rr = cmdiff(dx)
    plt = ggplot(data=rr, aes(x=birth, y=d))
    plt = plt + geom_line()
    plt = plt + geom_ribbon(aes(ymin=lcb, ymax=ucb), alpha=0.3, linetype=0)
    return(plt)
}

plt1 = cmplot("sex", dx, T)
plt2 = cmplot("reg", dx)
plt3 = cmplot("occ", dx)
plt4 = cmdiff_plot(dx)

pdf("lifespan_r_lowess.pdf")
print(plt0)
print(plt1)
print(plt2)
print(plt3)
print(plt4)
dev.off()
