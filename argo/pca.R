source("read.R")

library(ggplot2)
library(ggrastr)
library(dr)

m = dim(temp)[1] # Number of pressure points
n = dim(temp)[2] # Number of profiles

# The mean profile
tempmean = apply(temp, 1, mean)
psalmean = apply(psal, 1, mean)

# Center the profiles
tempc = temp - outer(tempmean, array(1, n))
psalc = psal - outer(psalmean, array(1, n))

# Get the principal components
cc = cov(t(tempc))
ee = eigen(cc, symmetric=T)
eigval = ee$values
eigvec = ee$vectors

# Flip each PC loading vector if it is mostly negative.
for (j in 1:dim(eigvec)[2]) {
    if (sum(eigvec[,j] < 0) > sum(eigvec[,j] > 0)) {
        eigvec[,j] = -eigvec[,j]
    }
}

# Scores for the dominant PC's
scores = t(tempc) %*% eigvec[,1:10]

pdf("pca_r.pdf")

# Plot the mean profile
da = data.frame(pressure=pressure, tempmean=tempmean)
plt = ggplot(aes(x=pressure, y=tempmean), data=da) + geom_line()
print(plt)

# Make several plots depicting the loading patterns.
for (j in 1:5) {
    # Plot the loadings
    da = data.frame(pressure=pressure, loading=eigvec[,j])
    plt = ggplot(aes(x=pressure, y=loading), data=da) + geom_line()
    plt = plt + labs(x="Pressure", y="Temperature loading")
    plt = plt + ggtitle(sprintf("PC %d", j))
    print(plt)

    # Plot the mean profile +/- multiples of the loadings
    da = data.frame()
    s = sd(scores[,j])
    for (f in c(-2, -1, 0, 1, 2)) {
        dx = data.frame(pressure=pressure, profile=tempmean+f*s*eigvec[,j], f=f)
        da = rbind(da, dx)
    }
    da$f = as.factor(da$f)
    plt = ggplot(aes(x=pressure, y=profile, color=f, group=f), data=da) + geom_line()
    plt = plt + labs(x="Pressure", y="Temperature")
    plt = plt + ggtitle(sprintf("PC %d", j))
    print(plt)
}

# Plot the j^th PC score against the k^th feature.
plot_pcscores = function(j, k) {

    xx = seq(min(Y[,k]), max(Y[,k]), length=100)
    m = lowess(Y[,k], scores[,j])
    mf = approxfun(m$x, m$y)
    f = abs(m$y - scores[,j])
    r = lowess(Y[,k], f)
    rf = approxfun(r$x, r$y)
    da = data.frame(x=xx, y=mf(xx), r=rf(xx))

    f = 2
    da$y1 = da$y - f*da$r
    da$y2 = da$y + f*da$r
    plt = ggplot(aes(x=x, y=y), data=da)
    plt = plt + labs(x=fn[k], y=sprintf("PC %d score", j))
    plt = plt + geom_ribbon(aes(x=x, ymin=y1, ymax=y2), fill="grey70")
    plt = plt + geom_line()
    print(plt)
}

fn = c("Latitude", "Longitude", "Day")
Y = cbind(lat, lon, day)

for (j in 1:3) {
    for (k in 1:3) {
        plot_pcscores(j, k)
    }
}

# Flip signs in two sets of CCA loadings to make them easier to understand in a
# plot.
flip = function(xcoef, ycoef) {
    for (j in 1:dim(xcoef)[2]) {
        if (mean(xcoef[,j] > 0) + mean(ycoef[,j] > 0) < 1) {
            xcoef[,j] = -xcoef[,j]
            ycoef[,j] = -ycoef[,j]
        }
    }
    return(list(xcoef=xcoef, ycoef=ycoef))
}

# Use a combined PCA/CCA approach to relate temperature and salinity.
# Produce plots for a sequence of dimensions for the PCA dimension
# reduction.
plot_cca = function() {
    X = t(tempc)
    Y = t(psalc)
    svx = svd(X)
    svy = svd(Y)

    for (q in c(1, 2, 5, 10)) {
        cc = cancor(svx$u[,1:q], svy$u[,1:q])
        rr = cc$cancorr

        # Passing nrow and ncol is needed to avoid
        # different behavior when q=1.
        ddx = diag(svx$d[1:q], q, q)
        ddy = diag(svy$d[1:q], q, q)

        xcoef = svx$v[,1:q] %*% solve(ddx, cc$xcoef)
        ycoef = svy$v[,1:q] %*% solve(ddy, cc$ycoef)
        ll = flip(xcoef, ycoef)
        xcoef = ll$xcoef
        ycoef = ll$ycoef
        da = data.frame(temp_weights=xcoef[,1], psal_weights=ycoef[,1], pressure=pressure)
        plt = ggplot(aes(x=pressure, y=temp_weights), data=da) + geom_line()
        plt = plt + labs(x="Pressure", y="Temperature coefficient")
        plt = plt + ggtitle(sprintf("CCA with q=%d components, r=%.2f", q, rr[1]))
        if (min(ycoef[,1]) > 0) {
            plt = plt + ylim(0, NA)
        }
        print(plt)
        plt = ggplot(aes(x=pressure, y=psal_weights), data=da) + geom_line()
        plt = plt + labs(x="Pressure", y="Salinity coefficient")
        plt = plt + ggtitle(sprintf("CCA with q=%d components, r=%.2f", q, rr[1]))
        if (min(ycoef[,1]) > 0) {
            plt = plt + ylim(0, NA)
        }
        print(plt)
    }
}

plot_cca()

# SIR is a dimension reduction regression (DR) method.  Here we will use
# it doidentify factors within the temperature data that predict latitude.
plot_sir = function() {

    X = t(tempc)

    # Due to high-dimensionality, we project the temperature data to
    # a limited number of PC's before using SIR to estimate the coefficients.
    q = 5
    svx = svd(X)
    dd = dr.compute(svx$u[,1:q], lat, array(1, length(lat)))
    ddx = diag(svx$d[1:q], q, q)
    b = svx$v[,1:q] %*% solve(ddx, dd$evectors)

    # Plot the loadings
    da = data.frame(pressure=c(pressure, pressure, pressure))
    da$dir = c(b[,1], b[,2], b[,3])
    oo = array(1, length(pressure))
    da$group = c(1*oo, 2*oo, 3*oo)
    da$group = as.factor(da$group)
    plt = ggplot(aes(x=pressure, y=dir, by=group, color=group), data=da) + geom_line()
    plt = plt + labs(x="Pressure", y="Temperature loading")
    plt = plt + ggtitle(sprintf("SIR with q=%d PC components", q))
    print(plt)

    # Plot the scores against latitude, longitude, and day.
    scores = X %*% b[,1:3]
    for (v in c("lat", "lon", "day")) {
        for (j in 1:3) {
            da = data.frame(lat=lat, lon=lon, day=day, s=scores[,j])
            plt = ggplot(aes(x=!!sym(v), y=s), data=da) + geom_point()
            plt = plt + labs(x=v, y=sprintf("SIR component %d", j))
            rasterize(plt, layers="Point")
            print(plt)
        }
    }
}

da = plot_sir()

dev.off()
