library(ggplot2)
source("read.R")

# See this reference for more information about the
# support point algorithm:
# https://arxiv.org/pdf/1609.01811.pdf

# Equation 22 in Mak et al.
update_support = function(X, Y) {
    p = dim(Y)[1]
    N = dim(Y)[2]
    n = dim(X)[2]
    XX = array(0, c(p, n))

    for (i in 1:n) {
        Dx = X[,i] - X
        DxN = sqrt(colSums(Dx^2))
        DxN[i] = Inf
        Dy = X[,i] - Y
        DyN = sqrt(colSums(Dy^2))
        q = sum(1 / DyN)
        XX[,i] = (Dx %*% (1/DxN)) * (N / n)
        XX[,i] = XX[,i] + Y %*% (1/DyN)
        XX[,i] = XX[,i] / q
    }

    return(XX)
}

# Calculate N support points for the data in Y.  The points
# are stored in the columns of Y.
support = function(Y, N, maxiter=1000) {

    p = dim(Y)[1]
    n = dim(Y)[2]
    X = array(rnorm(N*p), c(p, N))

    for (i in 1:maxiter) {
        X1 = update_support(X, Y)
        ee = norm(X1 - X)
        X = X1
        if (ee < 1e-8) {
            break
        }
    }

    return(X)
}

pdf("support_R.pdf")

# Plot support points for temperature and salinity separately.
for (j in 1:2) {

    if (j == 1) {
        x = temp
        ylab = "Temperature"
    } else {
        x = psal
        ylab = "Salinity"
    }

    # Make plots with different numbers of support points.
    for (npt in c(1, 2, 5, 10, 20)) {
        cat(sprintf("npt=%d\n", npt))
        X = support(x, npt, maxiter=200)

        # Put the support points in long form for plotting
        da = data.frame(x=array(X), g=kronecker(1:npt, array(1, 100)),
                        pressure=kronecker(array(1, npt), pressure))

        plt = ggplot(aes(x=pressure, y=x, group=g), data=da) + geom_line()
        plt = plt + labs(x="Pressure", y=ylab) + ggtitle(sprintf("%d support points", npt))
        print(plt)
    }
}

dev.off()
