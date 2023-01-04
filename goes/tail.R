library(dplyr)
library(ggplot2)

source("read.R")

df = get_goes(2017)

# Estimate the tail index.  If P(X > t) ~ 1/x^a, then
# the Hill estimator estimates a.
hill = function(x, p=0.001) {
    x = sort(x)
    n = length(x)
    k = ceiling(p*n)
    m = n - k
    lx = log(x[m:length(x)])
    lx = lx - lx[1]
    alpha = 1 / mean(lx)
    return(alpha)
}

# Make a Pareto plot using the upper p proportion
# of the data in x.
pareto_plot = function(x, title, p=0.01) {
    n = length(x[[1]])
    k = ceiling((1 - p) * n)
    q = seq(1, n) / (n + 1)
    q = q[k:length(q)]

    dp = data.frame()
    g = 1
    for (y in x) {
        y = sort(y)
        y = y[k:length(y)]
        dr = data.frame(x=log(y), y=log(1 - q))
        dr$g = g
        g = g + 1
        dp = rbind(dp, dr)
    }
    dp$g = as.factor(dp$g)

    plt = ggplot(aes(x=x, y=y, color=g, group=g), data=dp) + geom_line()
    plt = plt + ggtitle(title)
    plt = plt + labs(x="Observed log quantile", y="Log complementary probability")
    print(plt)
}

# A more ad-hoc measure of tail thickness.
tailratio = function(x, p=0.99, q=0.75, r=0.5, ref=qnorm) {
    qp = quantile(x, p)
    qq = quantile(x, q)
    qr = quantile(x, r)
    numer = (qp - qr) / (qq - qr)
    denom = (ref(p) - ref(r)) / (ref(q) - ref(r))
    return(numer / denom)
}

pdf("tail_r.pdf")

# Check the Hill estimator using Pareto data
cat("Pareto data:\n")
n = 1e6
for (b in c(1, 2, 3, 4)) {
    f1 = runif(n)^(-1/b)
    dd = list()
    for (k in 1:10) {
        dd[[k]] = runif(n)^(-1/b)
    }
    pareto_plot(dd, sprintf("Pareto data with b=%.1f", b))
    for (p in c(1e-5, 1e-4, 1e-3, 1e-2)) {
        alpha = hill(f1, p)
        cat(sprintf("%5d    %8.6f %8d %12.2f\n", b, p, ceiling(p*length(f1)), alpha))
    }
}

# Check the Hill estimator using non-Pareto data with a Pareto tail
cat("\nNon-Pareto data with Pareto tail:\n")
n = 1e6
for (b in c(1, 2, 3, 4)) {
    x = list()
    for (k in 1:10) {
        f1 = runif(n)^(-1/b)
        ii = sample(1:n, ceiling(n/2))
        f1[ii] = -log(runif(length(ii)))
        x[[k]] = f1
    }
    pareto_plot(x, sprintf("Pareto/exponential mixture with b=%.1f", b))
    f1 = x[[1]]
    for (p in c(1e-5, 1e-4, 1e-3, 1e-2)) {
        alpha = hill(f1, p)
        cat(sprintf("%5d    %8.6f %8d %12.2f\n", b, p, ceiling(p*length(f1)), alpha))
    }
}

# What does the Hill estimator do when the tail is not heavy?
n = 1e6
x = list()
for (k in 1:10) {
    x[[k]] = rnorm(n)
}
pareto_plot(x, "Gaussian data")
f1 = x[[1]]
cat("\nGaussian (light-tailed) data:\n")
for (p in c(1e-5, 1e-4, 1e-3, 1e-2)) {
    alpha = hill(f1, p)
    cat(sprintf("%8.6f %8d %12.2f\n", p, ceiling(p*length(f1)), alpha))
}

# Make Pareto plots of the GOES-flux data and first differences.
pareto_plot(list(as.vector(df$Flux1)), "GOES Flux-1 data")
pareto_plot(list(as.vector(df$Flux2)), "GOES Flux-2 data")
pareto_plot(list(as.vector(diff(df$Flux1))), "GOES Flux-1 data (differenced)")
pareto_plot(list(as.vector(diff(df$Flux2))), "GOES Flux-2 data (differenced)")

# Estimate tail parameters for the GOES-flux data and first differences.
for (d in c(F, T)) {
    f1 = as.vector(df$Flux1)
    if (d) {
        cat("\nX-ray flux data (differenced):\n")
        f1 = diff(f1)
    } else {
        cat("\nX-ray flux data:\n")
    }
    for (p in c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) {
        alpha = hill(f1, p)
        cat(sprintf("%8.6f %8d %12.2f\n", p, ceiling(p*length(f1)), alpha))
    }
}

# Make plots of the tail ratios
get_tailratios = function(x) {
    zr = NULL
    for (r in c(0.5, 0.75, 0.9)) {
        for (f in c(0.5, 0.9)) {
            q = r + (1 - r)*f
            for (g in c(0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999)) {
                p = q + (1 - q)*g
                tw = tailratio(x, p=p, q=q, r=r, ref=qexp)
                zr = rbind(zr, c(p, q, r, tw))
            }
        }
    }
    zr = data.frame(zr)
    colnames(zr) = c("p", "q", "r", "tw")
    return(zr)
}

plot_tailratios = function(zr, title) {
    zr$lp = -log(1-zr$p)
    plt = ggplot(aes(x=lp, y=tw, color=interaction(q, r), group=interaction(q, r)), data=zr) + geom_line()
    plt = plt + labs(x="-log(1-p)", y="Tail ratio")
    print(plt)
}

fl1 = df$Flux1
zr = get_tailratios(fl1)
plot_tailratios(zr, "Flux-1 tail ratios")

fl2 = df$Flux2
zr = get_tailratios(fl2)
plot_tailratios(zr, "Flux-2 tail ratios")

fl1d = diff(fl1)
zr = get_tailratios(fl1d)
plot_tailratios(zr, "Differenced flux-1 tail ratios")

fl2d = diff(fl2)
zr = get_tailratios(fl2d)
plot_tailratios(zr, "Differenced flux-2 tail ratios")

dev.off()
