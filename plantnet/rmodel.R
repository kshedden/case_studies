library(dplyr)
library(splines)
library(lme4)
library(ggplot2)

source("read.R")

pdf("rmodel_r.pdf")

# Each of these variables can be the outcome, and the other
# two are predictors.
va = c("decimalLatitude", "decimalLongitude", "elevation")

# Names for the variables in 'va' used for labeling plots.
titles = list(decimalLatitude="Latitude", decimalLongitude="Longitude", elevation="Elevation")

df = df[complete.cases(df),]

# Fit a multilevel regression with variable 'v' as the outcome using
# the data in dataframe 'df'.  The variables in 'va' that are not in
# 'v' are used as coviates (modeled with splines).
fitmodel = function(v, df) {

    # Control for the variables except for the outcome
    te = va[va != v]
    tx = sprintf("bs(%s, 5)", te)

    s = paste(tx, collapse=" + ")
    fml = sprintf("%s ~ year + %s + (1 + decade | scientificName)", v, s)
    fml = as.formula(fml)

    # Fit the mixed model and print the model summary.
    m0 = lmer(fml, data=df)
    cat(sprintf("n=%d observations", dim(model.matrix(m0))[1]))
    rr = ranef(m0)$scientificName
    cat(sprintf("n=%d species", length(rr)))
    print(summary(m0))
    cat("\n\n")

    return(list(mm=m0, rr=rr))
}

# Generate a plot of predicted trends for each species.
make_plots = function(mm, rr, va) {

    # The fitted value at the central location, used as an intercept.
    xx = model.matrix(mm)
    xm = colMeans(xx)
    xm[2] = 0
    icept = sum(xm * fixef(mm))

    # The year slope (in year units.)
    ys = fixef(mm)[2]

    # Plot these years
    yr = seq(2010, 2020, length.out=10)

    # The plotted years as decades
    yx = (yr - meanyear) / 10

    nspecies = dim(rr)[1]
    dp = array(0, c(10*nspecies, 3))
    ii = 0
    for (i in 1:nspecies) {
        dp[(ii+1):(ii+10), 1] = i
        dp[(ii+1):(ii+10), 2] = yr
        dp[(ii+1):(ii+10), 3] = icept + rr[i, 1] + ys*yr + rr[i, 2]*yx
        ii = ii + 10
    }
    dp = data.frame(dp)
    colnames(dp) = c("g", "x", "y")

    plt = ggplot(aes(x=x, y=y, group=g), data=dp) + geom_line(color="grey")
    plt = plt + labs(x="Year", y=titles[[v]])
    print(plt)
}

for (v in va) {
    fm = fitmodel(v, df)
    mm = fm$mm
    rr = fm$rr
    make_plots(mm, rr, v)
}

dev.off()
