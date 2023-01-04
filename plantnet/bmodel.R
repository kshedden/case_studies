# Assess the evidence that the mean latitude, longitude, or elevation
# of occurrences of a species are changing in time.  This is done by
# fitting regression models with least squares, in which one of latitude,
# longitude, or elevation is the dependent variable, and the other two are
# independent variables.  Two models are fit -- one including year
# and one not including year, and we assess the statistical significance
# of the year effect by comparing these two models.

library(readr)
library(dplyr)
library(lubridate)
library(splines)
library(ggplot2)

source("read.R")

# Generate basis functions for latitude, longitude, elevation,
# and day within year (seasonality).  Seasonality and longitude
# are circular variables.
setbasis = function(df) {

    # Basis functions for season.
    per = 365
    for (k in 1:4) {
        s = sprintf("sin_day_%d", k)
        df = df %>% mutate(!!s := sin(2*pi*dayOfYear/per))
        s = sprintf("cos_day_%d", k)
        df = df %>% mutate(!!s := cos(2*pi*dayOfYear/per))
        per = per / 2
    }

    # Basis functions for longitude.
    per = 180
    for (k in 1:4) {
        s = sprintf("sin_lon_%d", k)
        df = df %>% mutate(!!s := sin(2*pi*decimalLongitude/per))
        s = sprintf("cos_lon_%d", k)
        df = df %>% mutate(!!s := cos(2*pi*decimalLongitude/per))
        per = per / 2
    }

    # Basis functions for latitude.
    x = (df$decimalLatitude - 45) / 100
    for (k in 1:3) {
        s = sprintf("lat%d", k)
        df = df %>% mutate(!!s := x^k)
    }

    # Basis functions for elevation.
    x = (df$elevation - 100) / 1000
    for (k in 1:3) {
        s = sprintf("elv%d", k)
        df = df %>% mutate(!!s := x^k)
    }

    # Basis functions for year.
    x = (df$year - 2010) / 100
    for (k in 1:3) {
        s = sprintf("year%d", k)
        df = df %>% mutate(!!s := x^k)
    }

    return(df)
}

# Build formula terms for each model.
a = paste(sprintf("sin_day_%d", 1:4), collapse=" + ")
b = paste(sprintf("cos_day_%d", 1:4), collapse=" + ")
dayterms = sprintf("%s + %s", a, b)
a = paste(sprintf("sin_lon_%d", 1:4), collapse=" + ")
b = paste(sprintf("cos_lon_%d", 1:4), collapse=" + ")
lonterms = sprintf("%s + %s", a, b)
latterms = paste(sprintf("lat%d", 1:3), collapse=" + ")
elvterms = paste(sprintf("elv%d", 1:3), collapse=" + ")

# Get covariate terms to adjust for everything except the response.
get_covariate_terms = function(response) {
    terms = list()
    tx = list(dayterms, lonterms, latterms, elvterms)
    na = c("day", "decimalLongitude", "decimalLatitude", "elevation")
    for (j in 1:4) {
        if (response != na[j]) {
            terms[[1+length(terms)]] = tx[[j]]
        }
    }
    return(paste(terms, collapse=" + "))
}

# Fit models with year as a linear term, and additive main effects
# for all non-response variables.
fit_linear = function(response) {

    terms = get_covariate_terms(response)

    fml0 = as.formula(sprintf("%s ~ %s", response, terms))
    fml1 = as.formula(sprintf("%s ~ year + %s", response, terms))

    # Fit a model for each species.
    pva = df %>% group_by(scientificName) %>% group_modify(~{
        dv = setbasis(.x)
        m0 = lm(fml0, data=dv)
        m1 = lm(fml1, data=dv)
        pv = anova(m0, m1)[["Pr(>F)"]]
        dd = data.frame(pv=pv[2], coef=coef(m1)[2])
        return(dd)
    })

    return(pva)
}

# Fit a nonlinear model for the response variable (should be one of
# latitude, longitude, or elevation) predicted by other control
# variables, with and without year.  If year is included, it
# is modeled nonlinearly with a 5-degree spline.
fit_nonlin = function(response) {

    terms = get_covariate_terms(response)

    # Compare these two models.
    fml0 = as.formula(sprintf("%s ~ %s", response, terms))
    fml1 = as.formula(sprintf("%s ~ bs(year, 5) + %s", response, terms))

    # Fit a model for each species.
    pvb = df %>% group_by(scientificName) %>% group_modify(~{
        dv = setbasis(.x)
        m0 = lm(fml0, data=dv)
        m1 = lm(fml1, data=dv)
        pv = anova(m0, m1)[["Pr(>F)"]]
        dd = data.frame(pv=pv[2], coef=0)
        return(dd)
    })

    return(pvb)
}

make_plots = function(pva, plot_slopes, title) {
    # Quantile plot of p-values
    pvas = sort(pva$pv)
    n = length(pvas)
    dp = data.frame(x=seq(1/n, 1-1/n, length.out=n), y=pvas)
    plt = ggplot(aes(x=x, y=y), data=dp) + geom_line()
    plt = plt + geom_abline(intercept=0, slope=1)
    plt = plt + labs(x="Expected p-value", y="Observed p-value")
    plt = plt + ggtitle(title)
    print(plt)

    if (plot_slopes) {
        # Scatterplot year slopes against p-values
        plt = ggplot(aes(x=pva$pv, pva$coef), data=pva) + geom_point()
        plt = plt + labs(x="p-value", y="Year slope")
        plt = plt + ggtitle(title)
        print(plt)
    }
}

pdf("bmodel_r.pdf")

responses = c("decimalLatitude", "decimalLongitude", "elevation")
titles = c("Latitude", "Longitude", "Elevation")

for (j in 1:3) {
    pva = fit_linear(responses[j])
    make_plots(pva, TRUE, sprintf("%s (linear year)", titles[j]))
    pvb = fit_nonlin(responses[j])
    make_plots(pvb, FALSE, sprintf("%s (nonlinear year)", titles[j]))
}

dev.off()
