library(dplyr)
library(ggplot2)
library(readr)
library(dr)
source("read.R")

dx = df %>% dplyr::select(BPXSY1, RIAGENDR, RIDAGEYR, BMXWT, BMXHT, BMXBMI, BMXLEG, BMXARML, BMXARMC, BMXWAIST, BMXHIP)
dx = dx[complete.cases(dx),]

dx$RIAGENDRx = recode(dx$RIAGENDR, "F"=1, "M"=-1)

for (m in names(dx)) {
    if (is.numeric(dx[[m]])) {
        dx[[m]] = scale(dx[[m]], scale=F)
    }
}

m = dr(BPXSY1 ~ RIAGENDRx + RIDAGEYR + BMXWT + BMXHT + BMXBMI + BMXLEG + BMXARML + BMXARMC + BMXWAIST + BMXHIP, dx)

# Based on these tests we focus on the first two directions
print(dr.test(m))

scores = m$x %*% m$evectors

pdf("sbp_dr_r.pdf")

# Stratify on the j'th score, plot the mean of SBP with respect to the
# k'th score.
plotstrat = function(j, k) {
    dp = data.frame(strat=scores[,j], x=scores[,k], y=dx$BPXSY1)
    dp[,"stratq"] = ntile(dp[,"strat"], 5)
    xx = seq(min(scores[,k]), max(scores[,k]), length.out=100)
    dg = dp %>% group_by(stratq) %>% group_modify(~ {
        m = lowess(.x$x, .x$y)
        f = approxfun(m)
        data.frame(score=xx, sbp=f(xx))
    })
    dg$stratq = as.factor(dg$stratq)
    plt = ggplot(aes(x=score, y=sbp, color=stratq), data=dg)
    plt = plt + geom_line()
    plt = plt + labs(x=sprintf("Score %d", k), y="SBP (centered)")
    print(plt)
}

plotstrat(2, 1)
plotstrat(1, 2)

# Make a contour plot showing the mean SBP for each
# combination of score values (for the first two factors).
dp = data.frame(s1=scores[,1], s2=scores[,2], y=dx$BPXSY1)
rr = loess(y ~ s1 + s2, dp)
s1 = seq(min(scores[,1]), max(scores[,1]), length.out=100)
s2 = seq(min(scores[,2]), max(scores[,2]), length.out=100)
x1 = kronecker(s1, array(1, 100))
x2 = kronecker(array(1, 100), s2)
xx = data.frame(s1=x1, s2=x2)
xx$yy = predict(rr, xx)
plt = ggplot(aes(s1, s2, z=yy), data=xx) + geom_contour_filled()
print(plt)

# Plot the DV (SBP) against each score, or plot each score against every covariate.
dx$sex2 = recode(dx$RIAGENDR, "F"="EF", "M"="EM")
for (j in 1:2) {
    for (x in names(dx)) {
        if (x != "BPXSY1") {
            dp = data.frame(x=dx[[x]], y=scores[,j], sex=as.factor(dx$RIAGENDR),
                            sex2=as.factor(dx$sex2))
            plt = ggplot(aes(x=x, y=y, color=sex), data=dp)
            plt = plt + geom_point(alpha=0.2, fill=NA)
            plt = plt + geom_smooth(se=FALSE)
            plt = plt + labs(x=x, y=sprintf("Score %d", j))
        } else {
            dp = data.frame(x=scores[,j], y=dx[[x]])
            plt = ggplot(aes(x=x, y=y), data=dp) + geom_point(alpha=0.2, fill=NA)
            plt = plt + geom_smooth(se=FALSE)
            plt = plt + labs(y=x, x=sprintf("Score %d", j))
        }
        print(plt)
    }
}

dev.off()
