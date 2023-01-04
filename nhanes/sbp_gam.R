library(dplyr)
library(ggplot2)
library(mgcv)

pdf("sbp_gam_r.pdf")

source("read.R")

dx = df %>% dplyr::select(BPXSY1, RIAGENDR, RIDAGEYR, BMXWT, BMXHT, BMXBMI, BMXLEG, BMXARML, BMXARMC, BMXWAIST, BMXHIP)
dx = dx[complete.cases(dx),]

dx$RIAGENDRx = recode(dx$RIAGENDR, "F"=1, "M"=-1)
dx$RIAGENDR = as.factor(dx$RIAGENDR)

# Compare sexes at fixed BMI
plot1 = function(mm, bmi, ii) {
    dp = dx[1:100,]
    dp$RIDAGEYR = seq(18, 80, length.out=100)
    dp$BMXBMI = bmi

    dr = data.frame()
    for (sex in c("F", "M")) {
        dp$RIAGENDR = sex
        sbp = predict(mm, newdata=dp)
        d2 = data.frame(age=dp$RIDAGEYR, sbp=sbp, sex=sex)
        dr = rbind(dr, d2)
    }

    plt = ggplot(aes(x=age, y=sbp, color=sex, by=sex), data=dr) + geom_line()
    plt = plt + ggtitle(sprintf("Model %d", ii))
    plt = plt + labs(x="Age", y="SBP")
    print(plt)
}

# Compare BMI 30 versus 25 for one sex.
plot2 = function(mm, sex, ii) {
    dp = dx[1:100,]
    dp$RIDAGEYR = seq(18, 80, length.out=100)
    dp$RIAGENDR = sex

    dr = data.frame()
    for (bmi in c(25, 30)) {
        dp$BMXBMI = bmi
        sbp = predict(mm, newdata=dp)
        d2 = data.frame(age=dp$RIDAGEYR, sbp=sbp, bmi=dp$BMXBMI)
        dr = rbind(dr, d2)
    }
    dr$bmi = as.factor(dr$bmi)

    plt = ggplot(aes(x=age, y=sbp, color=bmi, by=bmi), data=dr) + geom_line()
    plt = plt + ggtitle(sprintf("Model %d (%s only)", ii, sex))
    plt = plt + labs(x="Age", y="SBP")
    print(plt)
}

plot_all = function(mm, ii) {
    plot(mm, se=TRUE)
    plot1(mm, 25, ii)
    plot1(mm, 30, ii)
    plot2(mm, "F", ii)
    plot2(mm, "M", ii)
}

# A simple additive model
f0 = as.formula("BPXSY1 ~ s(RIDAGEYR) + s(BMXBMI) + RIAGENDR")
m0 = gam(f0, data=dx)
plot_all(m0, 0)

# Allow the age effects to differ by sex
f1 = as.formula("BPXSY1 ~ s(RIDAGEYR, by=RIAGENDR) + BMXBMI")
m1 = gam(f1, data=dx)
plot_all(m1, 1)

# Allow the BMI effects to differ by sex
f2 = as.formula("BPXSY1 ~ RIDAGEYR + s(BMXBMI, by=RIAGENDR)")
m2 = gam(f2, data=dx)
plot_all(m2, 2)

# Allow both the age and BMI effects to differ by sex
f3 = as.formula("BPXSY1 ~ s(RIDAGEYR, by=RIAGENDR) + s(BMXBMI, by=RIAGENDR)")
m3 = gam(f3, data=dx)
plot_all(m3, 3)

# Allow age and BMI to relate non-additively to SBP
f4 = as.formula("BPXSY1 ~ te(RIDAGEYR, BMXBMI) + s(RIDAGEYR, by=RIAGENDR) + s(BMXBMI, by=RIAGENDR)")
m4 = gam(f4, data=dx)
plot_all(m4, 4)

# Allow the non-additive relationship between age and BMI to differ by sex
f5 = as.formula("BPXSY1 ~ te(RIDAGEYR, BMXBMI, by=RIAGENDR) + s(RIDAGEYR, by=RIAGENDR) + s(BMXBMI, by=RIAGENDR)")
m5 = gam(f5, data=dx)
plot_all(m5, 5)

# Add a completely random covariate to see how well the inference works
dx$E = rnorm(dim(dx)[1])
f6 = as.formula("BPXSY1 ~ te(RIDAGEYR, BMXBMI, by=RIAGENDR) + s(RIDAGEYR, by=RIAGENDR) + s(BMXBMI, by=RIAGENDR) + s(E)")
m6 = gam(f6, data=dx)

# Compare each successive pair of models.
ml = list(m0, m1, m2, m3, m4, m5, m6)
for (i in 1:(length(ml)-1)) {
    aa = anova(ml[[i]], ml[[i+1]], test="F")
    cat(sprintf("Model %d versus %d: %f\n", i-1, i, aa[["Pr(>F)"]][2]))
}

dev.off()
