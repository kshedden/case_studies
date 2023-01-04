library(splines)
library(ggplot2)
library(readr)
source("read.R")

pdf("sbp_lm_r.pdf")

# Data for plotting
dp = df[1:100,]
dp$RIDAGEYR = seq(18, 80, length.out=100)
dp$BPXSY1 = 0
dp$RIDRETH1 = "OH"

# Plot predicted SBP by sex
plot1 = function(mm, ii, dbands=FALSE) {
    df = data.frame()
    sig = summary(mm)$sigma
    for (sex in c("F", "M")) {
        dp[,"RIAGENDR"] = sex
        dp[,"BMXBMI"] = 25
        yh = predict(mm, dp)
        dx = data.frame(SBP=yh, sex=sex, age=dp[,"RIDAGEYR"])
        dx[,"SBP1"] = dx[,"SBP"] - sig
        dx[,"SBP2"] = dx[,"SBP"] + sig
        df = rbind(df, dx)
    }
    plt = ggplot(aes(x=RIDAGEYR, y=SBP, color=sex, by=sex), data=df) + geom_line()
    if (dbands) {
        for (sex in c("F", "M")) {
            plt = plt + geom_line(aes(y=SBP1), linetype="dotted")
            plt = plt + geom_line(aes(y=SBP2), linetype="dotted")
        }
    }
    print(plt)
}

# Plot predicted SBP by sex and BMI (25 versus 30)
plot2 = function(mm, ii) {
    df = data.frame()
    for (sex in c("F", "M")) {
        for (bmi in c(25, 30)) {
            dp[,"RIAGENDR"] = sex
            dp[,"BMXBMI"] = bmi
            yh = predict(mm, dp)
            dx = data.frame(SBP=yh, bmi=bmi, sex=sex, age=dp[,"RIDAGEYR"])
            df = rbind(df, dx)
        }
    }
    plt = ggplot(aes(x=RIDAGEYR, y=SBP, color=interaction(sex, bmi), by=interaction(sex, bmi)), data=df) + geom_line()
    print(plt)
}

# Compare females and males at fixed BMI, with confidence bands
plot3 = function(mm, ii) {
    dp[,"BMXBMI"] = 25
    yy = list()
    xm = list()
    for (sex in c("F", "M")) {
        dp[,"RIAGENDR"] = sex
        yy[[length(yy)+1]] = predict(mm, dp)
        xm[[length(xm)+1]] = model.matrix(mm, data=dp)
    }
    vc = vcov(mm)
    xd = xm[[1]] - xm[[2]]
    yd = yy[[1]] - yy[[2]]
    se = sqrt(diag(xd %*% vc %*% t(xd)))
    df = data.frame(yd=yd, se=se, age=dp[,"RIDAGEYR"])
    df[,"y1"] = df[,"yd"] + 2*df[,"se"]
    df[,"y2"] = df[,"yd"] - 2*df[,"se"]
    plt = ggplot(aes(x=RIDAGEYR, y=yd), data=df)
    plt = plt + geom_ribbon(aes(x=RIDAGEYR, ymin=y1, ymax=y2), fill="grey70")
    plt = plt + geom_line()
    plt = plt + labs(x="Age", y="Female SBP minus male SBP")
    print(plt)
}

# Compare BMI 25 to BMI 30, for females only
plot4 = function(mm, ii) {
    dp[,"RIAGENDR"] = "F"
    yy = list()
    xm = list()
    for (bmi in c(30, 25)) {
        dp[,"BMXBMI"] = bmi
        yy[[length(yy)+1]] = predict(mm, dp)
        xm[[length(xm)+1]] = model.matrix(mm, data=dp)
    }
    vc = vcov(mm)
    xd = xm[[1]] - xm[[2]]
    yd = yy[[1]] - yy[[2]]
    se = sqrt(diag(xd %*% vc %*% t(xd)))
    df = data.frame(yd=yd, se=se, age=dp[,"RIDAGEYR"])
    df[,"y1"] = df[,"yd"] + 2*df[,"se"]
    df[,"y2"] = df[,"yd"] - 2*df[,"se"]
    plt = ggplot(aes(x=RIDAGEYR, y=yd), data=df)
    plt = plt + geom_ribbon(aes(x=RIDAGEYR, ymin=y1, ymax=y2), fill="grey70")
    plt = plt + geom_line()
    plt = plt + labs(x="Age", y="SBP at BMI 30 vs. 25")
    print(plt)
}

all_plots = function(mm, ii) {
    plot1(mm, ii)
    plot1(mm, ii, dbands=TRUE)
    plot2(mm, ii)
    plot3(mm, ii)
    plot4(mm, ii)
}

# Very basic model
f0 = as.formula("BPXSY1 ~ RIDAGEYR + RIAGENDR + BMXBMI")
m0 = lm(f0, df)
all_plots(m0, 0)

# Allow age slopes to differ by sex
f1 = as.formula("BPXSY1 ~ RIDAGEYR * RIAGENDR + BMXBMI")
m1 = lm(f1, df)
all_plots(m1, 1)

# Allow BMI slopes to differ by sex
f2 = as.formula("BPXSY1 ~ RIDAGEYR + RIAGENDR * BMXBMI")
m2 = lm(f2, df)
all_plots(m2, 2)

# Allow age and BMI slopes to differ by sex
f3 = as.formula("BPXSY1 ~ (RIDAGEYR + BMXBMI) * RIAGENDR")
m3 = lm(f3, df)
all_plots(m3, 3)

# Full interactions among age, BMI, and sex
f4 = as.formula("BPXSY1 ~ RIDAGEYR * BMXBMI * RIAGENDR")
m4 = lm(f4, df)
all_plots(m4, 4)

# Basic additive model with splines
f5 = as.formula("BPXSY1 ~ bs(RIDAGEYR, 5) + RIAGENDR + BMXBMI")
m5 = lm(f5, df)
all_plots(m5, 5)

# Allow age trends to differ by sex
f6 = as.formula("BPXSY1 ~ bs(RIDAGEYR, 5) * RIAGENDR + BMXBMI")
m6 = lm(f6, df)
all_plots(m6, 6)

# Allow BMI trends to differ by sex
f7 = as.formula("BPXSY1 ~ bs(RIDAGEYR, 5) + RIAGENDR * BMXBMI")
m7 = lm(f7, df)
all_plots(m7, 7)

# Allow age and BMI trends to differ by sex
f8 = as.formula("BPXSY1 ~ (bs(RIDAGEYR, 5) + BMXBMI) * RIAGENDR")
m8 = lm(f8, df)
all_plots(m8, 8)

# Full interactions among nonlinear age, BMI, and sex
f9 = as.formula("BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR")
m9 = lm(f9, df)
all_plots(m9, 9)

# Full interactions among nonlinear age, BMI, and sex with additive control for ethnicity
f10 = as.formula("BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR + RIDRETH1")
m10 = lm(f10, df)
all_plots(m10, 10)

# Full interactions among nonlinear age, BMI, and sex with ethnicity x sex interactions
f11 = as.formula("BPXSY1 ~ (bs(RIDAGEYR, 5) * BMXBMI + RIDRETH1) * RIAGENDR")
m11 = lm(f11, df)
all_plots(m11, 11)

# Full interactions among nonlinear age, BMI, and sex, and between sex and ethnicity, and
# between linear age and ethnicity.
f12 = as.formula("BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR + (RIAGENDR + RIDAGEYR) * RIDRETH1")
m12 = lm(f12, df)
all_plots(m12, 12)

# Full interactions among nonlinear age, BMI, and sex, and between sex and ethnicity, and
# between nonlinear age and ethnicity.
f13 = as.formula("BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR + (RIAGENDR + bs(RIDAGEYR, 5)) * RIDRETH1")
m13 = lm(f13, df)
all_plots(m13, 13)

# Full interactions among nonlinear age, BMI, sex, and ethnicity.
f14 = as.formula("BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR * RIDRETH1")
m14 = lm(f14, df)
all_plots(m14, 14)

# Calculate AIC for all models
mm = list(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
maic = array(0, length(mm))
for (i in 1:length(mm)) {
    maic[i] = AIC(mm[[i]])
}
maic = maic - min(maic)

dev.off()
