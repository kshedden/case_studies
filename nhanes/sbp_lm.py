import numpy as np
import patsy
import statsmodels.api as sm
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from read import df

cm = matplotlib.cm.get_cmap("tab10")
pdf = PdfPages("sbp_lm_py.pdf")

dp = df.iloc[0:50, :].copy()
dp["RIDAGEYR"] = np.linspace(18, 80, 50)
dp["RIDRETH1"] = "MA"

# Plot predicted SBP by sex at fixed BMI.
def plot1(rr, ii, dbands=False, bmi=2):
    sigma = np.sqrt(rr.scale)
    plt.figure(figsize=(7.6, 5))
    plt.clf()
    plt.axes([0.12, 0.12, 0.7, 0.8])
    plt.grid(True)
    for (jj,sex) in enumerate(["F", "M"]):
        dp["RIAGENDR"] = sex
        dp["BMXBMI"] = bmi
        yh = rr.predict(exog=dp)
        plt.plot(dp.RIDAGEYR, yh, color=cm(jj/10), label={"F": "Female", "M": "Male"}[sex])
        if dbands:
            plt.plot(dp.RIDAGEYR, yh + sigma, ":", color=cm(jj/10))
            plt.plot(dp.RIDAGEYR, yh - sigma, ":", color=cm(jj/10))
    plt.xlabel("Age", size=15)
    plt.ylabel("SBP", size=15)
    plt.title("Model %d" % ii)
    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, loc="center right")
    leg.draw_frame(False)
    pdf.savefig()

# Plot predicted SBP by sex at two different BMI levels (25 versus 30)
def plot2(rr, ii, bmis = [25, 30]):
    plt.figure(figsize=(7.6, 5))
    plt.clf()
    plt.axes([0.12, 0.12, 0.7, 0.8])
    plt.grid(True)
    for sex in ["F", "M"]:
        for bmi in bmis:
            dp["RIAGENDR"] = sex
            dp["BMXBMI"] = bmi
            yh = rr.predict(exog=dp)
            plt.plot(dp.RIDAGEYR, yh, label="%s/%.0f" %
                     ({"F": "Female", "M": "Male"}[sex], bmi))
    plt.xlabel("Age", size=15)
    plt.ylabel("SBP", size=15)
    plt.title("Model %d" % ii)
    ha, lb = plt.gca().get_legend_handles_labels()
    leg = plt.figlegend(ha, lb, loc="center right")
    leg.draw_frame(False)
    leg.set_title("Sex/BMI")
    pdf.savefig()

# Compare females and males at fixed BMI, for each age, with confidence bands
def plot3(rr, ii, bmi=25):
    dp["BMXBMI"] = bmi
    yy, xm = [], []
    for sex in ["F", "M"]:
        dp["RIAGENDR"] = sex
        xx = patsy.dmatrix(rr.model.data.design_info, dp, return_type="dataframe")
        xm.append(xx)
        y = rr.predict(exog=dp)
        yy.append(y)
    xd = xm[0] - xm[1]
    vc = np.dot(xd, np.dot(rr.cov_params(), xd.T))
    se = np.sqrt(np.diag(vc))
    yd = yy[0] - yy[1]

    plt.figure(figsize=(7.6, 5))
    plt.clf()
    plt.axes([0.12, 0.12, 0.7, 0.8])
    plt.grid(True)
    plt.fill_between(dp.RIDAGEYR, yd-2*se, yd+2*se, color="grey")
    plt.plot(dp.RIDAGEYR, yd, color="black")
    plt.xlabel("Age", size=15)
    plt.ylabel("SBP difference", size=15)
    plt.title("Model %d" % ii)
    plt.title("SBP difference based on sex (F-M) at BMI=25")
    pdf.savefig()

# Compare BMI 25 to BMI 30, for one sex only
def plot4(rr, ii, sex="F"):
    dp["RIAGENDR"] = sex
    yy, xm = [], []
    for bmi in [30, 25]:
        dp["BMXBMI"] = bmi
        xx = patsy.dmatrix(rr.model.data.design_info, dp, return_type="dataframe")
        xm.append(xx)
        y = rr.predict(exog=dp)
        yy.append(y)
    xd = xm[0] - xm[1]
    vc = np.dot(xd, np.dot(rr.cov_params(), xd.T))
    se = np.sqrt(np.diag(vc))
    yd = yy[0] - yy[1]

    plt.figure(figsize=(7.6, 5))
    plt.clf()
    plt.axes([0.12, 0.12, 0.7, 0.8])
    plt.grid(True)
    plt.fill_between(dp.RIDAGEYR, yd-2*se, yd+2*se, color="grey")
    plt.plot(dp.RIDAGEYR, yd, color="black")
    plt.xlabel("Age", size=15)
    plt.ylabel("SBP", size=15)
    plt.title("Model %d" % ii)
    plt.title("SBP difference based on BMI (30-25) for females")
    pdf.savefig()

def plot_all(rr, ii):
    plot1(rr, ii, bmi=25)
    plot1(rr, ii, dbands=True)
    plot2(rr, ii, bmis=[25, 30])
    plot3(rr, ii, bmi=25)
    plot4(rr, ii, sex="F")

# Very basic model
f0 = "BPXSY1 ~ RIDAGEYR + RIAGENDR + BMXBMI"
m0 = sm.OLS.from_formula(f0, df)
r0 = m0.fit()
plot_all(r0, 0)

# Allow age slopes to differ by sex
f1 = "BPXSY1 ~ RIDAGEYR * RIAGENDR + BMXBMI"
m1 = sm.OLS.from_formula(f1, df)
r1 = m1.fit()
plot_all(r1, 1)

# Allow BMI slopes to differ by sex
f2 = "BPXSY1 ~ RIDAGEYR + RIAGENDR * BMXBMI"
m2 = sm.OLS.from_formula(f2, df)
r2 = m2.fit()
plot_all(r2, 2)

# Allow age and BMI slopes to differ by sex
f3 = "BPXSY1 ~ (RIDAGEYR + BMXBMI) * RIAGENDR"
m3 = sm.OLS.from_formula(f3, df)
r3 = m3.fit()
plot_all(r3, 3)

# Full interactions among age, BMI, and sex
f4 = "BPXSY1 ~ RIDAGEYR * BMXBMI * RIAGENDR"
m4 = sm.OLS.from_formula(f4, df)
r4 = m4.fit()
plot_all(r4, 4)

# Check AICs
print("AIC for models 0-4:")
print([x.aic for x in (r0, r1, r2, r3, r4)])

# Basic additive model with splines
f5 = "BPXSY1 ~ bs(RIDAGEYR, 5) + RIAGENDR + BMXBMI"
m5 = sm.OLS.from_formula(f5, df)
r5 = m5.fit()
plot_all(r5, 5)

# Allow age trends to differ by sex
f6 = "BPXSY1 ~ bs(RIDAGEYR, 5) * RIAGENDR + BMXBMI"
m6 = sm.OLS.from_formula(f6, df)
r6 = m6.fit()
plot_all(r6, 6)

# Allow BMI trends to differ by sex
f7 = "BPXSY1 ~ bs(RIDAGEYR, 5) + RIAGENDR * BMXBMI"
m7 = sm.OLS.from_formula(f7, df)
r7 = m7.fit()
plot_all(r7, 7)

# Allow age and BMI trends to differ by sex
f8 = "BPXSY1 ~ (bs(RIDAGEYR, 5) + BMXBMI) * RIAGENDR"
m8 = sm.OLS.from_formula(f8, df)
r8 = m8.fit()
plot_all(r8, 8)

# Full interactions among nonlinear age, BMI, and sex
f9 = "BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR"
m9 = sm.OLS.from_formula(f9, df)
r9 = m9.fit()
plot_all(r9, 9)

# Full interactions among nonlinear age, BMI, and sex with additive control for ethnicity
f10 = "BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR + RIDRETH1"
m10 = sm.OLS.from_formula(f10, df)
r10 = m10.fit()
plot_all(r10, 10)

# Full interactions among nonlinear age, BMI, and sex with ethnicity x sex interactions
f11 = "BPXSY1 ~ (bs(RIDAGEYR, 5) * BMXBMI + RIDRETH1) * RIAGENDR"
m11 = sm.OLS.from_formula(f11, df)
r11 = m11.fit()
plot_all(r11, 11)

# Full interactions among nonlinear age, BMI, and sex, and between sex and ethnicity, and
# between linear age and ethnicity.
f12 = "BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR + (RIAGENDR + RIDAGEYR) * RIDRETH1"
m12 = sm.OLS.from_formula(f12, df)
r12 = m12.fit()
plot_all(r12, 12)

# Full interactions among nonlinear age, BMI, and sex, and between sex and ethnicity, and
# between nonlinear age and ethnicity.
f13 = "BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR + (RIAGENDR + bs(RIDAGEYR, 5)) * RIDRETH1"
m13 = sm.OLS.from_formula(f13, df)
r13 = m13.fit()
plot_all(r13, 13)

# Full interactions among nonlinear age, BMI, sex, and ethnicity.
f14 = "BPXSY1 ~ bs(RIDAGEYR, 5) * BMXBMI * RIAGENDR * RIDRETH1"
m14 = sm.OLS.from_formula(f14, df)
r14 = m14.fit()
plot_all(r14, 14)

# Check AICs
print("AIC for models 0-9:")
maic = [x.aic for x in (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14)]
maic = np.asarray(maic)
maic -= maic.min()
print(maic)

pdf.close()


# Plot the basis functions
pdf = PdfPages("splines.pdf")
df = df.sort_values(by="RIDAGEYR")
for j in [3, 5, 10]:
    y = patsy.dmatrix("0 + bs(RIDAGEYR, %d)" % j, df)
    plt.clf()
    plt.grid(True)
    plt.title("%d dimensional cubic b-spline basis" % j)
    for k in range(j):
        plt.plot(df["RIDAGEYR"], y[:, k], "-", color="blue")
    plt.xlabel("RIDAGEYR", size=15)
    plt.ylabel("Basis function value", size=15)
    pdf.savefig()

pdf.close()
