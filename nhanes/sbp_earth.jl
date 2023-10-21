# # National Health and Nutrition Examination Survey (NHANES)

# Here we demonstrate EARTH by building regression models
# to study the variation in systolic blood pressure (SBP)
# in terms of age, BMI, sex, and ethnicity.

# See the [NHANES website](https://wwwn.cdc.gov/nchs/nhanes) for more
# information about these data.

using CategoricalArrays
using CSV
using DataFrames
using Earth
using CairoMakie
using Printf

rm("plots", recursive = true, force = true)
mkdir("plots")

function get_data()

    pa = "/home/kshedden/data/Teaching/nhanes/2017-2018"
    demog = open(joinpath(pa, "DEMO_J.csv.gz")) do io
        CSV.read(io, DataFrame)
    end

    bmx = open(joinpath(pa, "BMX_J.csv.gz")) do io
        CSV.read(io, DataFrame)
    end

    bpx = open(joinpath(pa, "BPX_J.csv.gz")) do io
        CSV.read(io, DataFrame)
    end

    demog = demog[:, [:SEQN, :RIAGENDR, :RIDAGEYR, :RIDRETH1]]
    bmx = bmx[:, [:SEQN, :BMXBMI]]
    bpx = bpx[:, [:SEQN, :BPXSY1]]

    da = leftjoin(demog, bmx, on=:SEQN)
    da = leftjoin(da, bpx, on=:SEQN)

    da = da[completecases(da), :]
    da = disallowmissing(da)
    da[!, :RIAGENDR] = replace(da[:, :RIAGENDR], 1=>"Male", 2=>"Female") .|> String
    da[!, :RIDRETH1] = replace(da[:, :RIDRETH1], 1=>"MA", 2=>"OH", 3=>"NHW", 4=>"NHB", 5=>"OR") .|> String
    da = filter(r->r.RIDAGEYR >= 18, da)
    return da
end

da = get_data()

# To use categorical variables in Earth they must be
# explicitly typed as CategoricalArray.
da[!, :RIDRETH1] = CategoricalArray(da[:, :RIDRETH1]);
da[!, :RIAGENDR] = CategoricalArray(da[:, :RIAGENDR]);

# Define the response variable as a float vector:
y = da[:, :BPXSY1];

# Construct the covariates as a named tuple:
X = (RIDAGEYR=da[:, :RIDAGEYR], BMXBMI=da[:, :BMXBMI], RIAGENDR=da[:, :RIAGENDR], RIDRETH1=da[:, :RIDRETH1]);

# Fit an additive model, limiting the order and degree of each
# term to 1.  Each term only involves a single covariate.
cfg = EarthConfig(; maxorder=1, maxdegree=1)
m1 = fit(EarthModel, X, y; config=cfg, verbose=true)

# Fit another model that allows nonlinear main effects and two-way
# interactions.
cfg = EarthConfig(; maxorder=2, maxdegree=1)
m2 = fit(EarthModel, X, y; config=cfg, verbose=true)

# Plot the fitted mean blood pressure by age at fixed levels of BMI and race,
# for females and for males.
function sbp_by_age(m, ifig; bmi=25, eth="NHB")
    dp = da[1:100, [:RIDAGEYR, :BMXBMI, :RIAGENDR, :RIDRETH1]]
    dp[1:50, :RIDAGEYR] = range(18, 80, 50)
    dp[51:100, :RIDAGEYR] = range(18, 80, 50)
    dp[1:50, :RIAGENDR] .= "Female"
    dp[51:100, :RIAGENDR] .= "Male"
    dp[:, :BMXBMI] .= bmi
    dp[:, :RIDRETH1] .= eth

    yh = predict(m, dp)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Age", ylabel="SBP",
              xlabelsize=18, ylabelsize=18, title="")
    lines!(ax, dp[1:50, :RIDAGEYR], yh[1:50], label="Female")
    lines!(ax, dp[51:100, :RIDAGEYR], yh[51:100], label="Male")
    axislegend()
    save(@sprintf("plots/%03d.pdf", ifig), fig)

    return ifig + 1
end

# Plot the fitted mean blood pressure by BMI at fixed levels of age and race,
# for females and for males.
function sbp_by_bmi(m, ifig; sex="Female", age=50, eth="NHB")
    dp = da[1:100, [:RIDAGEYR, :BMXBMI, :RIAGENDR, :RIDRETH1]]
    dp[1:50, :BMXBMI] = range(18, 40, 50)
    dp[51:100, :BMXBMI] = range(18, 40, 50)
    dp[1:50, :RIAGENDR] .= "Female"
    dp[51:100, :RIAGENDR] .= "Male"
    dp[:, :RIDAGEYR] .= age
    dp[:, :RIDRETH1] .= eth

    yh = predict(m, dp)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="BMI", ylabel="SBP",
              xlabelsize=18, ylabelsize=18, title="")
    lines!(ax, dp[1:50, :BMXBMI], yh[1:50], label="Female")
    lines!(ax, dp[51:100, :BMXBMI], yh[51:100], label="Male")
    axislegend()
    save(@sprintf("plots/%03d.pdf", ifig), fig)

    return ifig + 1
end

# The plot below shows the estimated conditional mean blood
# pressure values for non-hispanic black females, at three
# levels of BMI.

ifig = 0
ifig = sbp_by_age(m1, ifig)
ifig = sbp_by_age(m2, ifig)

ifig = sbp_by_bmi(m1, ifig)
ifig = sbp_by_bmi(m2, ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -r300 -sOutputFile=sbp_earth_jl.pdf $f`
run(c)
