using CSV
using DataFrames
using Statistics
using Distributions
using CairoMakie
using Printf
using Dates

include("generalizedpareto.jl")

# https://link.springer.com/article/10.1007/s42452-020-03199-8

rm("plots", recursive = true, force = true)
mkdir("plots")

target_dir = "/home/kshedden/data/Teaching/precip"

fname = "USW00094847.csv" # Detroit
#fname = "USW00012839.csv" # Miami

df = open(joinpath(target_dir, fname*".gz")) do io
    CSV.read(io, DataFrame)
end

df = df[:, [:DATE, :PRCP]]
df = df[completecases(df), :]
df = disallowmissing(df)

# Convert precipitation to millimiters
df[:, :PRCP] ./= 10

# Use this threshold for calculating exceedances
thresh = 5.0

# Median annual maximum
df[:, :year] = [year(x) for x in df[:, :DATE]]
annmax = median(combine(groupby(df, :year), :PRCP=>maximum)[:, 2])

# Return values x, p such that the slope of p on x
# estimates the shape parameter (tail index) of
# a Pareto distribution, or the rate parameter
# of an Exponential distribution.  The elements of
# 'z' should follow a distribution whose upper tail is
# either Pareto distribution or Exponential.
# The upper p0 fraction of the data in z are used
# to produce x, p.  The values in x are the order
# statistics of z in the Exponential case, and the log
# order statistics of z in the Pareto case. The values in
# p are derived from probability points.
function tail_shape(z; p0=0.1, family=:pareto)
    if !(family in [:exponential, :pareto])
        error("Unknown family $(family)")
    end
    p = 1 - p0
    z = sort(z)
    n = length(z)
    m = Int(round(p*n))
    x = z[m:end]
    if family == :pareto
        x = log.(x)
    end
    p = log.(1 .- (m:n) ./ (n+1))
    return x, p
end

function check_pareto(alpha; n=10000, p0=0.1)
    u = rand(Uniform(), n)
    z = u.^(-1/alpha) # Pareto values
    x, p = tail_shape(z; p0=p0, family=:pareto)
    alpha_hat = -cov(p, x) / var(x)
    println(alpha, " ", alpha_hat)
end

function check_exponential(mu; n=10000, p0=0.1)
    z = -mu*log.(rand(n)) # Exponential values
    x, p = tail_shape(z; p0=p0, family=:exponential)
    rate_hat = -cov(p, x) / var(x)
    mu_hat = 1 / rate_hat
    println(mu, " ", mu_hat)
end

# Use least squares regression in the tail of a quantile plot
# to estimate the shape parameter.
function fit_tail_reg(x, fig, ax; p0=0.99, family=:pareto)

    x, p = tail_shape(x; p0=p0, family=family)
    lines!(ax, x, p; color=:orange)

    # Estimate the tail index using a least squares fit to the order
    # statistics.
    alpha_hat = -cov(p, x) / var(x)
    icept = mean(p) + alpha_hat*mean(x)

    # The coordinates of the best-fit line
    xx = extrema(x)
    xx = [xx[1], xx[2]]
    yy = icept .- alpha_hat.*xx

    lines!(ax, xx, yy; color=:purple)

    return icept, alpha_hat
end

function plot_tails(z, p0, thresh, family, ifig)

    n = length(z)

    # Select only extreme values and translate back to the origin
    z = z[z .>= thresh]
    z .-= thresh

    # The number of selected observations
    m = Int(round(p0*n))

    xlabel = family == :pareto ? "log Q(p)" : "Q(p)"

    fig = Figure()
    ax = Axis(fig[1,1], ylabel="log(1-p)", xlabel=xlabel,
              xlabelsize=18, ylabelsize=18)
    icept, alpha = fit_tail_reg(z, fig, ax; p0=p0, family=family)
    ax.title = @sprintf("%s model, threshold=%.1f, top %.1f%% (n=%d), alpha=%.3f",
                        titlecase(String(family)), thresh, 100*p0, m, alpha)
    save(@sprintf("plots/%03d.pdf", ifig), fig)
    return ifig + 1
end

function mobs_return(z, mr; thresh=thresh, family=:exponential, gp=nothing)

    n = length(z)

    # Select only extreme values and translate back to the origin
    ix = z .>= thresh
    q = sum(ix) / n # proportion of values exceeding the threshold
    z = z[ix]
    z .-= thresh

    pr = 1 .- 1 ./ (q.*mr)

    if family == :exponential
        mn = mean(z)
        m0 = thresh .- mn*log.(1 .- pr)
    elseif family == :generalizedpareto
        m0 = thresh .+ quantile(gp, pr)
    else
        error("!!")
    end

    return m0
end

function hill(z; k=100)

    z = sort(z)
    z = log.(z)

    h = 0.0
    for j in 1:k-1
        h += z[end + 1 - j] - z[end - k + 1]
    end
    h /= k - 1

    return 1 / h
end

function plot_hill(z, ifig)

    kv = 20:5:500
    ta = [hill(z; k=k) for k in kv]

    fig = Figure()
    ax = Axis(fig[1,1], ylabel="Tail index estimate", xlabel="k",
              xlabelsize=18, ylabelsize=18)
    ax.title = "Hill estimate of the tail index"
    lines!(ax, kv, ta)
    save(@sprintf("plots/%03d.pdf", ifig), fig)
    return ifig + 1
end

function fit_gpar(z, alpha, ifig)

    # Select only extreme values and translate back to the origin
    z = z[z .>= thresh]
    z .-= thresh
    n = length(z)

    xi = 1 / alpha
    sigma = median(z) * xi / (2^xi - 1)

    gp = GeneralizedPareto(sigma, xi)

    pp = (1:n) / (n + 1)
    qq = [quantile(gp, p) for p in pp]
    z = sort(z)

    fig = Figure()
    ax = Axis(fig[1,1], ylabel="Order statistics", xlabel="GP quantiles",
              xlabelsize=18, ylabelsize=18)
    ax.title = @sprintf("alpha=%.2f", alpha)
    lines!(ax, qq, z)
    save(@sprintf("plots/%03d.pdf", ifig), fig)
    ifig += 1

    fig = Figure()
    ax = Axis(fig[1,1], ylabel="Sqrt order statistics", xlabel="Sqrt GP quantiles",
              xlabelsize=18, ylabelsize=18)
    ax.title = @sprintf("alpha=%.2f", alpha)
    lines!(ax, sqrt.(qq), sqrt.(z))
    save(@sprintf("plots/%03d.pdf", ifig), fig)
    ifig += 1

    return gp, ifig
end

function mle_analysis(z, ifig)
    z = z[z .> thresh] .- thresh
    mle = fit(GeneralizedPareto, z, μ=0.; improved=false)

    n = length(z)
    pp = (1:n) ./ (n+1)
    qq = [quantile(mle, p) for p in pp]
    z = sort(z)

    fig = Figure()
    ax = Axis(fig[1,1], ylabel="Order statistics", xlabel="GP quantiles (MLE)",
              xlabelsize=18, ylabelsize=18)
    ax.title = @sprintf("MLE: %s", mle)
    lines!(ax, qq, z)
    save(@sprintf("plots/%03d.pdf", ifig), fig)
    ifig += 1

    return mle, ifig
end

function main(ifig)

    # Quantile plots
    for family = [:pareto, :exponential]
        for p0 in [0.5, 0.1, 0.05, 0.01]
            ifig = plot_tails(df[:, :PRCP], p0, thresh, family, ifig)
        end
    end

    mle, ifig = mle_analysis(df[:, :PRCP], ifig)

    ifig = plot_hill(df[:, :PRCP], ifig)

    gp = []
    for alpha in [3, 4, 5, 6]
        gp1, ifig = fit_gpar(df[:, :PRCP], alpha, ifig)
        push!(gp, gp1)
    end

    yr = [1, 10, 100, 500, 1000]

    cfg = vcat((:exponential,nothing), [(:generalizedpareto,g) for g in gp], (:generalizedpareto, mle))

    for (f,g) in cfg
        println("$(f) $(g)")
        mr = mobs_return(df[:, :PRCP], 365 .* yr; family=f, thresh=thresh, gp=g)
        rr = DataFrame(Years=yr, MR=mr)
        println(rr)
    end

    return ifig
end

ifig = 0
ifig = main(ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -r300 -sOutputFile=tails.pdf $f`
run(c)

