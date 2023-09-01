using CSV
using DataFrames
using Statistics
using Distributions
using CairoMakie
using LinearAlgebra
using Printf
using Dates

# A paper using extreme value techniques to study rainfall in Brazil:
# https://link.springer.com/article/10.1007/s42452-020-03199-8

rm("plots", recursive = true, force = true)
mkdir("plots")

target_dir = "/home/kshedden/data/Teaching/precip"

fname = "USW00094847.csv" # Detroit
#fname = "USW00012839.csv" # Miami

# Estimate the parameters of a generalized Pareto distribution
# using the empirical Bayes method of Zhang and Stephens.
# https://www.jstor.org/stable/pdf/40586625.pdf
function gp_estimate(z)

    n = length(z)
    xstar = quantile(z, 0.25)
    m = ceil(20 + sqrt(n))
    xmax = maximum(z)

    tgrid = 1/xmax .+ (1 .- sqrt.(m ./ ((1:m) .- 0.5))) / (3 * xstar)

    function profile(theta)
        k = -mean(log, 1 .- theta*z)
        return n*(log(theta/k) + k - 1)
    end

    ltg = [profile(t) for t in tgrid]
    ltg .-= maximum(ltg)
    Ltg = exp.(ltg)
    Ltg ./= sum(Ltg)
    theta_hat = dot(Ltg, tgrid)
    k_hat = -mean(log, 1 .- theta_hat*z)
    sigma_hat = k_hat / theta_hat

    return GeneralizedPareto(sigma_hat, -k_hat)
end

df = open(joinpath(target_dir, fname*".gz")) do io
    CSV.read(io, DataFrame)
end

df = df[:, [:DATE, :PRCP]]
df = df[completecases(df), :]
df = disallowmissing(df)

# Convert precipitation to millimeters
df[:, :PRCP] ./= 10

# Use this threshold for calculating exceedances
thresh = 5.0

# Median annual maximum
df[:, :year] = [year(x) for x in df[:, :DATE]]
annmax = median(combine(groupby(df, :year), :PRCP=>maximum)[:, 2])

# Returns values x, p such that the slope of p on x
# estimates the shape parameter (tail index) of
# a distribution with power-law tails, or the rate parameter
# of a distribution with exponential tails.
# The upper p0 fraction of the data in z are used
# to produce x, p.  The values in x are the order
# statistics of z in the Exponential case, and the log
# order statistics of z in the Pareto case. The values in
# p are derived from probability points.
function tail_shape(z; p0=0.1, family=:powerlaw)
    if !(family in [:exponential, :powerlaw])
        error("Unknown family $(family)")
    end
    p = 1 - p0
    z = sort(z)
    n = length(z)
    m = Int(round(p*n))
    x = z[m:end]
    if family == :powerlaw
        x = log.(x)
    end
    p = log.(1 .- (m:n) ./ (n+1))
    return x, p
end

# Use least squares regression in the tail of a quantile plot
# to estimate the shape parameter.
function fit_tail_reg(x, fig, ax; p0=0.99, family=:powerlaw)

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

    # Exceedances
    z = z[z .>= thresh]
    z .-= thresh

    # The number of selected observations
    m = Int(round(p0*n))

    xlabel = family == :powerlaw ? "log Q(p)" : "Q(p)"

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

# Use the upper k order statistics to estimate
# the tail index of a heavy-tailed distribution
# using the Hill estimator.
function hill(z; k=300)

    z = sort(z)
    z = log.(z)

    h = 0.0
    for j in 1:k-1
        h += z[end + 1 - j] - z[end + 1 - k]
    end
    h /= k - 1

    return 1 / h
end

# Plot the Hill estimate of the tail index for a range
# of values of the tuning parameter k.
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

    # Exceedances
    z = z[z .>= thresh]
    z .-= thresh
    n = length(z)

    # Use quantile matching at the median to estimate
    # the scale parameter.
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

function eb_analysis(z, ifig)

    # Exceedances
    z = z[z .> thresh] .- thresh

    # Empirical Bayes estimate of Zhang and Stephens.
    eb = gp_estimate(z)

    n = length(z)
    pp = (1:n) ./ (n+1)
    qq = [quantile(eb, p) for p in pp]
    z = sort(z)

    fig = Figure()
    ax = Axis(fig[1,1], ylabel="Order statistics", xlabel="GP quantiles (EB)",
              xlabelsize=18, ylabelsize=18)
    ax.title = @sprintf("EB: %s", eb)
    lines!(ax, qq, z)
    save(@sprintf("plots/%03d.pdf", ifig), fig)
    ifig += 1

    return eb, ifig
end

function main(ifig)

    # Quantile plots
    for family = [:powerlaw, :exponential]
        for p0 in [0.5, 0.1, 0.05, 0.01]
            ifig = plot_tails(df[:, :PRCP], p0, thresh, family, ifig)
        end
    end

    eb, ifig = eb_analysis(df[:, :PRCP], ifig)

    ifig = plot_hill(df[:, :PRCP], ifig)

    gp = []
    for alpha in [3, 4, 5, 6]
        gp1, ifig = fit_gpar(df[:, :PRCP], alpha, ifig)
        push!(gp, gp1)
    end

    yr = [1, 10, 100, 500, 1000]

    cfg = vcat((:exponential,nothing), [(:generalizedpareto,g) for g in gp], (:generalizedpareto, eb))

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
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -r300 -sOutputFile=tails_julia.pdf $f`
run(c)

