using LombScargle, UnicodePlots

include("read.jl")

# The time values are in units of seconds, frequency is 1/period.  So a freuency of 1e-6 is a period of
# 1e6, which is 1e6/(60*60) hours.
plan = LombScargle.plan(df[:, :Time], df[:, :Flux1], oversampling=5, minimum_frequency=1e-6, maximum_frequency=2*1e-4)
prgm = lombscargle(plan)

f, p = freqpower(prgm)

# Get the period of the maximum power in units of hours.
mf = findmaxfreq(prgm)
mph = (1/mf) / (60*60)

# Plot the periodogram, the horizontal axis is cycles/second.
plt = lineplot(f, p, canvas=BlockCanvas)
pw = fapinv(prgm, 0.05)
lineplot!(plt, f, pw*ones(length(p)))
println(plt)
