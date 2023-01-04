using DataFrames, Statistics, LinearAlgebra, Loess
using UnicodePlots

include("read.jl")

# Create a matrix of observed variables that describe
# the location and time at which each profile was
# obtained.
n = length(lat)
Y = zeros(n, 3)
Y[:, 1] = lat
Y[:, 2] = lon
Y[:, 3] = date

# Get the PC's of the profiles.
tempc = copy(temp)
tempmean = mean(temp, dims=2)[:]
for j in 1:size(temp,1)
    tempc[j, :] = temp[j, :] .- tempmean[j]
end
cc = cov(tempc')
eg = eigen(Symmetric(cc))

# Reorder the PC's so that the dominant factors
# are first.
ii = sortperm(eg.values, rev=true)
pcw = eg.values[ii]
pcv = eg.vectors[:, ii]

# For interpretability flip the PC's that are
# mostly negative.
for j in 1:size(pcv, 2)
    if sum(pcv[:, j] .< 0) > sum(pcv[:, j] .>= 0)
        pcv[:, j] .*= -1
    end
end

# Get the PC scores for temperature
scores = tempc' * pcv[:, 1:20]

# Plot the j^th PC score against the k^th feature.
function plot_pcscores(j, k)

    # Plot the mean profile
    plt = lineplot(pressure, tempmean; xlabel="Pressure", ylabel="Mean temperature")
    println(plt)

    # Plot the PC loadings
    fn = ["Latitude", "Longitude", "Day"]
    plt = lineplot(pressure, pcv[:, j]; xlabel="Pressure", ylabel="PC $(j) loading")
    println(plt)

    # Plot the conditional mean PC score against an observed variable
    xx = range(extrema(Y[:, k])..., length=100)
    m = loess(Y[:, k], scores[:, j])
    yh = predict(m, Y[:, k])
    resid = scores[:, j] - yh
    r = loess(Y[:, k], abs.(resid), degree=1)
    yy = predict(m, xx)
    yr = predict(r, xx)

    f = 2
    ymx = maximum(yy + f*yr)
    ymn = minimum(yy - f*yr)
    plt = lineplot(xx, yy; ylabel="PC $(j) score", xlabel=fn[k], color=:red, ylim=[ymn, ymx])
    plt = lineplot!(plt, xx, yy-f*yr, color=:white)
    plt = lineplot!(plt, xx, yy+f*yr, color=:white)
    println(plt)
end

for k in 1:3
    for j in 1:3
        plot_pcscores(j, k)
    end
end
