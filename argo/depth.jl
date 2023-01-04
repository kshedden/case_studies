using LinearAlgebra, UnicodePlots, Statistics

include("read.jl")

"""
Calculate the spatial depth of column i of x relative
to the other columns of x.  See A.2 of this reference
for a definition:

https://projecteuclid.org/journals/statistical-science/volume-32/issue-4/On-a-General-Definition-of-Depth-for-Functional-Data/10.1214/17-STS625.pdf
"""
function depth(i, x)

    p, n = size(x)
    d = zeros(p)
    u = zeros(p)

    for j in 1:size(x, 2)
        if j != i
            d .= x[:, i] - x[:, j]
            w = norm(d)
            if w > 0
                u .+= d / w
            end
        end
    end

    return 1 - norm(u / n)
end

# PCA of the temperature data
tempc = copy(temp)
for i in 1:size(tempc, 1)
    tempc[i, :] .-= mean(tempc[i, :])
end
cc = cov(tempc')
a,b = eigen(Symmetric(cc))
scores = tempc' * b[:, end-1:end]

# Get all the depths of the temperature data
dp = [depth(i, temp) for i in 1:size(temp, 2)]
ii = sortperm(dp)

qq = quantile(dp, 0.9)
i0 = dp .< qq
i1 = dp .>= qq
plt = scatterplot(scores[i0, 1], scores[i0, 2], color=:white)
scatterplot!(plt, scores[i1, 1], scores[i1, 2], color=:red)
println(plt)
