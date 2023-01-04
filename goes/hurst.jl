using Statistics

include("read.jl")

nn = [4, 8, 16, 32, 64, 128]

ti = df[:, :Time]
fl = df[:, :Flux1]

# Estimate the Hurst parameter.
function hurst(ti, fl, dr, d)

	r = zeros(length(dr), 2)
	for (j, m) in enumerate(dr)

		# Generate a matrix of non-overlapping blocks of
		# size m.
		_, flx = make_blocks(ti, fl, m, d)

		# Calculate the sample mean of each block.
		bm = mean(flx, dims=1)[:]

		# Take the sample variance of the block means.
		r[j, :] = [m, var(bm)]
	end

	# Estimate the Hurst exponent from the variances of
	# the sample means.
	rl = log.(r)
	b = cov(rl[:, 1], rl[:, 2]) / var(rl[:, 1])
	a = mean(rl[:, 2]) - b*mean(rl[:, 1])

	return 1 + b/2
end

# As a check, estimate the Hurst parameter for 
# IID normal data (the true value is 1/2). 
fx = randn(length(ti))
h0 = hurst(ti, fx, nn, 0)
println(h0)

# As another check, simulate correlated data
# with short-range dependence (the true value
# is 1/2).
fx = randn(length(ti))
r = 0.5
for i in 2:length(ti)
	fx[i] = r*fx[i-1] + sqrt(1 - r^2)*fx[i]
end
h0 = hurst(ti, fx, nn, 0)
println(h0)

# Estimate the Hurst Parameter for the GOES data.
h0 = hurst(ti, 1e8*fl, nn, 0)
h1 = hurst(ti, 1e8*fl, nn, 1)
h2 = hurst(ti, 1e8*fl, nn, 2)
println([h0, h1, h2])
