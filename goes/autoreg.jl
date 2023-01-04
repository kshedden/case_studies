#=
Use ridge regression to predict X-ray flux at a given time
distance into the future, using a block of consecutive values.
=#

using Statistics, LinearAlgebra, PyPlot, Printf

rm("plots", recursive=true, force=true)
mkdir("plots")

include("read.jl")

ti = df[:, :Time]
fl = df[:, :Flux1]

# Use blocks of size m, and use the first q observations to predict
# the final observation.
m = 1000
q = 200

# The time points of the predictor information relative to the time
# being predicted.
tax = range(-2*m, -2*(m-q))[1:q]
tax = tax / 60

_, flx = make_blocks(ti, fl, m, 0)

flx = log.(1e-8 .+ flx)

# Construct a design matrix and response vector.
x = flx[1:q, :]'
y = flx[end, :]

# Center the data
y .-= mean(y)
for j in 1:size(x, 2)
	x[:, j] .-= mean(x[:, j])
end

# Regress y on x using ridge regression, with penalty parameter f.
function ridge(x, y, f)
	u, s, v = svd(x)
	b = v * Diagonal(s ./ (s.^2 .+ f)) * u' * y
	return b
end

# Consider how the regression coefficients look for various values
# of the penalty parameter.
function doridge(ifig)
	for f in Float64[1, 10, 100, 1000, 10000]
		b = ridge(x, y, f)
		PyPlot.clf()
		PyPlot.grid(true)
		PyPlot.plot(tax, b, "-")
		PyPlot.title(@sprintf("f=%.0f", f))
		PyPlot.ylabel("Coefficient", size=15)
		PyPlot.xlabel("Minutes before current time", size=15)
		PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
		ifig += 1
	end

	return ifig
end

ifig = 0
ifig = doridge(ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dAutoRotatePages=/None -dNOPAUSE -dBATCH -dSAFER -sOutputFile=autoreg_jl.pdf $f`
run(c)
