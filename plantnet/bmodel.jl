#=
Assess the evidence that the mean latitude, longitude, or elevation
of occurrences of a species are changing in time.  This is done by
fitting regression models with least squares, in which one of latitude,
longitude, and elevation is the dependent variable, and the other two are 
independent variables.  Two models are fit -- one including year
and one not including year, and we assess the statistical significance
of the year effect.
=#

using DataFrames, CSV, PyPlot, Dates, GLM, Distributions, Printf
using LinearAlgebra

rm("plots", force=true, recursive=true)
mkdir("plots")

include("read.jl")

# Generate basis functions for latitude, longitude, elevation,
# and day within year (seasonality).  Seasonality and longitude
# are circular variables.
setbasis = function(df)

	# Basis functions for day (seasonality).
    per = 365
    for k in 1:4
        df[!, "sin_day_$(k)"] = sin.(2*pi*df[:, :dayOfYear]/per)
        df[!, "cos_day_$(k)"] = cos.(2*pi*df[:, :dayOfYear]/per)
        per /= 2
    end

    # Basis functions for longitude.
    per = 180
    for k in 1:4
        df[!, "sin_lon_$(k)"] = sin.(2*pi*df[:, :decimalLongitude]/per)
        df[!, "cos_lon_$(k)"] = cos.(2*pi*df[:, :decimalLongitude]/per)
        per /= 2
    end

    # Basis functions for latitude.
    x = (df[:, :decimalLatitude] .- 45) / 100
    for k in 1:3
        df[!, "lat$(k)"] = x.^k
    end

    # Basis functions for latitude.
    x = (df[:, :elevation] .- 100) ./ 1000
    for k in 1:3
        df[!, "elv$(k)"] = x.^k
    end

	# Basis functions for year.
	x = (df[:, :year] .- 2010) ./ 100
	for k in 1:3
		df[!, "year$(k)"] = x.^k
	end
    
    return df
end

# These are the terms that will be used to build a regression model.
dayterms = [term("sin_day_$(k)") for k in 1:4]
dayterms = vcat(dayterms, [term("cos_day_$(k)") for k in 1:4])
lonterms = [term("sin_lon_$(k)") for k in 1:4]
lonterms = vcat(lonterms, [term("cos_lon_$(k)") for k in 1:4])
latterms = [term("lat$(k)") for k in 1:3]
elvterms = [term("elv$(k)") for k in 1:3]

# Adjust for everything except the response
function get_covariate_terms(response)
	terms = []
	v = [:day, :decimalLongitude, :decimalLatitude, :elevation]
	for (j,te) in enumerate([dayterms, lonterms, latterms, elvterms])
		if response != v[j]
			push!(terms, te...)
		end
	end
	return terms
end

# Fit a linear model for the response variable (should be one of
# latitude, longitude, or elevation) predicted by other control
# variables, with and without year.  If year is included, it
# is modeled linearly.
function fit_linear(response)

	terms = get_covariate_terms(response)

	fml0 = term(response) ~ sum(terms)
	fml1 = term(response) ~ term(:year) + sum(terms)

	pva = []
	for (k, dv) in pairs(groupby(df, :scientificName))
    	dv = setbasis(dv)
    	m0 = lm(fml0, dv)
    	m1 = lm(fml1, dv)
    	lrt = 2*(loglikelihood(m1) - loglikelihood(m0))
    	dof = length(coef(m1)) - length(coef(m0))
    	pv = 1 - cdf(Chisq(dof), lrt)
    	push!(pva, [pv, coef(m1)[2]])
	end
	pva = hcat(pva...)
	return pva
end

# Fit a nonlinear model for the response variable (should be one of
# latitude, longitude, or elevation) predicted by other control
# variables, with and without year.  If year is included, it
# is modeled nonlinearly with a cubic polynomial.
function fit_nonlin(response)

	terms = get_covariate_terms(response)

	fml0 = term(response) ~ sum(terms)
	fml1 = term(response) ~ term(:year) + term(:year2) + term(:year3) + sum(terms)

	pvb = []
	for (k, dv) in pairs(groupby(df, :scientificName))
    	dv = setbasis(dv)
    	m0 = lm(fml0, dv)
    	m1 = lm(fml1, dv)
    	lrt = 2*(loglikelihood(m1) - loglikelihood(m0))
    	dof = length(coef(m1)) - length(coef(m0))
    	pv = 1 - cdf(Chisq(dof), lrt)
    	push!(pvb, [pv, 0])
	end
	pvb = hcat(pvb...)
	return pvb
end

function make_plots(pva, ifig, plot_slopes, title)

	# Quantile/quantile plot of p-values
	PyPlot.clf()
	PyPlot.grid(true)
	pvas = sort(pva[1, :])
	n = length(pvas)
	x = range(1/n, 1-1/n, n)
	PyPlot.plot([0, 1], [0, 1], color="black")
	PyPlot.plot(x, pvas, "o", mfc="none")
	PyPlot.xlabel("Expected p-value", size=15)
	PyPlot.ylabel("Observed p-value", size=15)
	PyPlot.title(title)
	PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
	ifig += 1

	if plot_slopes
		# Scatterplot year slopes against p-values.
		PyPlot.clf()
		PyPlot.grid(true)
		PyPlot.plot(pva[1, :], pva[2, :], "o", mfc="none")
		PyPlot.xlabel("p-value", size=15)
		PyPlot.ylabel("Year slope", size=15)
		PyPlot.title(title)
		PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
		ifig += 1
	end

	return ifig
end

function run_all()
	ifig = 0
	title = Dict(:decimalLatitude=>"Latitude", :decimalLongitude=>"Longitude", :elevation=>"Elevation")
	for v in [:decimalLatitude, :decimalLongitude, :elevation]
		pva = fit_linear(v)
		pvb = fit_nonlin(v)
		ifig = make_plots(pva, ifig, true, title[v])
		ifig = make_plots(pvb, ifig, false, title[v])
	end
	return ifig
end

ifig = run_all()

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dAutoRotatePages=/None -dNOPAUSE -dBATCH -dSAFER -sOutputFile=bmodel_jl.pdf $f`
run(c)
