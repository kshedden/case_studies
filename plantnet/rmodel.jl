#=
Use mixed linear models to understand how the ocurrences of different
plant species as reported in the plantnet database are changing in
time.  Each mixed model uses one location variable (latitude,
longitude, elevation) as the response and controls for the other two.
Control variables are modeled with cubic polynomials to capture
nonlinearity.  The models include fixed and random effects for year to
capture the time trend.  The random slopes are unique to each species.
In addition to the model output we generate a plot of the predicted
linear trend for each species based on the mixed model BLUPS (best
linear unbiased predictors).
=#

using PyPlot, MixedModels, Printf, Statistics, StatsModels

rm("plots", force=true, recursive=true)
mkdir("plots")

include("read.jl")

# Each of these variables can be the outcome, and the other
# two are predictors.
va = [:decimalLatitude, :decimalLongitude, :elevation]

# Names for the variables in 'va' used for labeling plots.
titles = Dict(:decimalLatitude=>"Latitude", :decimalLongitude=>"Longitude", :elevation=>"Elevation")

# Build polynomial terms
for v in va
    z = df[:, v]
    z = (z .- mean(skipmissing(z))) ./ std(skipmissing(z))
    for j in 1:3
        w = Symbol(@sprintf("%s%d", string(v), j))
        df[:, w] = z.^j
    end
end

function fitmodel(v)

    # Control for the variables except for the outcome
    te = [x for x in va if x != v]
    tx = []
    for u in te # covariates
        for j in 1:3 # polynomial terms
            push!(tx, @sprintf("%s%d", string(u), j))
        end
    end
    tx = Symbol.(tx)

    # Ideally this would work but for some reason it does not
    fx = term(v) ~ term(:year) + sum(term.(tx)) + (term(1) + term(:decade) | term(:scientificName))

    # This is a workaround to the line above.
    s = join(string.(tx), " + ")
    c = "@formula($(v) ~ year + $(s) + (1 + decade | scientificName))"
    f0 = eval(Meta.parse(c))

    # Fit the mixed model and print the model summary.
    m0 = fit(MixedModel, f0, df)
    println(@sprintf("n=%d observations", nobs(m0)))
    rr = ranef(m0)
    println(@sprintf("n=%d species", size(rr[1], 2)))
    println(m0)
    println("\n\n")

    return m0, rr
end

# Generate a plot of predicted trends for each species.
function make_plots(mm, rr, va, ifig)

    # The fitted value at the central location, used as an intercept.
    xx = modelmatrix(mm)
    xm = mean(xx, dims=1)[:]
    xm[2] = 0
    icept = xm' * coef(mm)

    # The year slope (in year units.)
    ys = coef(mm)[2]

    rr = rr[1]

    # Plot these years
    yr = collect(range(2010, 2020, 20))

    # The plotted years as decades
    yx = (yr .- mean(yr)) / 10

    PyPlot.clf()
    PyPlot.grid(true)
    for i in 1:size(rr, 2)
        PyPlot.plot(yr, icept + rr[1, i] .+ ys.*yr + rr[2, i]*yx, "-", color="grey", alpha=0.4)
    end
    PyPlot.xlabel("Year", size=15)
    PyPlot.ylabel(titles[va], size=15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    return ifig + 1
end

function run_all(ifig)
    ifig = 0
    for v in va
        mm, rr = fitmodel(v)
        ifig = make_plots(mm, rr, v, ifig)
    end
    return ifig
end

ifig = 0
ifig = run_all(ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dAutoRotatePages=/None -dNOPAUSE -dBATCH -dSAFER -sOutputFile=rmodel_jl.pdf $f`
run(c)
