#=
Examine the lifespans of notable people using the BHHT data.

The main tool here is a form of local regression known as
"LOESS", which is a form of local polynomial regression.
=#

using DataFrames, CSV, DataFrames, UnicodePlots, Loess

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"

# Load the dataset.
df = open(joinpath(pa, "cross-verified-database.csv.gz")) do io
    CSV.read(io, DataFrame)
end

# Create a lifespan variable (years of life).
df[:, :lifespan] = df[:, :death] - df[:, :birth]

# The loess routines used below require the data to be Float-type
# and completely non-missing.
dx = df[:, [:birth, :lifespan, :gender, :level1_main_occ]]
dx = dx[completecases(dx), :]
dx[!, :birth] = Float64.(dx[:, :birth])
dx[!, :lifespan] = Float64.(dx[:, :lifespan])

# To avoid censoring, exclude people born after 1920.
dx = filter(r -> r.birth >= 1500 && r.birth <= 1920, dx)

# Group the data by the values of the nominal variable 'vname',
# then use loess to estimate the conditional mean lifespan
# given age for each group.
function analyze(dx, vname)

    plt = nothing

    # Drop cases from groups that are too small to yield informative
    # results.
    dx = if vname == :gender
        filter(r->r[vname] in ["Female", "Male"], dx)
        else
        filter(r->!(r[vname] in ["Other", "Missing"]), dx)
    end

    for dd in groupby(dx, vname)

        vn = first(dd[:, vname])
        println(vn, " ", size(dd, 1))

        # Fit the loess model.
        dd = sort(dd, :birth)
        ll = loess(dd.birth, dd.lifespan; degree=1)

        # Use the fitted loess to predict values for plotting.
        xx = range(extrema(dd.birth)..., 100)
        yy = predict(ll, xx)

        if isnothing(plt)
            plt = lineplot(xx, yy, name = vn, ylim=[50, 80])
        else
            lineplot!(plt, xx, yy, name = vn)
        end
    end

    println(plt)
end

analyze(dx, :gender)
analyze(dx, :level1_main_occ)
