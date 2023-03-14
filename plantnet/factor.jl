using CSV, DataFrames, Statistics, LinearAlgebra
using PyPlot, Printf

mdates = PyPlot.matplotlib.dates

ifig = 0

rm("plots", recursive=true, force=true)
mkdir("plots")

pa = "/home/kshedden/mynfs/data/Teaching/plantnet"

df = open(joinpath(pa, "plants_occurrences.csv.gz")) do io
    CSV.read(io, DataFrame)
end

dz = open(joinpath(pa, "plants_locations.csv.gz")) do io
    CSV.read(io, DataFrame)
end

# Pivot the data so that the species are in the columns
# and the dates are in the rows.
dx = df[:, [:Date, :scientificName, :nobs]]
dx = filter(r->r.Date >= "2018-01-01", dx)
dx = unstack(dx, :Date, :scientificName, :nobs)
dates = dx[:, :Date]
dx = select(dx, Not(:Date))
species = names(dx)
dx = dx[:, sortperm(species)]
species = names(dx)
dx = Matrix(dx)

# Make sure that the occurrence data and meta-data are ordered
# the same way.
@assert all(dz[:, :scientificName] .== species)

# Variance stabilizing transformation
dx = sqrt.(dx)

speciesmeans = mean(dx, dims=1)[:]
datemeans = mean(dx, dims=2)[:]

# Double center the data
dx .-= mean(dx)
for j in 1:size(dx, 2)
    dx[:, j] .-= mean(dx[:, j])
end
for i in 1:size(dx, 1)
    dx[i, :] .-= mean(dx[i, :])
end

function plotmeans(ifig)

    # Plot the mean at each date (this is the mean number of
    # occurrences summed over all species on a specific date).
    PyPlot.clf()
    PyPlot.grid(true)
    PyPlot.plot(dates, datemeans, "-")
    PyPlot.gca().xaxis.set_major_locator(mdates.YearLocator(5))
    for x in PyPlot.gca().xaxis.get_ticklabels()
        x.set_rotation(-90)
    end
    PyPlot.xlabel("Date", size=15)
    PyPlot.ylabel("Date mean (mean occurrences/species)", size=15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    # The species have no order so plot the order statistics
    # of the species means.
    PyPlot.clf()
    PyPlot.grid(true)
    PyPlot.plot(sort(speciesmeans), "-")
    PyPlot.xlabel("Sorted position", size=15)
    PyPlot.ylabel("Species mean (mean occurrences/day)", size=15)
    PyPlot.title("Species mean order statistics")
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    return ifig
end

ifig = plotmeans(ifig)

# Factor the matrix once, then sort so that the
# species scores are increasing for the first
# factor.
u, s, v = svd(dx)
ii = sortperm(v[:, 1])
dx = dx[:, ii]
dz = dz[ii, :]

# Factor the matrix again.
u, s, v = svd(dx)

function plotscree(ifig)

    # A scree plot (singular value versus position).
    PyPlot.clf()
    PyPlot.grid(true)
    PyPlot.plot(s, "-")
    PyPlot.xlabel("SVD component", size=15)
    PyPlot.ylabel("Singular value", size=15)
    PyPlot.title("Scree plot")
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    # The scree plot in log/log space.
    s1 = s[s.>1e-8]
    PyPlot.clf()
    PyPlot.grid(true)
    PyPlot.plot(log.(1:length(s1)), log.(s1), "-")
    PyPlot.xlabel("Log SVD component", size=15)
    PyPlot.ylabel("Log singular value", size=15)
    PyPlot.title("Scree plot (log space)")
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    return ifig
end

ifig = plotscree(ifig)

function make_plots(ifig)
    for j in 1:10

        # Plot the left-side (date) loadings for component j.
        PyPlot.clf()
        PyPlot.axes([0.13, 0.2, 0.8, 0.7])
            PyPlot.grid(true)
        PyPlot.plot(dates, u[:, j], "-")
        PyPlot.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
        for x in PyPlot.gca().xaxis.get_ticklabels()
            x.set_rotation(-90)
        end
        PyPlot.xlabel("Date", size=15)
        PyPlot.ylabel("Date factor $(j) score", size=15)
        PyPlot.title("Factor $(j)")
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1

        # Plot the right side (date) loadings for component j.
        PyPlot.clf()
        PyPlot.axes([0.15, 0.2, 0.8, 0.7])
        PyPlot.grid(true)
        PyPlot.title("Factor $(j)")
        PyPlot.ylabel("Species factor", size=15)
        PyPlot.xlabel("Species", size=15)
        PyPlot.plot(v[:, j], "-")
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1

        # Plot the species scores against the mean latitude at which
        # each species is observed.
        PyPlot.clf()
        PyPlot.axes([0.15, 0.2, 0.8, 0.7])
        PyPlot.grid(true)
        PyPlot.plot(dz[:, :decimalLatitude], v[:, j], "o", mfc="none")
        PyPlot.xlabel("Mean latitude", size=15)
        PyPlot.ylabel("Factor $(j) score", size=15)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1

        # Plot the species scores against the mean longitude at which
        # each species is observed.
        PyPlot.clf()
        PyPlot.axes([0.15, 0.2, 0.8, 0.7])
        PyPlot.grid(true)
        PyPlot.plot(dz[:, :decimalLongitude], v[:, j], "o", mfc="none")
        PyPlot.xlabel("Mean longitude", size=15)
        PyPlot.ylabel("Factor $(j) score", size=15)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1

        # Plot the species scores against the mean elevation at which
        # each species is observed.
        PyPlot.clf()
        PyPlot.grid(true)
        PyPlot.plot(dz[:, :elevation], v[:, j], "o", mfc="none")
        PyPlot.xlabel("Mean elevation", size=15)
        PyPlot.ylabel("Factor $(j) score", size=15)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1
    end

    return ifig
end

ifig = make_plots(ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dAutoRotatePages=/None -dNOPAUSE -dBATCH -dSAFER -sOutputFile=factor_jl_plots.pdf $f`
run(c)
