using DataFrames
using CSV
using DataFrames
using UnicodePlots
using MultivariateStats
using Printf
using StatsBase
using UnicodePlots
using PyPlot

function variable_plot(mca::MCA; x = 1, y = 2, ordered=[], kwargs...)

    fig = PyPlot.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.grid(true)

    # Set up the colormap
    cm = get(kwargs, :cmap, PyPlot.get_cmap("tab10"))

    # Set up the axis limits
    G = variable_coords(mca)
    G = DataFrame(Variable=G.Variable, Level=G.Level, Coord1=G.Coord[:, x], Coord2=G.Coord[:, y])

    mn1 = 1.2*minimum(G.Coord1)
    mx1 = 1.2*maximum(G.Coord1)
    mn2 = 1.2*minimum(G.Coord2)
    mx2 = 1.2*maximum(G.Coord2)

    xlim = get(kwargs, :xlim, [mn1, mx1])
    ylim = get(kwargs, :ylim, [mn2, mx2])
    ax.set_xlim(xlim...)
    ax.set_ylim(ylim...)

    for (j,g) in enumerate(groupby(G, :Variable))

        vname = first(g[:, :Variable])
        if vname in ordered
            PyPlot.plot(g[:, :Coord1], g[:, :Coord2], "-", color=cm(j))
        end

        for (k, v) in enumerate(eachrow(g))
            lb = "$(vname)-$(v.Level)"
            ax.text(g[k, :Coord1], g[k, :Coord2], lb, color=cm(j), ha = "center", va = "center")
        end
    end

    PyPlot.xlabel("Dimension $(x)")
    PyPlot.ylabel("Dimension $(y)")

    return fig
end

# Load the dataset.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"
df = open(joinpath(pa, "cross-verified-database.csv.gz")) do io
    CSV.read(io, DataFrame)
end

dx = df[:, [:birth, :un_region, :gender, :level1_main_occ]]
dx = rename(dx, :un_region=>:reg)
dx = rename(dx, :gender=>:sex)
dx = rename(dx, :level1_main_occ=>:occ)

dx = dx[completecases(dx), :]
dx = filter(r->r.birth >= 1500, dx)
dx = filter(r->r.occ != "Other", dx)
dx = filter(r->r.occ != "Missing", dx)
dx = filter(r->r.sex != "Other", dx)

dx[!, :birth] = Float64.(dx[:, :birth])

# Era of birth, round birth year to the nearest 50 years
dx[:, :era] = round.(2*dx[:, :birth]; digits=-2)./2
dx = select(dx, Not(:birth))
dx[!, :era] = [@sprintf("%4d", x) for x in dx[:, :era]]

dx = disallowmissing(dx)
for x in names(dx)
    dx[!, x] = String.(dx[:, x])
end

m = fit(MCA, dx; d=3)

# Demonstrate that people with similar object scores tend to have
# similar values on the analysis variables.
F = object_coords(m.C)
ii = sortperm(F[:, 1])
for _ in 1:10
    jj = sample(1:length(ii)-1)
    println("Similar pair:")
    println(dx[[ii[jj], ii[jj+1]], :])
    println("")
end

# Demonstrate that people with different object scores tend to have
# dissimilar values on the analysis variables.
for _ in 1:10
    jj = sample(1:1000)
    jj = [jj, length(ii) - jj]
    println("Dissimilar pair:")
    println(dx[ii[jj], :])
    println("")
end

# Plot the category scores
p = variable_plot(m; x=1, y=2, ordered=["era"])
p.savefig("bhht_mca_12.pdf")
p = variable_plot(m; x=1, y=3, ordered=["era"])
p.savefig("bhht_mca_13.pdf")
p = variable_plot(m; x=2, y=3, ordered=["era"])
p.savefig("bhht_mca_23.pdf")

# Add an articial variable that is independent of all other
# variables, to demonstrate that MCA places such variables
# at the origin.
dx[:, :junk] = sample(["A", "B", "C"], size(dx, 1))
m2 = fit(MCA, dx; d=3)
p = variable_plot(m2; x=1, y=2, ordered=["era"])
p.savefig("bhht_mca_12_junk.pdf")
dx = select(dx, Not(:junk))

# Try to show how well we are separating the objects.
F = object_coords(m2)
for k in 1:size(dx, 2)

    vn = m2.vnames[k]
    z = []
    u = unique(dx[:, k])
    sort!(u)
    for v in u
        ii = findall(dx[:, k] .== v)
        push!(z, F[ii, 1])
    end

    PyPlot.clf()
    PyPlot.title(vn)
    PyPlot.boxplot(z, labels=u)
    PyPlot.gca().set_ylabel("Component 1 object score")
    PyPlot.savefig("obj_scores_$(vn).pdf")
end
