"""
Examine the lifespans of notable people using the BHHT data.

This analysis uses survival analysis methods, allowing us
to use information from still-living people.
"""

using HazardRegression
using SurvivalPseudo
using DataFrames
using CSV
using UnicodePlots
using Printf

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/data/Teaching/bhht"

# Load the dataset.  Use the latin-1 encoding since there is some non-UTF
# data in the file.  Add "nrows=100000" when developing to reduce the run
# time (but use the complete data to get final results).
df = CSV.read(joinpath(pa, "cross-verified-database.csv.gz"), DataFrame)

# Create a lifespan variable (years of life).  It will be missing for people who are currently living.
df[!, :lifespan] = df[:, :death] - df[:, :birth]

# Exclude people born before 1500, there is too little data to gain a meaningful
# understanding of the trends in lifespan prior to this year.
dx = filter(r->!ismissing(r.birth) && r.birth >= 1500, df)
censor_year = maximum(skipmissing(dx[:, :death]))
dx = select(dx, [:birth, :lifespan, :gender, :un_region])

# There are a small number of people with missing or "Other" gender but it
# is too small of a sample to draw conclusions.
dx = filter(r->!ismissing(r.gender) && r.gender in ["Female", "Male"], dx)

# Censor
dx[!, :clifespan] = [ismissing(r.lifespan) ? censor_year - r.birth : r.lifespan for r in eachrow(dx)]
dx[!, :died] = [ismissing(x) ? 0 : 1 for x in dx.lifespan]

# Now we can drop all rows with missing data
dx = select(dx, Not(:lifespan))
dx = dx[completecases(dx), :]
dx = disallowmissing(dx)
dx = filter(r->r.clifespan > 0, dx)

# A categorical variable indicating the century in which a person was born.
dx[!, :era] = floor.((dx[:, :birth] .- 1500) / 100)

# Plot the survival functions for people born in each century
marginal_survival = function()
    fst = true
    plt = nothing
    for g in groupby(dx, :era)
        sf = SurvivalFunction(zeros(length(g.clifespan)), g.clifespan, g.died)
        era = first(g[:, :era])
        lab = @sprintf("%d", 1500+100*era)
        if fst
            plt = lineplot(sf.utime, sf.surv, name=lab, height=20, width=60,
                           xlim=(0, 100), xlabel="Age", ylabel="Proportion alive")
            fst = false
        else
            lineplot!(plt, sf.utime, sf.surv, name=lab)
        end
    end
    return plt
end

plt = marginal_survival()
println(plt)
