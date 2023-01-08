"""
Examine the lifespans of notable people using the BHHT data.

This analysis uses survival analysis methods, allowing us
to use information from still-living people.
"""

using SurvivalAnalysis
using DataFrames
using CSV
using UnicodePlots
using Printf

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"

# Load the dataset.  Use the latin-1 encoding since there is some non-UTF
# data in the file.  Add "nrows=100000" when developing to reduce the run
# time (but use the complete data to get final results).
df = open(joinpath(pa, "cross-verified-database.csv.gz")) do io
    CSV.read(io, DataFrame) #, limit=100000)
end

# Create a lifespan variable (years of life).  It will be missing for people who are currently living.
df[!, :lifespan] = df[:, :death] - df[:, :birth]

# Exclude people born before 1500, there is too little data to gain a meaningful
# understanding of the trends in lifespan prior to this year.
dx = filter(r->!ismissing(r.birth) && r.birth >= 1500, df)
dx = select(dx, [:birth, :lifespan, :gender, :un_region, :level1_main_occ])

# There are a small number of people with missing or "Other" gender but it
# is too small of a sample to draw conclusions.
dx = filter(r->!ismissing(r.gender) && r.gender in ["Female", "Male"], dx)

# Drop uninformative occupation codes.
dx = filter(r->!ismissing(r.level1_main_occ) && !(r.level1_main_occ in ["Missing", "Other"]), dx)

# Censor at 2022
censor_year = 2022
dx[!, :clifespan] = [ismissing(r.lifespan) ? censor_year - r.birth : r.lifespan for r in eachrow(dx)]
dx[!, :died] = [ismissing(x) ? 0 : 1 for x in dx.lifespan]

# Now we can drop all rows with missing data
dx = select(dx, Not(:lifespan))
dx = dx[completecases(dx), :]
dx = disallowmissing(dx)

# A categorical variable indicating the century in which a person was born.
dx[!, :era] = floor.((dx[:, :birth] .- 1500) / 100)

# Plot the survival functions for people born in each century
marginal_survival = function()
    fst = true
    plt = nothing
    for g in groupby(dx, :era)
        sf = kaplan_meier(@formula(Srv(clifespan, died) ~ 1), g)
        y = predict(sf, dx[1:1, :]).survival_matrix
        era = first(g[:, :era])
        lab = @sprintf("%d", 1500+100*era)
        if fst
            plt = lineplot(y.time, y.survival[:], name=lab, height=20, width=60,
                           xlim=(0, 100), xlabel="Age", ylabel="Proportion alive")
            fst = false
        else
            lineplot!(plt, y.time, y.survival[:], name=lab)
        end
    end
    return plt
end

plt = marginal_survival()
println(plt)
