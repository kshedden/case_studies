# Examine factors associated with birth count variation among US
# counties using Principal Components Regression and Poisson GLM/GEE.

using GLM, GEE, Statistics, UnicodePlots, LinearAlgebra, Printf

include("prep.jl")

da = leftjoin(births, pop, on = :FIPS)
da = leftjoin(da, rucc, on = :FIPS)
da[:, :logPop] = log.(da[:, :Population])

da = da[completecases(da), :]
da = disallowmissing(da)
da = sort(da, [:FIPS, :year])

# Calculate the mean and variance within each county to
# assess the mean/variance relationship.
mv = combine(groupby(births, :FIPS), :Births => mean, :Births => var)
scatterplot(log.(mv[:, :Births_var]), log.(mv[:, :Births_mean]))

# GLM, not appropriate since we have repeated measures on counties
fml = @formula(Births ~ logPop + RUCC_2013)
m0 = glm(fml, da, Poisson())

# GEE accounts for correlated data.
m1 = gee(fml, da, da[:, :FIPS], Poisson())

# GEE with log population as offset instead of covariate
fml = @formula(Births ~ RUCC_2013)
m2 = gee(fml, da, da[:, :FIPS], Poisson(); offset=da[:, :logPop])

# Use gamma family with non-canonical log link function to better match
# the mean/variance relationship
m3 = gee(fml, da, da[:, :FIPS], Gamma(), IndependenceCor(), LogLink(); offset=da[:, :logPop])

# Process the demographic data, -- replace missing values with 0
# and transform with square root to stabilize the variance.
for c in names(demog)
    if c != "FIPS"
        demog[:, c] = replace(demog[:, c], missing => 0)
        demog[:, c] = sqrt.(demog[:, c])
        demog[:, c] .-= mean(demog[:, c])
    end
end

# Get factors (principal components) from the demographic data
u, s, v = svd(Matrix{Float64}(demog[:, 2:end]))

# The proportion of explained variance.
pve = s .^ 2
pve ./= sum(pve)

# Put the PC's into a dataframe and merge into the dataframe for modeling.
demog_f = DataFrame(:FIPS => demog[:, :FIPS])
for k = 1:100
    demog_f[:, @sprintf("pc%02d", k)] = u[:, k]
end
da = leftjoin(da, demog_f, on = :FIPS)

# GLM, not appropriate since we have repeated measures on counties
npc = 5
fml =
    term(:Births) ~
        term(:logPop) + term(:RUCC_2013) + sum([term(@sprintf("pc%02d", j)) for j = 1:npc])
m4 = glm(fml, da, Poisson())
m5 = gee(fml, da, da[:, :FIPS], Poisson(), offset=da[:, :logPop])

# Include this number of factors in subsequent models
npc = 20

# Convert the coefficients back to the original coordinates
function convert_coef(c, npc)
    return v[:, 1:npc] * (c ./ s[1:npc])
end

# Restructure the coefficients so that the age bands are
# in the columns.
function restructure(c)
    cc = copy(na)
    cc[:, :coef] = c
    cc = unstack(cc, [:Race, :Origin, :Sex], :Age, :coef)
    return cc
end

# This function fits a Poisson GLM to the data using 'npc' principal components
# as explanatory variables.
function fitmodel(npc)
    # A GEE using log population as an offset
    fml = term(:Births) ~ sum([term(@sprintf("pc%02d", j)) for j = 1:npc])
    m = gee(fml, da, da[:, :FIPS], Poisson(), offset = da[:, :logPop])

    # We need this for score-testing
    m0 = gee(fml, da, da[:, :FIPS], Poisson(), offset = da[:, :logPop]; dofit=false)

    # Convert the coefficients back to the original coordinates
    c = convert_coef(coef(m)[2:end], npc)

    # Restructure the coefficients so that the age bands are
    # in the columns.
    c = restructure(c)

    return c, m, m0
end

# Fit models with these numbers of PCs.
pcs = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

models = []
for npc in pcs
    c, m, m0 = fitmodel(npc)
    push!(models, (m, m0))
end

for i in 2:length(models)
    smaller = models[i-1][1]
    bigger = models[i][2]
    st = scoretest(bigger.model, smaller.model)
    println(st)
end
