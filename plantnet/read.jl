using DataFrames, CSV, Dates

# The raw data file constructed by prep.jl should be available at this path.
pa = "/home/kshedden/myscratch/plantnet"

# Load the raw data
df = open(joinpath(pa, "short.csv.gz")) do io
    CSV.read(io, DataFrame)
end

# The data are heavily skewed toward the more recent years, optionally
# restrict the analysis to these years.
firstyear = 2010
df = filter(r->r.Date >= Date(firstyear, 1, 1), df)

# Generate some time variables
df[!, :year] = year.(df[:, :Date])
df[!, :dayOfYear] = dayofyear.(df[:, :Date])

# There are very few records from southern hemisphere
df = filter(r->r.decimalLatitude >= 0, df)

# Elevation is very skewed, to prevent the extreme values
# from dominating the fit, clamp at 3000 meters.
df[!, :elevation] = clamp.(df[:, :elevation], -Inf, 3000)

# Centered and scaled year (one unit of decode corresponds to 10
# years of time).
meanyear = mean(df[:, :year])
df[:, :decade] = (df[:, :year] .- meanyear) / 10
