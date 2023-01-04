using CSV, DataFrames, Dates, Statistics

# The raw file from plantnet should be located here.
pa = "/home/kshedden/myscratch/plantnet"

# This is the raw data file.  If your file name does not
# match this needs to be changed.
fn = "0140072-220831081235567.csv.gz"

fn = joinpath(pa, fn)
df = open(fn) do io
    CSV.read(io, DataFrame)
end

# Keep the 200 most common species.
ds = combine(groupby(df, :scientificName), nrow)
sort!(ds, :nrow, rev=true)
species = ds[1:200, :scientificName]
species = Set(species)
df = filter(r->r.scientificName in species, df)
df[:, :Date] = Date.(df[:, :year], df[:, :month], df[:, :day])
df = df[:, [:scientificName, :Date, :elevation, :decimalLatitude, :decimalLongitude]]

# Calculate the circular mean
function circmean(x)
    x = pi * x / 180
    s = mean(sin, x)
    c = mean(cos, x)
    return 180 * atan(s, c) / pi
end

# Mean latitude, longitude, and elevation for each species.
dz = combine(groupby(df, :scientificName), :decimalLatitude=>mean=>:decimalLatitude,
                  :elevation=>(x->mean(skipmissing(x)))=>:elevation,
                  :decimalLongitude=>circmean=>:decimalLongitude)

dfa = combine(groupby(df, [:Date, :scientificName]), nrow=>:nobs)
dr = minimum(df[:, :Date]):Day(1):maximum(df[:, :Date])
db = DataFrame(Iterators.product(species, dr))
db = rename(db, 1=>:scientificName, 2=>:Date)
dc = leftjoin(db, dfa, on=[:Date, :scientificName])
dc[:, :nobs] = replace(dc[:, :nobs], missing=>0)
dc[:, :year] = year.(dc[:, :Date])
dc[:, :month] = month.(dc[:, :Date])
dc[:, :dayOfYear] = dayofyear.(dc[:, :Date])

df = sort(dc, [:scientificName, :Date])
CSV.write(joinpath(pa, "plants_occurrences.csv.gz"), df, compress=true)

dz = sort(dz, :scientificName)
CSV.write(joinpath(pa, "plants_locations.csv.gz"), dz, compress=true)
