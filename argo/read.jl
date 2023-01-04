using CSV, DataFrames

dpath = "/home/kshedden/data/Teaching/argo/julia"

lat = open(joinpath(dpath, "lat.csv.gz")) do io
    CSV.read(io, DataFrame)
end
lat = lat[:, 1]

lon = open(joinpath(dpath, "lon.csv.gz")) do io
    CSV.read(io, DataFrame)
end
lon = lon[:, 1]

date = open(joinpath(dpath, "date.csv.gz")) do io
    CSV.read(io, DataFrame)
end
date = date[:, 1]
date = [x.value for x in date - minimum(date)]

pressure = open(joinpath(dpath, "pressure.csv.gz")) do io
    CSV.read(io, DataFrame)
end
pressure = pressure[:, 1]

temp = open(joinpath(dpath, "temp.csv.gz")) do io
    CSV.read(io, DataFrame)
end
temp = Matrix(temp)

psal = open(joinpath(dpath, "psal.csv.gz")) do io
    CSV.read(io, DataFrame)
end
psal = Matrix(psal)

