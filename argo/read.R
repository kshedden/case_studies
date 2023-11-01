library(readr)

# There is no R version of prep.py, you must run the
# Python or julia prep script before using read.R.

# This path must match the path in prep.py.
qpath = "/home/kshedden/data/Teaching/argo/python"

Sys.setenv(VROOM_CONNECTION_SIZE=1000000)

lat = read_csv(file.path(qpath, "lat.csv.gz"))
lat = as.vector(lat[,1])$Column1

lon = read_csv(file.path(qpath, "lon.csv.gz"))
lon = as.vector(lon[,1])$Column1

date = read_csv(file.path(qpath, "date.csv.gz"))
date = as.vector(date[,1])$Column1
day = date - min(date)

pressure = read_csv(file.path(qpath, "pressure.csv.gz"))
pressure = as.vector(pressure[,1])$Column1

temp = read_csv(file.path(qpath, "temp.csv.gz"))
temp = as.matrix(temp)

psal = read_csv(file.path(qpath, "psal.csv.gz"))
psal = as.matrix(psal)
