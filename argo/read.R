library(readr)
library(dplyr)
library(lubridate)

# There is no R version of prep.py, you must run the
# Python or julia prep script before using read.R.

# This path must match the path in prep.py.
qpath = "/home/kshedden/data/Teaching/argo/python"

# The latitude where each profile was recorded
lat = read_csv(file.path(qpath, "lat.csv.gz"))
lat = as.vector(lat[,1])$Column1

# The longitude where each profile was recorded
lon = read_csv(file.path(qpath, "lon.csv.gz"))
lon = as.vector(lon[,1])$Column1

# The day of year on which a profile was recorded
date = read_csv(file.path(qpath, "date.csv.gz"))
date = as.vector(date[,1])$Column1
day = yday(date)

# The pressure (in decibars) for each point along a profile.
pressure = read_csv(file.path(qpath, "pressure.csv.gz"))
pressure = as.vector(pressure[,1])$Column1

# The temperature profiles
temp = read_csv(file.path(qpath, "temp.csv.gz"), col_names=F)
temp = as.matrix(temp) %>% t

# The salinity profiles
psal = read_csv(file.path(qpath, "psal.csv.gz"), col_names=F)
psal = as.matrix(psal) %>% t
