# Run this script first to download the data.  You will need to configure
# some of the variables below to match the directory layout on your machine.

# This file provides some information about each station, and allows you
# to match the file names to the station descriptions.
#  https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt

library(R.utils)

base_url = "https://www.ncei.noaa.gov/data/daily-summaries/access"

# Configure this to point to a valid directory on your computer.
target_dir = "/home/kshedden/data/Teaching/precip"

# Download files for these locations.
files = c("USW00094847.csv", "USW00012839.csv")

for (f in files) {
    url = file.path(base_url, f)
    target = file.path(target_dir, f)
    download.file(url, target)
    gzip(target, overwrite=TRUE)
}
