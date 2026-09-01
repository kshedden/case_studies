import os
import pandas as pd

# Download the NCDC data, also available from this site:
#  https://www.ncei.noaa.gov/access/search/data-search/daily-summaries

# This file provides some information about each station, and allows you
# to match the file names to the station descriptions.
#  https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt

# Change this to point to a valid path on your computer.
target_dir = "/home/kshedden/data/Teaching/precip"

# Each file corresponds to a specific geographic location.
files = ["USW00094847.csv", "USW00012839.csv",]

base_url = "https://www.ncei.noaa.gov/data/daily-summaries/access"

for f in files:
    df = pd.read_csv(os.path.join(base_url, f))
    df.to_csv(os.path.join(target_dir, f + ".gz"))
