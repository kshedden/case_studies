import os
import pandas as pd

# Download the NCDC data, also available from this site:
#  https://www.ncei.noaa.gov/access/search/data-search/daily-summaries

base_url = "https://www.ncei.noaa.gov/data/daily-summaries/access"

target_dir = "/home/kshedden/data/Teaching/precip"

files = ["USW00094847.csv", "USW00012839.csv",]

for f in files:
    df = pd.read_csv(os.path.join(base_url, f))
    df.to_csv(os.path.join(target_dir, f + ".gz"))
