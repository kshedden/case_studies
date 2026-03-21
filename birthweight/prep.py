"""
Extract birth record data from NCHS archives.  The raw data are in
fixed width format.  This script converts these to CSV files,
including birth weight and several predictors of birth weight that
seem to be consistently coded across years.

The data files and documentation are the "Birth data files" from here:

https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm

The same data files seem to also be available here:

https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/

The 1983 RUCC codes can be obtained here:

https://www.ers.usda.gov/data-products/rural-urban-continuum-codes
"""

import pandas as pd
from pathlib import Path
import requests
import shutil
import os
import re
import gzip
import pandas as pd
import numpy as np

# Directory in which to place the data files
spath = Path("/home/kshedden/data/Teaching/birthweight/births")

years = [1971, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1999, 1991, 1992]

# Column specification for 1969-1971.  The first number of each
# 2-tuple is the starting position (1-based counting), and the
# second number is the width of the field.
colspec_1969_1971={
    "year": (1,1),
    "state": (13,2),
    "county": (15,3),
    "smsa": (22, 3),
    "sex": (35,1),
    "dadrace": (37,1),
    "momrace": (38,1),
    "momage": (41,2),
    "birthorder": (58,2),
    "dadage": (69,2),
    "birthweight": (73,4),
    "plurality": (81,1),
    "interval": (117, 3),
    "popsize": (21, 1)}

# The column specification for 1981 is the same as 1969-1971 for the fields we are using
colspec_1981_1988={
    "year": (1,1),
    "state": (13,2),
    "county": (15,3),
    "smsa": (22,3),
    "sex": (35,1),
    "dadrace": (37,1),
    "momrace": (38,1),
    "momage": (41,2),
    "birthorder": (58,2),
    "dadage": (69,2),
    "birthweight": (73,4),
    "plurality": (81,1),
    "interval": (117, 3),
    "popsize": (21, 1)}

# Column specification for 1991 is different from above
colspec_1989_1992={
    "year": (1,4),
    "state": (32,2),
    "county": (34,3),
    "msa": (54,4),
    "sex": (189,1),
    "dadrace": (160,2),
    "momrace": (80,2),
    "momage": (70,2),
    "birthorder": (103,2),
    "dadage": (154,2),
    "birthweight": (193,4),
    "plurality": (201,1),
    "interval": (128, 3),
    "popsize": (26, 1)}

# Convert a column specification for use by Pandas.
def make_colspec(colspec):
    cs, cn = [], []
    for (k,v) in colspec.items():
        cn.append(k)
        cs.append([v[0]-1, v[0]+v[1]-1])
    return cs, cn

# Load the appropriate column specification for a given year.
# NOTE: Only tested for 1971, 1981, 1991, 192
def get_colspec(year):
    if 1969 <= year <= 1971:
        return colspec_1969_1971
    elif 1981 <= year <= 1988:
        return colspec_1981_1988
    elif 1989 <= year <= 1992:
        return colspec_1989_1992
    else:
        1/0

# Recode numerical race codes to strings, collapsing to four categories.
# This also works for 1981-1989
def recode_race_1971_1988(da):
    return da.replace({0: "Asian", 1: "White", 2: "Black", 3: "Native", 4: "Asian",
                       5: "Asian", 6: "Asian", 7: np.nan, 8: "Asian", 9: np.nan})

# Recode numerical race codes to strings, collapsing to four categories.
def recode_race_1989_1992(da):
    return da.replace({1: "White", 2: "Black", 3: "Native", 4: "Asian", 5: "Asian",
                       6: "Asian", 7: "Asian", 8: "Asian", 9: np.nan, 18: "Asian",
                       28: "Asian", 38: "Asian", 48: "Asian", 58: "Asian", 68: "Asian",
                       78: "Asian", 99: np.nan})

# Modify some of the fields to make the values consistent across years.
def clean_births(da, year):

    da["sex"] = da["sex"].replace({1: "male", 2: "female"})
    da["birthorder"] = da["birthorder"].replace({99: np.nan})
    da["dadage"] = da["dadage"].replace({99: np.nan})
    da["interval"] = da["interval"].replace({888: np.nan, 999: np.nan, 777: np.nan})
    da["birthweight"] = da["birthweight"].replace({9999: np.nan})

    if 1969 <= year <= 1971:
        da["year"] = da["year"].replace({9: 1969, 0: 1970, 1: 1971})
        da["momrace"] = recode_race_1971_1988(da["momrace"])
        da["dadrace"] = recode_race_1971_1988(da["dadrace"])
    elif 1981 <= year <= 1988:
        da["year"] += 1980
        da["momrace"] = recode_race_1971_1988(da["momrace"])
        da["dadrace"] = recode_race_1971_1988(da["dadrace"])
    elif 1989 <= year <= 1992:
        da["momrace"] = recode_race_1989_1992(da["momrace"])
        da["dadrace"] = recode_race_1989_1992(da["dadrace"])

    return da

# Convert the fixed width files to CSV, select and process variables for consistency across years.
def get_births(year, spath):
    colspec = get_colspec(year)
    cs, cn = make_colspec(colspec)
    da = pd.read_fwf(spath / f"Natl{year}.pub.gz", colspecs=cs, compression="gzip", header=None)
    da.columns = cn
    da = clean_births(da, year)
    return da

def download(year, spath):

    # Download the zip archive
    url_pat = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/"
    url = f"{url_pat}/Nat{year}.zip"
    response = requests.get(url)
    target = spath / f"Nat{year}.zip"
    open(target, "wb").write(response.content)

    # Unpack the zip archive and remove the archive file
    shutil.unpack_archive(target, spath)
    target.unlink()

    # Gzip the data file
    # NOTE: casing of the raw fwf filenames is irregular, they all
    # seem to end in 'pub'
    fnames = os.listdir(spath)
    r = re.compile(str(year) + "(.pub|.txt|)$", re.IGNORECASE)
    fnames = [f for f in fnames if r.search(f)]
    assert len(fnames) == 1
    fname = fnames[0]
    source = spath / fname
    target = spath / f"Natl{year}.pub.gz"
    with open(source, "rb") as fin:
        with gzip.open(target, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    source.unlink()

for year in years:
    download(year, spath)

for year in years:
    da = get_births(year, spath)
    da.to_csv(spath / f"{year}.csv.gz", index=False)
