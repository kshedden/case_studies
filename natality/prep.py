# Prepare the birth count and other related data for analysis.

import pandas as pd
import numpy as np
import os, gzip

# Path to the data files
pa = "/home/kshedden/data/Teaching/natality"

# Create a long form version of the birth data.
dl = []
for y in range(2011, 2021):
    da = pd.read_csv(os.path.join(pa, "%4d.txt.gz" % y), delimiter="\t",
                     dtype={"County Code": object})
    da = da[["County", "County Code", "Births"]]
    da["year"] = y
    dl.append(da)
births = pd.concat(dl)
births = births.rename({"County Code": "FIPS"}, axis=1)
births["Births"] = pd.to_numeric(births.Births, errors="coerce")
births = births.dropna()
births = births[~births["County"].str.contains("Unidentified")]

# Subset the demographics file to 2016.  This is a large file
# and takes time, so only run this once and save the result
# for future runs.
if False:
    f = os.path.join(pa, "us.1990_2020.19ages.txt.gz")
    g = os.path.join(pa, "2016ages.txt.gz")
    with gzip.open(g, "w") as out:
        with gzip.open(f) as inp:
            for line in inp:
                if line.startswith(b"2016"):
                    out.write(line)

# Read the demographics for 2016.  It is a fixed-width format
# file.
x = [1, 5, 7, 9, 12, 14, 15, 16, 17, 19, 27]
cs = [(x[i]-1, x[i+1]-1) for i in range(len(x)-1)]
with gzip.open(os.path.join(pa, "2016ages.txt.gz")) as io:
    demog = pd.read_fwf(io, colspecs=cs, header=None)
demog.columns = ["Year", "State", "StateFIPS", "CountyFIPS", "Registry",
                 "Race", "Origin", "Sex", "Age", "Population"]

# Create a FIPS code that matches the FIPS code in the birth data
demog["FIPS"] = ["%02d%03d" % (x, y) for (x, y) in zip(demog.StateFIPS, demog.CountyFIPS)]
demog = demog[["FIPS", "Race", "Origin", "Sex", "Age", "Population"]]

# Recode some variables to more interpretable text labels
# See http://seer.cancer.gov/popdata/popdic.html for code information
demog["Sex"] = demog["Sex"].replace([2, 1], ["F", "M"]) # Female/Male
demog["Origin"] = demog["Origin"].replace([0, 1], ["N", "H"]) # Non-Hispanic/Hispanic
demog["Race"] = demog["Race"].replace([1, 2, 3, 4], ["W", "B", "N", "A"]) # White/Black/Native/Asian

# The overall population per county
pop = demog.groupby("FIPS")["Population"].sum().reset_index()

# Pivot to put the age bands in the columns
demog = demog.pivot(index="FIPS", columns=["Race", "Origin", "Sex", "Age"], values="Population")

# Age group labels
age_groups = ["0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
              "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79",
              "80-84", "85+"]

na = demog.columns.tolist()
demog.columns = ["%s_%s_%s_%d" % tuple(x) for x in na]

# Replace missing demographic values with 0 and transform
# with square root to stabilize the variance.
demog = demog.fillna(0)
demog = np.sqrt(demog)

# Get the Rural/Urban Continuity Codes (RUCC)
rucc = pd.read_excel(os.path.join(pa, "ruralurbancodes2013.xls"), sheet_name=None)
rucc = rucc["Rural-urban Continuum Code 2013"]
rucc = rucc[["FIPS", "RUCC_2013"]]
rucc["FIPS"] = ["%05d" % x for x in rucc.FIPS]
