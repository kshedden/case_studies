import os
import numpy as np
import pandas as pd

# The raw file from plantnet should be located here.
pa = "/home/kshedden/myscratch/plantnet"

# This is the raw data file.  If your file name does not
# match it this needs to be changed.
fn = "0140072-220831081235567.csv.gz"

df = pd.read_csv(os.path.join(pa, fn), delimiter="\t")

# Keep the 200 most common species.
ds = df.groupby("scientificName").size()
ds = ds.sort_values(ascending=False)
species = ds.index[0:200].values
df = df.loc[df.scientificName.isin(set(species)), :]
df["Date"] = pd.to_datetime(df[["year", "month", "day"]])

df = df[["scientificName", "Date", "elevation", "decimalLatitude", "decimalLongitude"]]
df.to_csv(os.path.join(pa, "short.csv.gz"), index=None)

# Calculate the circular mean.
def circmean(x):
    x = np.pi * x / 180
    s = np.sin(x).mean()
    c = np.cos(x).mean()
    return 180 * np.arctan2(s, c) / np.pi

# Mean latitude, longitude, and elevation for each species.
dz = df.groupby("scientificName").agg({"decimalLatitude": np.mean, "elevation": np.mean, 
                                       "decimalLongitude": circmean})

# Create a long-form dataframe with a row for every day and for
# every species.
dfa = df.groupby(["Date", "scientificName"]).size().reset_index()
dfa = dfa.rename({0: "nobs"}, axis=1)
dr = pd.date_range(df["Date"].min(), df["Date"].max())
species = pd.Series(list(species), name="scientificName")
db = pd.DataFrame(dr).merge(species, how="cross")
db.columns = ["Date", "scientificName"]
dc = db.merge(dfa, on=["Date", "scientificName"], how="left")
dc["nobs"] = dc["nobs"].fillna(0)

# Include some additional date information.
dc["year"] = [x.year for x in dc["Date"]]
dc["month"] = [x.month for x in dc["Date"]]
dc["dayOfYear"] = [x.dayofyear for x in dc["Date"]]

df = dc.sort_values(["scientificName", "Date"])
df.to_csv(os.path.join(pa, "plants_occurrences.csv.gz"), index=None)

dz = dz.reset_index()
dz = dz.sort_values(by="scientificName")
dz.to_csv(os.path.join(pa, "plants_locations.csv.gz"), index=None)
