import numpy as np
import pandas as pd
import os

# The raw data file constructed by prep.jl should be available at this path.
pa = "/home/kshedden/myscratch/plantnet"

# Load the raw data
df = pd.read_csv(os.path.join(pa, "short.csv.gz"))
df["Date"] = pd.to_datetime(df["Date"])

# The data are heavily skewed toward the more recent years, optionally
# restrict the analysis to these years.
firstyear = 2010
df = df.loc[df.Date >= pd.to_datetime("%4d-1-1" % firstyear), :]

# Generate some time variables
df["year"] = [x.year for x in df.Date]
df["dayOfYear"] = [x.dayofyear for x in df.Date]

# There are very few records from southern hemisphere
# so remove them.
df = df.loc[df.decimalLatitude >= 0, :]

# Elevation is very skewed, to prevent the extreme values
# from dominating the fit, clamp at 3000 meters.
df["elevation"] = np.clip(df["elevation"], -np.inf, 3000)

# Centered and scaled year (one unit of decode corresponds to 10
# years of time).
meanyear = df["year"].mean()
df["decade"] = (df["year"] - meanyear) / 10
