from pathlib import Path
import pandas as pd
import numpy as np

# This script loads the Argo data into several arrays and vectors.  It is
# a separate script from the notebooks so that it can be shared by several
# notebooks.

# This must match the variable "qpath" in the 'prep.py'
# script.
qpath = Path("/home/kshedden/data/Teaching/argo/python")

# Latitude and longitude
lat = np.loadtxt(qpath / "lat.csv.gz", skiprows=1)
lon = np.loadtxt(qpath / "lon.csv.gz", skiprows=1)

# The day of year that the profile was measured.
date = np.loadtxt(qpath / "date.csv.gz", dtype="str", skiprows=1)
date = pd.to_datetime(date)
day = date.dayofyear

# The pressure values (in decibars) for each point along a profile.
pressure = np.loadtxt(qpath / "pressure.csv.gz", skiprows=1)

# The temperature and salinity values, each column is a profile
# obtained on a particular date and location.
temp = np.loadtxt(qpath / "temp.csv.gz", delimiter=",").T
psal = np.loadtxt(qpath / "psal.csv.gz", delimiter=",").T
