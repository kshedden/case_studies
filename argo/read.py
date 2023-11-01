from pathlib import Path
import pandas as pd
import numpy as np

# This script loads the data into several arrays
# and vectors.  It is a separate script so that
# it can be shared by several analysis scripts.

# This must match the variable "qpath" in the 'prep.py'
# script.
qpath = Path("/home/kshedden/data/Teaching/argo/python")

lat = np.loadtxt(qpath / "lat.csv.gz", skiprows=1)
lon = np.loadtxt(qpath / "lon.csv.gz", skiprows=1)

date = np.loadtxt(qpath / "date.csv.gz", dtype="str", skiprows=1)
date = pd.to_datetime(date)
day = date - date.min()
day = np.asarray([x.days for x in day])

temp = np.loadtxt(qpath / "temp.csv.gz", skiprows=1, delimiter=",")
pressure = np.loadtxt(qpath / "pressure.csv.gz", skiprows=1)
psal = np.loadtxt(qpath / "psal.csv.gz", skiprows=1, delimiter=",")
