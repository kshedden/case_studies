import pandas as pd
import numpy as np
import os

# This script loads the data into several arrays
# and vectors.  It is a separate script so that
# it can be shared by several analysis scripts.

# This must match the variable "qpath" in the 'prep.py'
# script.
qpath = "/home/kshedden/data/Teaching/argo/python"

lat = np.loadtxt(os.path.join(qpath, "lat.csv.gz"), skiprows=1)
lon = np.loadtxt(os.path.join(qpath, "lon.csv.gz"), skiprows=1)

date = np.loadtxt(os.path.join(qpath, "date.csv.gz"), dtype="str", skiprows=1)
date = pd.to_datetime(date)
day = date - date.min()
day = np.asarray([x.days for x in day])

temp = np.loadtxt(os.path.join(qpath, "temp.csv.gz"), skiprows=1, delimiter=",")
pressure = np.loadtxt(os.path.join(qpath, "pressure.csv.gz"), skiprows=1)
psal = np.loadtxt(os.path.join(qpath, "psal.csv.gz"), skiprows=1, delimiter=",")
