"""
This script takes the raw ARGO data and interpolates it onto a
uniform grid of depths.  Profiles are excluded if they do not
cover the range (minpress, maxpress), or if there are fewer than
100 measured points in the raw profile.
"""

from pathlib import Path
import os, datetime
from scipy.interpolate import interp1d
from netCDF4 import Dataset
import numpy as np

# Path to the NetCDF files here, this path must be the same
# as 'tpath' in the get_data.py or get_data.R script.  After running
# prep.py, the contents of tpath can be deleted.
tpath = Path("/scratch/stats_dept_root/stats_dept1/kshedden/argo/python")
dpath = tpath / Path("argo/raw")

# Store the files produced by this script here.  This must
# agree with the path set in the 'read.py' script.
qpath = Path("/home/kshedden/data/Teaching/argo/python")
os.makedirs(qpath, exist_ok=True)

# Retain only profiles that span this range of pressures (these are dbar
# values and 1 dbar is close to 1 meter).
minpress = 20
maxpress = 1500

# Use a grid of this size
ngrid = 100

# Interpolate all data onto this pressure grid
pressure = np.linspace(minpress, maxpress, ngrid)

def clean_range(x):
    mn = x.getncattr("valid_min")
    mx = x.getncattr("valid_max")
    z = np.asarray(x[:])
    z = np.where(z < mn, 9999, z)
    z = np.where(z > mx, 9999, z)
    return z

def get_raw(fn):
    ncf = Dataset(fn, "r")
    lat = ncf.variables["LATITUDE"]
    lon = ncf.variables["LONGITUDE"]
    pres = ncf.variables["PRES_ADJUSTED"]
    temp = ncf.variables["TEMP_ADJUSTED"]
    psal = ncf.variables["PSAL_ADJUSTED"]

    temp = clean_range(temp)
    psal = clean_range(psal)
    pres = clean_range(pres)

    return lat, lon, pres, temp, psal

def interp_profile(pres, temp, psal):
    ii = (pres != 9999) & (temp != 9999) & (psal != 9999)
    if sum(ii) < 100:
        # Exclude profiles with fewer than 100 measured values
        return None, None

    pres = pres[ii]
    temp = temp[ii]
    psal = psal[ii]
    if pres.min() >= minpress or pres.max() <= maxpress:
        # Exclude profiles that do not cover the specified pressure range
        return None, None

    ii = np.argsort(pres)
    pres = pres[ii]
    temp = temp[ii]
    psal = psal[ii]

    temp1 = interp1d(pres, temp)(pressure)
    psal1 = interp1d(pres, psal)(pressure)

    return temp1, psal1

def get_profiles():

    nskip = 0
    lat, lon, pres, temp, psal, date = [], [], [], [], [], []

    for root, dirs, files in os.walk(dpath):
        for file in files:
            if not file.endswith(".nc"):
                continue

            year = int(file[0:4])
            month = int(file[4:6])
            day = int(file[6:8])
            print(file)

            dt = datetime.date(year, month, day)
            lat1, lon1, pres1, temp1, psal1 = get_raw(Path(root) / file)

            for j in range(pres1.shape[0]):
                temp2, psal2 = interp_profile(pres1[j, :], temp1[j, :], psal1[j, :])
                if temp2 is not None:
                    lat.append(lat1[j])
                    lon.append(lon1[j])
                    temp.append(temp2)
                    psal.append(psal2)
                    date.append(dt)
                else:
                    nskip += 1

    lat = np.asarray(lat)
    lon = np.asarray(lon)
    date = np.asarray(date)
    temp = np.vstack(temp).T
    psal = np.vstack(psal).T

    return lat, lon, date, temp, psal

lat, lon, date, temp, psal = get_profiles()

date = [x.isoformat() for x in date]

# Save all the files
np.savetxt(qpath / "lat.csv.gz", lat, header="Column1", comments="")
np.savetxt(qpath / "lon.csv.gz", lon, header="Column1", comments="")
np.savetxt(qpath / "date.csv.gz", date, fmt="%s", header="Column1", comments="")
np.savetxt(qpath / "pressure.csv.gz", pressure, header="Column1", comments="")
#header = ["Column%d" % j for j in range(temp.shape[1])]
#header = ",".join(header)
np.savetxt(qpath / "temp.csv.gz", temp.T, delimiter=",", comments="", fmt="%.3f")
np.savetxt(qpath / "psal.csv.gz", psal.T, delimiter=",", comments="", fmt="%.3f")
