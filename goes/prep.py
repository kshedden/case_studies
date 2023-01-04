import os
import numpy as np
from astropy.io import fits
import pandas as pd

# Location of the fits files obtained using get_data.py
tpath = "/scratch/stats_dept_root/stats_dept1/kshedden/goes/python"

# Location of the csv file produced by this script
qpath = "/home/kshedden/data/Teaching/goes"

# Create a single long-form data file for the given year.
def prep_year(year):

    xl = []
    tpa = os.path.join(tpath, str(year))
    for root, dirs, files in os.walk(tpa):

        for f in files:
            yr = f[4:8]
            yr = int(yr)
            if yr != year:
                continue
            month = int(f[8:10])
            day = int(f[10:12])

            if f.endswith(".fits"):
                fn = os.path.join(root, f)
                try:
                    with fits.open(fn) as io:
                        fld = np.asarray(io["FLUXES"].data)
                        tim = fld["Time"][0,:]
                        flx = fld["Flux"][0,:,:]
                        x = np.hstack((tim[:, None], flx))
                        x = pd.DataFrame(x, columns=["Time", "Flux1", "Flux2"])
                        x["Year"] = year
                        x["Month"] = month
                        x["Day"] = day
                        x = x[["Year", "Month", "Day", "Time", "Flux1", "Flux2"]]
                        xl.append(x)
                except:
                    print("Failed to process: %s" % f)

            xx = pd.concat(xl, axis=0)

    xx = xx.sort_values(by=["Year", "Month", "Day", "Time"])
    xx.to_csv("%s/goes%4d.csv.gz" % (qpath, year), index=None)

prep_year(2017)
prep_year(2019)
