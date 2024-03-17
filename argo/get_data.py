"""
Download ARGO data
"""

import os, requests

# Store the raw files in this location.  These files
# can be deleted after the "prep.py" script is run.
tpath = "/scratch/stats_dept_root/stats_dept1/kshedden/argo/python"

# The file name template.
target = "%s/argo/raw/OOOO_ocean/YYYY/MM" % tpath

# Get the data from this web site
url = "https://data-argo.ifremer.fr/geo/OOOO_ocean/YYYY/MM/FF"

def getdata(ocean, year, month, day):

    fl = "%4d%02d%02d_prof.nc" % (year, month, day)
    url1 = url.replace("OOOO", ocean)
    url1 = url1.replace("YYYY", "%4d" % year)
    url1 = url1.replace("MM", "%02d" % month)
    url1 = url1.replace("DD", "%02d" % day)
    url1 = url1.replace("FF", fl)

    target1 = target.replace("OOOO", ocean)
    target1 = target1.replace("YYYY", "%4d" % year)
    target1 = target1.replace("MM", "%02d" % month)

    os.makedirs(target1, exist_ok=True)
    target1 = "%s/%s" % (target1, fl)
    response = requests.get(url1)
    open(target1, "wb").write(response.content)

def get_year(ocean, year):
    for month in range(1, 13):
        print("%s..." % month)
        if month in [4, 6, 9, 11]:
            days = 30
        elif month == 2:
            days = 29 if year % 4 == 0 else 28
        else:
            days = 31

        for day in range(1, days+1):
            print("%s." % day, end="", flush=True)
            getdata(ocean, year, month, day)
        print("")

get_year("pacific", 2020)
