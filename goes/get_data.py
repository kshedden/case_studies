import os, requests
from calendar import monthrange

# Set this path to a directory where raw files can be stored.
# These files can be deleted after running the prep.py script.
tpath = "/scratch/stats_dept_root/stats_dept1/kshedden/goes/python"

url = "https://umbra.nascom.nasa.gov/goes/fits"

# Get data for a specific year, from a specific GOES satellite.
# Browse the url above to see which satellites produce data in
# each year.
def getdata(year, go):
    os.makedirs("%s/%04d" % (tpath, year), exist_ok=True)
    print("%d: " % year, end="")
    for m in range(1, 13):
        print("%d.." % m, sep="", end="", flush=True)
        nday = monthrange(year, m)[1]
        for d in range(1, nday+1):
            fname = "go%02d%04d%02d%02d.fits" % (go, year, m, d)
            url1 = "%s/%04d/%s" % (url, year, fname)
            response = requests.get(url1)
            target = "%s/%04d/%s" % (tpath, year, fname)
            open(target, "wb").write(response.content)
    print("..DONE\n", flush=True)


getdata(2017, 13)
getdata(2019, 14)
