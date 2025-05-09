# Download ARGO data.

library(stringr)
library(lubridate)

# Store the raw files in this location.  These files
# can be deleted after the "prep.py" script is run.
tpath = "/scratch/stats_dept_root/stats_dept1/kshedden/argo/R"
dir.create(tpath, showWarnings=FALSE)

# The file name template.
target = sprintf("%s/argo/raw/OOOO_ocean/YYYY/MM", tpath)

# Get the data from this web site
url = "https://data-argo.ifremer.fr/geo/OOOO_ocean/YYYY/MM/FF"

getdata = function(ocean, year, month, day) {

    fl = sprintf("%4d%02d%02d_prof.nc", year, month, day)
    url1 = str_replace(url, "OOOO", ocean)
    url1 = str_replace(url1, "YYYY", sprintf("%4d", year))
    url1 = str_replace(url1, "MM", sprintf("%02d", month))
    url1 = str_replace(url1, "DD", sprintf("%02d", day))
    url1 = str_replace(url1, "FF", fl)

    target1 = str_replace(target, "OOOO", ocean)
    target1 = str_replace(target1, "YYYY", sprintf("%4d", year))
    target1 = str_replace(target1, "MM", sprintf("%02d", month))

    dir.create(target1, recursive=T, showWarnings=FALSE)
    target1 = sprintf("%s/%s", target1, fl)

    download.file(url1, target1, mode="wb")
}

get_year = function(ocean, year) {
    for (month in 1:12) {
        cat(sprintf("%s...\n", month))
        mdays = days_in_month(make_datetime(year, month, 1))

        for (day in (1:mdays)) {
            cat(sprintf("  %s.\n", day))
            getdata(ocean, year, month, day)
        }
        cat("\n")
    }
}

#get_year("atlantic", 2020)
get_year("pacific", 2020)
