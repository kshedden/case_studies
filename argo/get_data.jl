using Printf

# Store the raw files in this location.  These files
# can be deleted after the "prep.jl" script is run.
tpath = "/scratch/stats_dept_root/stats_dept1/kshedden/argo/julia"

# The file name template.
target = "$(tpath)/argo/raw/OOOO_ocean/YYYY/MM"

# Get the data from this web site
url = "https://data-argo.ifremer.fr/geo/OOOO_ocean/YYYY/MM/FF"

function getdata(ocean, year, month, day)

    fl = @sprintf("%4d%02d%02d_prof.nc", year, month, day)
    url1 = replace(url, "OOOO"=>ocean, "YYYY"=>@sprintf("%4d", year),
                   "MM"=>@sprintf("%02d", month), "DD"=>@sprintf("%02d", day),
                   "FF"=>fl)
    target1 = replace(target, "OOOO"=>ocean, "YYYY"=>@sprintf("%4d", year),
                      "MM"=>@sprintf("%02d", month))

    mkpath(target1)
    target1 = @sprintf("%s/%s", target1, fl)
    download(url1, target1)
end

function get_year(ocean, year)
    for month in 1:12
        print("$(month)...")
        days = if month in [4, 6, 9, 11]
            30
        elseif month == 2
            year % 4 == 0 ? 29 : 28
        else
            31
        end
        for day in 1:days
            print("$(day).")
            getdata(ocean, year, month, day)
        end
        println("")
    end
end

get_year("atlantic", 2020)
