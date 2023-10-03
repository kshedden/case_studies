using Printf
using Downloads
using ReadStatTables
using DataFrames
using CSV
using CodecZlib

pa = "/home/kshedden/data/Teaching/nhanes"
mkpath(pa)

files = ["DEMO", "BMX", "BPX", "BIOPRO", "DPQ"]

# Standard biochemistry starts in 2005
# Depression screener starts in 2005
years = collect(2005:2:2017)

function do_download(years)
    # Download the files.
    for y in years
        di = @sprintf("%4d-%4d", y, y+1)
        mkpath(joinpath(pa, di))
        letter = "ABCDEFGHIJ"[1 + div(y - 1999, 2)]
        for f in files
            g = "$(f)_$(letter).XPT"

            # Prior to 2005, BIOPRO has different names
            if y == 1999 && f == "BIOPRO"
                g = "LAB18_$(letter).XPT"
            elseif y < 2005 && f == "BIOPRO"
                g = "L40_$(letter).XPT"
            end

            # The first wave doesn't follow the naming pattern
            g = replace(g, "_A"=>"")
            h = joinpath(di, g)
            s = "https://wwwn.cdc.gov/Nchs/Nhanes/$(h)"
            println(s)
            Downloads.download(s, joinpath(pa, h))
            r = read(joinpath(pa, h), 10)
            if startswith(String(r), "<!DOCTYPE")
                warning(@sprintf("Unable to read: %s\n", joinpath(pa, h)))
                rm(joinpath(pa, h))
            end
        end
    end
end

do_download(years)

# Convert the files from SAS transport (XPT) to csv.
for (root, dirs, files) in walkdir(pa)
    for file in files
        f = joinpath(root, file)
        if !endswith(lowercase(f), ".xpt")
            continue
        end
        println(f)
        da = readstat(f) |> DataFrame
        open(GzipCompressorStream, replace(f, ".XPT"=>".csv.gz"), "w") do io
            CSV.write(io, da)
        end
    end
end
