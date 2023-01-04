using Printf, Interpolations, UnicodePlots, NetCDF, Dates, Tables, CodecZlib, CSV

# Store the NetCDF files here.  This rquires quite a bit of
# space but can be deleted as soon as this script runs to
# completion.
tpath = "/scratch/stats_dept_root/stats_dept1/kshedden/argo/julia"
dpath = "$(tpath)/argo/raw"

# Store the files produced by this script here
qpath = "/home/kshedden/data/Teaching/argo/julia"
mkpath(qpath)

minpress = 100
maxpress = 1500

# Interpolate all data onto this pressure grid
pressure = collect(range(minpress, maxpress, 100))

function clean_range(x, fn, vn)
    mn = ncgetatt(fn, vn, "valid_min")
    mx = ncgetatt(fn, vn, "valid_max")
    x[x .< mn] .= 9999
    x[x .> mx] .= 9999
    return x
end

function get_raw(fn)
    lat = ncread(fn, "LATITUDE")
    lon = ncread(fn, "LONGITUDE")
    pres = ncread(fn, "PRES_ADJUSTED")
    temp = ncread(fn, "TEMP_ADJUSTED")
    psal = ncread(fn, "PSAL_ADJUSTED")

    temp = clean_range(temp, fn, "TEMP")
    psal = clean_range(psal, fn, "PSAL")
    pres = clean_range(pres, fn, "PRES")

    return lat, lon, pres, temp, psal
end

function interp_profile(pres, temp, psal)
    ii = (pres .!= 9999) .& (temp .!= 9999) .& (psal .!= 9999)
    if sum(ii) < 100
        return nothing, nothing
    end
    pres = pres[ii]
    temp = temp[ii]
    psal = psal[ii]
    if minimum(pres) >= minpress || maximum(pres) <= maxpress
        return nothing, nothing
    end
    ii = sortperm(pres)
    pres = pres[ii]
    temp = temp[ii]
    psal = psal[ii]

    temp1 = linear_interpolation(pres, temp).(pressure)
    psal1 = linear_interpolation(pres, psal).(pressure)
    return temp1, psal1
end

function get_profiles()

    nskip = 0
    lat, lon, pres, temp, psal, date = [], [], [], [], [], []

    for (root, dirs, files) in walkdir(dpath)
        for file in files
            println(file)
            year = parse(Int, file[1:4])
            month = parse(Int, file[5:6])
            day = parse(Int, file[7:8])
            dt = Date(year, month, day)

            lat1, lon1, pres1, temp1, psal1 = get_raw(joinpath(root, file))

            for j in 1:size(pres1, 2)
                temp2, psal2 = interp_profile(pres1[:, j], temp1[:, j], psal1[:, j])
                if !isnothing(temp2)
                    push!(lat, lat1[j])
                    push!(lon, lon1[j])
                    push!(temp, temp2)
                    push!(psal, psal2)
                    push!(date, dt)
                else
                    nskip += 1
                end
            end
        end
    end
    println("nskip=$(nskip)")

    return lat, lon, date, hcat(temp...), hcat(psal...)
end

lat, lon, date, temp, psal = get_profiles()

# The Atlantic ocean is mostly west of 20 degrees longitude.
ii = findall(lon .< 20)
lat = lat[ii]
lon = lon[ii]
date = date[ii]
temp = temp[:, ii]
psal = psal[:, ii]

open(GzipCompressorStream, joinpath(qpath, "lat.csv.gz"), "w") do io
    CSV.write(io, Tables.table(lat))
end

open(GzipCompressorStream, joinpath(qpath, "lon.csv.gz"), "w") do io
    CSV.write(io, Tables.table(lon))
end

open(GzipCompressorStream, joinpath(qpath, "date.csv.gz"), "w") do io
    CSV.write(io, Tables.table(date))
end

open(GzipCompressorStream, joinpath(qpath, "pressure.csv.gz"), "w") do io
    CSV.write(io, Tables.table(pressure))
end

open(GzipCompressorStream, joinpath(qpath, "temp.csv.gz"), "w") do io
    CSV.write(io, Tables.table(temp))
end

open(GzipCompressorStream, joinpath(qpath, "psal.csv.gz"), "w") do io
    CSV.write(io, Tables.table(psal))
end
