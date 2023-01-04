using DataFrames, CSV

year = 2019

# Location of the csv file produced by this script
qpath = "/home/kshedden/data/Teaching/goes"

df = open(joinpath(qpath, "goes2019.csv.gz")) do io
    CSV.read(io, DataFrame)
end

function make_blocks(ti, fl, m, d)

    q = div(length(ti), m)
    n = q * m
    ti = ti[1:n]
    fl = fl[1:n]

    tix = reshape(ti, (m, div(n, m)))
    flx = reshape(fl, (m, div(n, m)))

    # Time difference within block
    td = tix[end, :] - tix[1, :]

    # Exclude the blocks that contain skips
    ii = abs.(td .- median(td)) .< 1
    tix = tix[:, ii]
    flx = flx[:, ii]

    for j in 1:d
        flx = diff(flx, dims=1)
    end

    return tix, flx
end

