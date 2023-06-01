using Printf
using DataFrames
using CSV

letter = Dict(2007=>"E", 2009=>"F", 2011=>"G", 2013=>"H", 2015=>"I", 2017=>"J")

pa = "/home/kshedden/data/Teaching/nhanes"

function read_data(year)

    fn = ["DEMO", "BMX", "BPX"]
    dl = []
    for f in fn
        di = @sprintf("%4d-%4d", year, year+1)
        letter = "ABCDEFGHIJ"[1 + div(year - 1999, 2)]
        letter = year > 1999 ? "_$(letter)" : ""
        dx = open(joinpath(pa, di, "$(f)$(letter).csv.gz")) do io
            CSV.read(io, DataFrame)
        end
        push!(dl, dx)
    end

    df = dl[1]
    df = leftjoin(df, dl[2], on=:SEQN)
    df = leftjoin(df, dl[3], on=:SEQN)

    df = filter(r->r.RIDAGEYR>=18, df)

    return df
end

function read_multi()

    va = [:SEQN, :RIDAGEYR, :RIAGENDR, :BPXSY1, :BMXBMI, :BMXHT, :BMXWT, :BMXLEG,
          :BMXARML, :BMXARMC, :BMXWAIST]

    da = []
    for y in 1999:2:2017
        df = read_data(y)
        df = df[:, va]
        df[!, :YOB] = y .- df[:, :RIDAGEYR]
        push!(da, df)
    end

    da = vcat(da...)
    return da
end
