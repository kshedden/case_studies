using Printf
using DataFrames
using CSV

letter = Dict(2007=>"E", 2009=>"F", 2011=>"G", 2013=>"H", 2015=>"I", 2017=>"J")

pa = "/home/kshedden/data/Teaching/nhanes"

id_vars = [:SEQN]
mhealth_vars = [:DPQ010, :DPQ020, :DPQ030, :DPQ040, :DPQ050, :DPQ060, :DPQ070, :DPQ080, :DPQ090,]
demog_vars = [:RIDAGEYR, :RIAGENDR]
body_vars = [:BMXBMI, :BMXHT, :BMXWT, :BMXLEG, :BMXARML, :BMXARMC, :BMXWAIST]
biochem_vars = [:LBXSAL, :LBXSATSI, :LBXSASSI, :LBXSAPSI, :LBXSBU, :LBXSCA, :LBXSCH,
                :LBXSC3SI, :LBXSCR, :LBXSGTSI, :LBXSGL, :LBXSIR, :LBXSLDSI, :LBXSPH, :LBDSTBSI,
                :LBXSTP, :LBXSTR, :LBXSUA, :LBXSNASI, :LBXSKSI, :LBXSCLSI, :LBXSOSSI, :LBXSGB]

function read_data(year)

    fn = ["DEMO", "BIOPRO", "BMX", "BPX", "DPQ"]
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
    for j in 2:length(dl)
        df = leftjoin(df, dl[j], on=:SEQN)
    end

    df = filter(r->r.RIDAGEYR>=18, df)
    df[:, :year] .= year

    return df
end

function read_multi()

    va = vcat(id_vars, [:BPXSY1], demog_vars, mhealth_vars, body_vars, biochem_vars)

    da = []
    for y in 2005:2:2017
        df = read_data(y)
        df = df[:, va]
        df[!, :YOB] = y .- df[:, :RIDAGEYR]
        push!(da, df)
    end

    da = vcat(da...)
    return da
end
