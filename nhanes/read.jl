using DataFrames

pa = "/home/kshedden/data/Teaching/nhanes"

fn = ["DEMO_J", "BMX_J", "BPX_J"]
dl = []
for f in fn
    dx = open(joinpath(pa, f*".csv.gz")) do io
        CSV.read(io, DataFrame)
    end
    push!(dl, dx)
end

df = dl[1]
df = leftjoin(df, dl[2], on=:SEQN)
df = leftjoin(df, dl[3], on=:SEQN)
