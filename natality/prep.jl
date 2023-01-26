using DataFrames, CSV, Printf, CodecZlib

# Path to the data files
pa = "/home/kshedden/mynfs/data/Teaching/natality"

# Create a long form version of the births
dl = []
for y = 2016:2020
    fn = joinpath(pa, @sprintf("%4d.txt.gz", y))
    tm = Dict("County Code" => String, "Births" => Float64)
    da = open(fn) do io
        CSV.read(io, DataFrame, delim = "\t", types = tm, silencewarnings = true)
    end
    da = da[:, ["County", "County Code", "Births"]]
    da = da[completecases(da), :]
    da[:, :year] .= y
    push!(dl, da)
end
births = vcat(dl...)
births = rename(births, "County Code" => "FIPS")

# Subset the demographics file to 2016
if false
    f = joinpath(pa, "us.1990_2020.19ages.txt.gz")
    g = joinpath(pa, "2016ages.txt.gz")
    open(GzipCompressorStream, g, "w") do out
        open(GzipDecompressorStream, f) do inp
            for line in eachline(inp)
                if startswith(line, "2016")
                    write(out, line)
                    write(out, "\n")
                end
            end
        end
    end
end

# Read the demographics for 2016.  It is a fixed-width format
# file.  Unfortunately Julia does not currently have a good
# fixed width file reader so we have to do a lot of tedious
# processing here.
x = [1, 5, 7, 9, 12, 14, 15, 16, 17, 19, 27]
cs = [(x[i], x[i+1] - 1) for i = 1:length(x)-1]
demog = (
    Year = Int[],
    State = String[],
    StateFIPS = String[],
    CountyFIPS = String[],
    Registry = String[],
    Race = String[],
    Origin = String[],
    Sex = String[],
    Age = String[],
    Population = Float64[],
)
open(GzipDecompressorStream, joinpath(pa, "2016ages.txt.gz")) do io
    for line in eachline(io)
        for j = 1:length(demog)
            v = line[cs[j][1]:cs[j][2]]
            T = eltype(demog[j])
            if T <: AbstractString
                push!(demog[j], v)
            elseif T <: Number
                push!(demog[j], parse(T, v))
            end
        end
    end
end
demog = DataFrame(demog)

# Create a FIPS code that matches the FIPS code in the birth data
demog[:, :FIPS] = ["$x$y" for (x, y) in zip(demog.StateFIPS, demog.CountyFIPS)]
demog = demog[:, [:FIPS, :Race, :Origin, :Sex, :Age, :Population]]

# Recode some variables to more interpretable text labels
demog[:, :Sex] = replace(demog[:, :Sex], "1" => "M", "2" => "F")
demog[:, :Origin] = replace(demog[:, :Origin], "0" => "N", "1" => "H")
demog[:, :Race] = replace(demog[:, :Race], "1" => "W", "2" => "B", "3" => "N", "4" => "A")

# The overall population per county
pop = combine(groupby(demog, :FIPS), :Population => sum)
pop = rename(pop, :Population_sum => :Population)

# Pivot to put the age bands in the columns
demog[:, :group] = [
    "$(w)_$(x)_$(y)_$(z)" for (w, x, y, z) in
    zip(demog[:, :Race], demog[:, :Origin], demog[:, :Sex], demog[:, :Age])
]
demog = select(demog, Not([:Race, :Origin, :Sex, :Age]))
demog = unstack(demog, :group, :Population)

na = names(demog)[2:end]
na = [split(x, "_") for x in na]
na = DataFrame(
    Race = [x[1] for x in na],
    Origin = [x[2] for x in na],
    Sex = [x[3] for x in na],
    Age = [x[4] for x in na],
)

# Get the Rural/Urban Continuity Codes (RUCC)
rucc = open("rucc2013.csv.gz") do io
    CSV.read(io, DataFrame, types = Dict(:FIPS => String))
end
