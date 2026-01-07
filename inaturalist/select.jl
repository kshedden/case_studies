using CSV
using DataFrames
using Dates

# Select these data
const kingdom = "Plantae"
const class = "Pinopsida"
#const class = "Polypodiopsida"

# Source path
pa = "/scratch/stats_dept_root/stats_dept1/kshedden/inaturalist"
fn = joinpath(pa, "observations.csv")

# Target path
gn = joinpath("/home/kshedden/data/Teaching/inaturalist/$(kingdom)_$(class).csv.gz")

vnames = ["scientificName", "decimalLatitude", "decimalLongitude", "eventDate", "class"]

for (ix, chunk) in enumerate(CSV.Chunks(fn; ntasks=100))
    chunk = DataFrame(chunk)
    chunk = select(chunk, vnames)
    chunk = filter(r->!ismissing(r[:class]) && r[:class] == class, chunk)
    chunk[:, :eventDate] = [split(x, "T")[1] for x in chunk[:, :eventDate]]
    chunk[!, :eventDate] = Date.(chunk[:, :eventDate])
    CSV.write(gn, chunk; compress=true, append=ix>1)
end

