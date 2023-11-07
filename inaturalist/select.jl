using CodecZlib

# Select these data
const kingdom = "Plantae"
#const class = "Pinopsida"
const class = "Polypodiopsida"

# Source path
pa = "/scratch/stats_dept_root/stats_dept1/kshedden/inaturalist"
fn = joinpath(pa, "0037643-231002084531237.csv.gz")

# Target path
gn = joinpath("/home/kshedden/data/Teaching/inaturalist/$(kingdom)_$(class).csv.gz")

vnames = ["species", "decimalLatitude", "decimalLongitude", "elevation", "eventDate"]

open(GzipDecompressorStream, fn) do io
	open(GzipCompressorStream, gn, "w") do out

		line = readline(io)
		row = split(line, '\t')

		# The positions of the rows of interest
		ii = [findfirst(v->v==x, row) for x in vnames]

		# Write out the header
		write(out, join(row[ii], ","))
		write(out, "\n")

		for r in eachline(io)
			row = split(r, '\t')
			if (row[4] != kingdom) || (row[6] != class)
				continue
			end
			write(out, join(row[ii], ","))
			write(out, "\n")
		end
		flush(out)
	end
end
