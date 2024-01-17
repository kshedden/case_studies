# Run this script first to download the data.  You will need to configure
# some of the variables below to match the directory layout on your machine.

using CSV
using DataFrames

base_url = "https://www.ncei.noaa.gov/data/daily-summaries/access"

# Configure this.
target_dir = "/home/kshedden/data/Teaching/precip"

# Download files for these locations.
files = ["USW00094847.csv", "USW00012839.csv",]

function main(files)
    for f in files
        df = download(joinpath(base_url, f))
        CSV.write(joinpath(target_dir, "$(f).gz"))
    end
end

main(files)
