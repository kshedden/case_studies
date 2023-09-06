library(R.utils)

base_url = "https://www.ncei.noaa.gov/data/daily-summaries/access"

target_dir = "/home/kshedden/data/Teaching/precip"

files = c("USW00094847.csv", "USW00012839.csv")

for (f in files) {
    url = file.path(base_url, f)
    target = file.path(target_dir, f)
    download.file(url, target)
    gzip(target, overwrite=TRUE)
}
