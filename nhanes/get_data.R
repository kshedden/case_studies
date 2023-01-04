library(haven)
library(readr)
library(R.utils)

urls = c("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT",
         "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BMX_J.XPT",
         "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPX_J.XPT")

pa = "/home/kshedden/data/Teaching/nhanes"
mkdirs(pa)

for (url in urls) {
    tail = strsplit(url, "/")[[1]]
    tail = tail[length(tail)]
    target = sprintf("%s/%s", pa, tail)
    download.file(url, target)
}

for (file in list.files(pa)) {
    if (!endsWith(file, ".XPT")) {
        next
    }
    f = file.path(pa, file)
    da = read_xpt(f)
    write_csv(da, gsub(".XPT", ".csv.gz", f))
}
