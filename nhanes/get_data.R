library(haven)
library(readr)
library(R.utils)

urls = c("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt",
         "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BMX_J.xpt",
         "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BPX_J.xpt")

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
    print(file.info(f))
    da = read_xpt(f)
    write_csv(da, gsub(".XPT", ".csv.gz", f))
}
