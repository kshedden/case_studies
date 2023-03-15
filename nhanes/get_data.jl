using ReadStatTables, DataFrames, CSV

urls = [
    "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT",
    "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BMX_J.XPT",
    "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPX_J.XPT",
    ]

pa = "/home/kshedden/mynfs/data/Teaching/nhanes"
mkpath(pa)

for url in urls
    _, tail = splitdir(url[9:end])
    target = joinpath(pa, tail)
    download(url, target)
end

for (root, dirs, files) in walkdir(pa)
    for file in files
        fx = splitext(file)
        ext = last(fx)
        if lowercase(ext) != ".xpt"
            continue
        end
        tb = readstat(joinpath(root, file))
        da = DataFrame(tb)
        CSV.write(joinpath(pa, replace(file, ext=>".csv.gz")), da)
    end
end
