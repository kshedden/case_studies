import os, requests
import pandas as pd

year = 2015

code = {2009: "F", 2011: "G", 2013: "H", 2015: "I", 2017: "J"}[year]


urls = [
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/YYYY/DataFiles/DEMO_CC.xpt",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/YYYY/DataFiles/BMX_CC.xpt",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/YYYY/DataFiles/BPX_CC.xpt",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/YYYY/DataFiles/BIOPRO_CC.xpt",
    ]

pa = "/home/kshedden/data/Teaching/nhanes"
os.makedirs(pa, exist_ok=True)

for urlx in urls:
    url = urlx.replace("YYYY", str(year)).replace("CC", code)
    response = requests.get(url)
    _, tail = os.path.split(url[8:])
    target = os.path.join(pa, tail)
    open(target, "wb").write(response.content)

for root, dirs, files in os.walk(pa):
    for file in files:
        if not file.endswith(".xpt"):
            continue
        da = pd.read_sas(os.path.join(root, file))
        da.to_csv(os.path.join(root, file.replace(".xpt", ".csv.gz")), index=None)
