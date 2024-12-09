import os, requests
import pandas as pd


urls = [
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BMX_J.xpt",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BPX_J.xpt",
    ]

pa = "/home/kshedden/data/Teaching/nhanes"
os.makedirs(pa, exist_ok=True)

for url in urls:
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
