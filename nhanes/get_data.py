import os, requests
import pandas as pd

urls = [
    "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT",
    "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BMX_J.XPT",
    "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPX_J.XPT",
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
        if not file.endswith(".XPT"):
            continue
        da = pd.read_sas(os.path.join(root, file))
        da.to_csv(os.path.join(root, file.replace(".XPT", ".csv.gz")), index=None)
