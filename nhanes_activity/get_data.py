import os, requests
import pandas as pd

pa = "/home/kshedden/data/Teaching/nhanes"

url = "https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/PAXHR_H.XPT"

#response = requests.get(url)
_, tail = os.path.split(url[8:])
target = os.path.join(pa, "2013-2014", tail)
#open(target, "wb").write(response.content)

da = pd.read_sas(target)
da["SEQN"] = da["SEQN"].astype(int)
da["PAXDAYH"] = pd.to_numeric(da["PAXDAYH"])
da["PAXDAYWH"] = pd.to_numeric(da["PAXDAYWH"])
da.to_csv(target.replace(".XPT", ".csv.gz"), index=None)
