import pandas as pd
import os

pa = "/home/kshedden/data/Teaching/nhanes/2017-2018"

fn = ["DEMO_J.csv.gz", "BMX_J.csv.gz", "BPX_J.csv.gz"]

da = []
for f in fn:
    dx = pd.read_csv(os.path.join(pa, f))
    da.append(dx)

df = pd.merge(da[0], da[1], how="left", on="SEQN")
df = pd.merge(df, da[2], how="left", on="SEQN")

df["RIAGENDR"] = df["RIAGENDR"].replace([1, 2], ["M", "F"])
df["RIDRETH1"] = df["RIDRETH1"].replace([1, 2, 3, 4, 5], ["MA", "OH", "NHW", "NHB", "Other"])

df = df.loc[df.RIDAGEYR >= 18, :]
