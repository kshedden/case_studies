import pandas as pd
from pathlib import Path
import numpy as np
from scipy import stats

# Directory in which to place the data files
spath = Path("/home/kshedden/data/Teaching/birthweight/births")

# Years to include
years = range(1981, 1990)

# Load the birth data
dfx = []
for year in years:
    df = pd.read_csv(spath / f"{year}.csv.gz", dtype={"county": str, "smsa": str, "popsize": str})
    dfx.append(df)

# Combine the birth data into one long dataframe
df = pd.concat(dfx, axis=0)
df["loc"] = ["%02d%s" % (a, b) for (a, b) in zip(df["state"], df["county"])]

df["popsize"] = df["popsize"].replace({"9": 4, "Z": np.nan})
df = df.loc[df["popsize"].notnull()] # Drop births to foreign mothers
df["popsize"] = df["popsize"].astype(int)

# Restrict to first born single births
df = df.query("birthorder==1 & plurality==1")

f1 = lambda x: stats.entropy(x.value_counts().values)
f2 = lambda x: np.mean(x != "White")

dr = df.groupby(["year", "loc"]).agg({"momrace": [f1, f2]})
dr.columns = dr.columns.to_flat_index()
dr = dr.rename(columns={("momrace", "<lambda_0>"): "momrace_entropy", ("momrace", "<lambda_1>"): "momrace_nonwhite"})
dr = dr.reset_index()

f1 = lambda x: np.sum(x <= 1500)
f2 = lambda x: np.sum((x > 1500) & (x <= 2500))
f3 = lambda x: np.sum(x > 2500)
f4 = lambda x: pd.isnull(x).mean()
f5 = lambda x: x.dropna().mean()

da = df.groupby(["year", "loc", "sex", "momrace"]).agg({"birthweight": [f1, f2, f3], "momage": np.mean, "dadage": [f4, f5], "popsize": "median"})
da.columns = da.columns.to_flat_index()
da = da.reset_index()
da = da.rename(columns={("birthweight", "<lambda_0>"): "vlb", ("birthweight", "<lambda_1>"): "lb", ("birthweight", "<lambda_2>"): "nb",
                ("popsize", "median"): "popsize_median", ("momage", "mean"): "momage_mean",
                ("dadage", "<lambda_0>"): "dadmissing_mean", ("dadage", "<lambda_1>"): "dadage_mean"})

da = pd.merge(da, dr, on=["year", "loc"], how="left")

da.to_csv(spath / "aggregated.csv.gz", compression="gzip", index=False)
