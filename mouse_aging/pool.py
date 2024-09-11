# Create a json dictionary containing all data from the Mendeley link below.
# Before running this script, first download the five zip archives from
# the Mendeley link, unzip them, and place the extracted directories into
# a directory called "raw_data", which is contained in the same directory
# as this script.

# https://www.nature.com/articles/s41467-022-34515-y
# https://data.mendeley.com/preview/ypz9zyc9rp?a=09b16f74-4581-48f7-94af-469e01757949

import os
import pandas as pd
import json
import gzip

di = {}

for root, dirs, files in os.walk("raw_data"):
    if root == "raw_data":
        continue
    ro = root.replace("raw_data/", "")
    di[ro] = {}
    for f in files:
        df = pd.read_csv(os.path.join(root, f))
        di[ro][f] = df.to_csv(index=None)

with gzip.open("mouse_data.json.gz", "wt") as io:
    json.dump(di, io)
