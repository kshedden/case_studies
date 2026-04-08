# United States birth data

This case study uses data on births in the United States, available [here](https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm).
Before around 1988, the US National Center for Health Statistics released data on each individual birth in the US, including
birth weight, parental information, and the county where the birth took place.  More recently, privacy concerns have
led to less data being released.  Therefore, here we focus on the older datasets.

The [prep.py](https://github.com/kshedden/case_studies/birthweights/prep.py) script downloads and processes the data files for several years.  The agg.py script aggregates these data into
counts.
