[Plantnet](https://plantnet.org/en) is a platform and mobile app that allows people
to submit photographs of plants, which may then be identified either by
AI or by human contributors.  

The data are available [here](https://www.gbif.org/publisher/da86174a-a605-43a4-a5e8-53d484152cd3),
and specifically we will use the smaller set of more curated data 
[here](https://www.gbif.org/dataset/7a3679ef-5582-4aaa-81f0-8c2545cafc81).  Use the download
link and select the "GBIF annotated archive" file (a direct link is
[here](https://www.gbif.org/occurrence/download?dataset_key=7a3679ef-5582-4aaa-81f0-8c2545cafc81)).
You should select the "simple" download, and may need to register on the site to complete
the download.

The data we have here records the location (latitude, longitude,
elevation) and date for each observation of each plant species. 

After downloading the data file as described above, run one of the "prep"
scripts to prepare analysis files.  You will need to adjust some of the
paths in the prep script.
