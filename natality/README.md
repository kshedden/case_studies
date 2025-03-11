# US Natality Data

To access the birth count data, visit
[https://wonder.cdc.gov/natality-current.html](https://wonder.cdc.gov/natality-current.html).
You will need to click "I agree" at the bottom of the page to accept
the terms of usage.  This will take you to the page where you can
download the data.

On the data download page, in section 1, change the option for "Group
Results By" to "County", then in section 4, select year 2011.  Then
click "Send" to get the data, and on the next page click "Export" to
download, and rename the downloaded file to "2011.txt".

Repeat the above steps for 2012, 2013, ..., 2020, to get 10 years of
data.  Then compress all the files using gzip.

Next visit
[https://seer.cancer.gov/popdata/download.html](https://seer.cancer.gov/popdata/download.html)
to obtain demographic data about the counties.  The relevant file is
under the table labeled "County-Level Population Files - 20 Age Groups", in the column labeled
"1990-2023 4 Expanded Races by Origin", in the row labeled "All States
Combined".  You should download the file in gzip/csv (.gz) format.
[Here](https://seer.cancer.gov/popdata/yr1990_2023.20ages/us.1990_2023.20ages.adjusted.txt.gz) is
a direct link.

The layout of the county demographics data is at
[https://seer.cancer.gov/popdata/popdic.html](https://seer.cancer.gov/popdata/popdic.html).

Then visit [https://www.ers.usda.gov/data-products/rural-urban-continuum-codes.aspx](https://www.ers.usda.gov/data-products/rural-urban-continuum-codes.aspx)
and download the 2013 "Rural-Urban Continuum Codes".

We will use the [Area Deprivation index](https://www.neighborhoodatlas.medicine.wisc.edu), which can be obtained
[here](https://www.dropbox.com/scl/fo/jo5ika5l84kr2wtw708vx/h?rlkey=8iy0w0mggomucx3hlcfmht88s&st=h2cnjmll&dl=0).

Once you have downloaded all these files, you can run one of the 'prep.X' scripts to
prepare the data.  As usual, you will have to review these files to understand what they
do, and set some configuration options for your local system.
