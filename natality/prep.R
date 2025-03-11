# Prepare the birth count and other related data for analysis.

library(readr)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)

# Path to the data files
pa = "/home/kshedden/data/Teaching/natality"

# Subset the demographics file to 2016.  This only needs to
# be run once.  Change the next line to read 'if (TRUE)' to
# run this section of code.
if (FALSE) {
    f = file.path(pa, "us.1990_2023.20ages.adjusted.txt.gz")
    g = file.path(pa, "2016ages.txt.gz")
    inc = gzfile(f)
    otc = gzfile(g, "w")

    process_chunk = function(chunk, pos) {
        chunk = chunk[startsWith(chunk, "2016")]
        writeLines(chunk, otc)
    }

    read_lines_chunked(f, chunk_size=100000, callback=process_chunk)

    close(inc)
    close(otc)
}

# Create a long form version of the births.
dl = list()
for (y in 2011:2020) {
    fn = file.path(pa, sprintf("%4d.txt.gz", y))
    ct = cols(Births=col_double())
    # This produces warnings
    da = read_tsv(fn, col_types=ct, show_col_types=F)
    da = da %>% select("County Code", "County", Births)
    da = da[complete.cases(da),]
    da$year = y
    dl[[length(dl)+1]] = da
}
births = bind_rows(dl)
births = rename(births, "FIPS"="County Code")
births = births %>% filter(!str_detect(County, "Unidentified"))

# Read the demographics for 2016.  It is a fixed-width format
# file.
start = c(1, 5, 7, 9, 12, 14, 15, 16, 17, 19)
end = c(4, 6, 8, 11, 13, 14, 15, 16, 18, 26)
na = c("Year", "State", "StateFIPS", "CountyFIPS", "Registry", "Race", "Origin", "Sex", "Age", "Population")
pos = fwf_positions(start, end, na)
fn = file.path(pa, "2016ages.txt.gz")
ct = cols(Population=col_number())
demog = read_fwf(fn, pos, ct, show_col_types=FALSE)

# Create a FIPS code that matches the FIPS code in the birth data
demog = demog %>% mutate(FIPS=paste(StateFIPS, CountyFIPS, sep=""))
demog = demog %>% select(FIPS, Race, Origin, Sex, Age, Population)

# Recode some variables to more interpretable text labels
demog = demog %>% mutate(Sex=recode(Sex, `2`="F", `1`="M"))
demog = demog %>% mutate(Origin=recode(Origin, `0`="N", `1`="H"))
demog = demog %>% mutate(Race=recode(Race, `1`="W", `2`="B", `3`="N", `4`="A"))

# The overall population per county
pop = demog %>% group_by(FIPS) %>% summarize(Population=sum(Population))

# Pivot to put the age bands in the columns
demog = demog %>% pivot_wider(id_cols="FIPS", names_from=c("Race","Origin","Sex","Age"), values_from="Population")

# Replace missing values with zero
demog = demog %>% mutate_at(vars(-"FIPS"), ~ replace_na(., 0))

# Get the Rural/Urban Continuity Codes (RUCC)
rucc = read_excel(file.path(pa, "ruralurbancodes2013.xls"), sheet="Rural-urban Continuum Code 2013")
rucc = rucc %>% select(FIPS, RUCC_2013)

# Get the ADI values, need to aggregate from zip code to county code
fname = "US_2022_ADI_Census_Block_Group_v4_0_1.csv.gz"
adi = read_csv(file.path(pa, fname), col_types=cols_only(FIPS=col_character(), ADI_NATRANK=col_double()))
adi = adi %>% mutate(FIPS5=substr(FIPS, 1, 5))
adi = drop_na(adi)
adi = adi %>% group_by(FIPS5) %>% summarize(ADI_NATRANK = median(ADI_NATRANK))
