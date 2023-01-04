# Prepare the birth count and other related data for analysis.

library(readr)
library(dplyr)
library(tidyr)
library(readxl)

# Path to the data files
pa = "/home/kshedden/data/Teaching/natality"

# Subset the demographics file to 2016.  This only needs to
# be run once.  Change the next line to read 'if (TRUE)' to
# run this section of code.
if (FALSE) {
    f = sprintf("%s/us.1990_2020.19ages.txt.gz", pa)
    g = sprintf("%s/2016ages.txt.gz", pa)
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

# Create a long form version of the births
dl = list()
for (y in 2016:2020) {

    # If you did not gzip the data files, then remove the
    # .gz below.
    fn = sprintf("%s/%4d.txt.gz", pa, y)
    ct = cols(Births=col_double())
    da = read_tsv(fn, col_types=ct, show_col_types=F)
    da = da[,c("County Code", "Births")]
    da = da[complete.cases(da),]
    da$year = y
    dl[[length(dl)+1]] = da
}
births = rbind(dl[[1]], dl[[2]], dl[[3]], dl[[4]], dl[[5]])
births = rename(births, "FIPS"="County Code")

# Read the demographics for 2016.  It is a fixed-width format
# file.
start = c(1, 5, 7, 9, 12, 14, 15, 16, 17, 19)
end = c(4, 6, 8, 11, 13, 14, 15, 16, 18, 26)
na = c("Year", "State", "StateFIPS", "CountyFIPS", "Registry", "Race", "Origin", "Sex", "Age", "Population")
pos = fwf_positions(start, end, na)
fn = sprintf("%s/2016ages.txt.gz", pa)
ct = cols(Population=col_number())
demog = read_fwf(fn, pos, ct, show_col_types=FALSE)

# Create a FIPS code that matches the FIPS code in the birth data
demog = mutate(demog, FIPS=paste(StateFIPS, CountyFIPS, sep=""))
demog = select(demog, FIPS, Race, Origin, Sex, Age, Population)

# Recode some variables to more interpretable text labels
demog = mutate(demog, Sex=recode(Sex, `2`="F", `1`="M"))
demog = mutate(demog, Origin=recode(Origin, `0`="N", `1`="H"))
demog = mutate(demog, Race=recode(Race, `1`="W", `2`="B", `3`="N", `4`="A"))

# The overall population per county
pop = group_by(demog, FIPS) %>% summarize(Population=sum(Population))

# Pivot to put the age bands in the columns
demog = pivot_wider(demog, id_cols="FIPS", names_from=c("Race","Origin","Sex","Age"), values_from="Population")

# Get the Rural/Urban Continuity Codes (RUCC)
rucc = read_excel(sprintf("%s/ruralurbancodes2013.xls", pa), sheet="Rural-urban Continuum Code 2013")
rucc = select(rucc, FIPS, RUCC_2013)
