library(dplyr)
library(readr)
library(lubridate)

# The raw file from plantnet should be located here.
pa = "/home/kshedden/myscratch/plantnet"

# This is the raw data file.  If your file name does not
# match it this needs to be changed.
fn = "0140072-220831081235567.csv.gz"

fn = file.path(pa, fn)
df = read_delim(fn, delim="\t")

# Keep the 200 most common species.
ds = df %>% group_by(scientificName) %>% summarize(nrow=n()) %>% arrange(desc(nrow)) %>% head(200)
species = ds$scientificName
ii = df$scientificName %in% species
df = df[ii,]
df = df %>% mutate(Date = make_datetime(year, month, day))
df = df %>% select(scientificName, Date, elevation, decimalLatitude, decimalLongitude)
write_csv(df, file.path(pa, "short.csv.gz"))

circmean = function(x) {
    x = pi * x / 180
    s = mean(sin(x))
    c = mean(cos(x))
    return(180 * atan2(s, c) / pi)
}

# Mean latitude, longitude, and elevation for each species.
dz = df %>% group_by(scientificName) %>% summarize(decimalLatitude=mean(decimalLatitude),
              elevation=mean(elevation, na.rm=TRUE), decimalLongitude=circmean(decimalLongitude))

# Create a long-form dataframe with a row for every day and for
# every species.
dfa = df %>% group_by(Date, scientificName) %>% summarize(nobs=n())
dt = seq(min(df$Date), max(df$Date), by="day")
db = expand.grid(species, dt)
colnames(db) = c("scientificName", "Date")
db = left_join(db, dfa, on=c("scientificName", "Date"))
db$nobs[is.na(db$nobs)] = 0

# Include some additional date information.
db = db %>% mutate(year = year(Date), month = month(Date), dayOfYear = yday(Date))

df = db %>% arrange(scientificName, Date)
write_csv(df, file.path(pa, "plants_occurrences.csv.gz"))

dz = dz %>% arrange(scientificName)
write_csv(dz, file.path(pa, "plants_locations.csv.gz"))
