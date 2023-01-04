library(readr)
library(dplyr)
library(lubridate)

# The raw data created by prep.R should be available at this path.
pa = "/home/kshedden/myscratch/plantnet"

# Load the raw data.
df = read_csv(file.path(pa, "short.csv.gz"))

# Create some additional date variables.
df = df %>% mutate(year = year(Date))
df = df %>% mutate(dayOfYear = yday(Date))

# There are very few records from the southern hemisphere.
df = df %>% filter(decimalLatitude >= 0)

# The data are heavily skewed toward the more recent years, optionally
# restrict the analysis to these years.
firstyear = 2016
df = df %>% filter(year >= firstyear)

meanyear = mean(df$year)
df = df %>% mutate(decade = (year - meanyear) / 10)
