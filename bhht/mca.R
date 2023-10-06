library(dplyr)
library(FactoMineR)
library(tidyr)
library(readr)

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"

da = read_csv(file.path(pa, "cross-verified-database.csv.gz"))
da = da %>% rename(occ=level1_main_occ, reg=un_region, sex=gender)

dx = da %>% select(birth, occ, sex, reg)
dx = dx %>% drop_na()
dx = dx %>% filter(birth >= 1500)

# Create a "century of birth" variable
dx = dx %>% mutate(bcen = round(birth, -2))
dx = dx %>% select(bcen, occ, sex, reg)
dx = dx %>% mutate(bcen = as.factor(bcen))

# Remove small groups that are difficult to interpret
# due to low precision or high confounding.
dx = dx %>% filter(occ != "Other" & occ != "Missing" & sex != "Other")

mm = MCA(dx, ncp=3, graph=FALSE)

pdf("bhht_r_mca.pdf")

plt = plot(mm, axes=c(1, 2), invisible="ind")
print(plt)

plt = plot(mm, axes=c(1, 3), invisible="ind")
print(plt)
dev.off()
