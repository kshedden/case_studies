library(dplyr)
library(FactoMineR)
library(readr)

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"

da = read_csv(file.path(pa, "cross-verified-database.csv.gz"))
da = rename(da, "occ"="level1_main_occ")
da = rename(da, "reg"="un_region")
da = rename(da, "sex"="gender")

dx = da[, c("birth", "occ", "sex", "reg")]
dx = dx[complete.cases(dx),]
dx = dx[dx$birth >= 1500,]

# Create a "century of birth" variable
dx$bcen = round(dx$birth, -2)
dx = dx[, c("bcen", "occ", "sex", "reg")]
dx$bcen = as.factor(dx$bcen)

# Remove small groups that are difficult to interpret
# due to low precision or high confounding.
dx = filter(dx, occ != "Other")
dx = filter(dx, occ != "Missing")
dx = filter(dx, sex != "Other")

mm = MCA(dx, ncp=3, graph=FALSE)

pdf("bhht_r_mca.pdf")

plt = plot(mm, axes=c(1, 2), invisible="ind")
print(plt)
dev.off()

plt = plot(mm, axes=c(1, 3), invisible="ind")
print(plt)
dev.off()
