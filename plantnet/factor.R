library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

pa = "/home/kshedden/mynfs/data/Teaching/plantnet"
df = read_csv(file.path(pa, "plants_occurrences.csv.gz"))
dz = read_csv(file.path(pa, "plants_locations.csv.gz"))

pdf("factor_r_plots.pdf")

# Pivot the data so that the species are in the columns
# and the dates are in the rows.
dx = df %>% select(Date, scientificName, nobs) %>% pivot_wider(Date, names_from=scientificName, values_from=nobs)
dx = dx %>% filter(Date >= "2018-01-01")
dates = dx$Date
dx = dx %>% select(-Date)
species = names(dx)
dx = as.matrix(dx)

# Variance stabilizing transformation
dx = sqrt(dx)

speciesmeans = colMeans(dx)
datemeans = rowMeans(dx)

# Double center the data
dx = dx - mean(dx)
n = dim(dx)[1]
p = dim(dx)[2]
dx = dx - outer(array(1, n), colMeans(dx))
dx = dx - outer(rowMeans(dx), array(1, p))

# Plot the means
pd = data.frame(x=dates, y=datemeans)
plt = ggplot(aes(x=x, y=y), data=pd) + geom_line()
plt = plt + labs(x="Date", y="Mean")
print(plt)

stopifnot(all(species == dz$scientificName))

# Factor the matrix once, then sort so that the
# species scores are increasing for the first
# factor.
usv = svd(dx)
ii = order(usv$v[, 1])
dx = dx[, ii]
dz = dz[ii,]

# Factor the matrix again.
usv = svd(dx)

# Scree plot
pd = data.frame(x=seq(length(usv$d)), s=usv$d)
plt = ggplot(aes(x=x, y=s), data=pd) + geom_line()
plt = plt + labs(x="SVD component", y="Singular value")
plt = plt + ggtitle("Scree plot")
print(plt)

# Log/log scree plot
s = usv$d
s1 = s[s > 1e-8]
pd = data.frame(x=seq(length(s1)), s=s1)
plt = ggplot(aes(x=log(x), y=log(s)), data=pd) + geom_line()
plt = plt + labs(x="log SVD component", y="log Singular value")
plt = plt + ggtitle("Scree plot (log space)")
print(plt)

for (j in 1:10) {
	# Plot the date factor scores
	n = length(dates)
	dp = data.frame(x=dates, y=usv$u[,j])
	plt = ggplot(aes(x=x, y=y), data=dp) + geom_line()
	plt = plt + labs(x="Date", y=sprintf("Date factor %d score", j))
	plt = plt + ggtitle(sprintf("Factor %d", j))
	print(plt)

	# Plot the species factor scores
	dp = data.frame(x=seq(1, length(usv$d)), y=usv$v[,j])
	plt = ggplot(aes(x=x, y=y), data=dp) + geom_line()
	plt = plt + labs(x="Species", y=sprintf("Species factor %d score", j))
	plt = plt + ggtitle(sprintf("Factor %d", j))
	print(plt)

	# Scatterplot the species factor scores against the mean
	# latitude for the species.
	dp = data.frame(x=dz$decimalLatitude, y=usv$v[,j])
	plt = ggplot(aes(x=x, y=y), data=dp) + geom_point()
	plt = plt + labs(x="Latitude", y=sprintf("Factor %d score", j))
	print(plt)

	# Scatterplot the species factor scores against the mean
	# longitude for the species.
	dp = data.frame(x=dz$decimalLongitude, y=usv$v[,j])
	plt = ggplot(aes(x=x, y=y), data=dp) + geom_point()
	plt = plt + labs(x="Longitude", y=sprintf("Factor %d score", j))
	print(plt)

	# Scatterplot the species factor scores against the mean
	# elvation for the species.
	dp = data.frame(x=dz$elevation, y=usv$v[,j])
	plt = ggplot(aes(x=x, y=y), data=dp) + geom_point()
	plt = plt + labs(x="Elevation", y=sprintf("Factor %d score", j))
	print(plt)
}

dev.off()
