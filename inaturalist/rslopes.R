library(readr)
library(lme4)
library(dplyr)
library(splines)
library(tidyr)
library(ggplot2)
library(ggrastr)
library(ggfortify)
library(mapdata)

#pclass = "Pinopsida"
pclass = "Polypodiopsida"

pa = sprintf("/home/kshedden/data/Teaching/inaturalist/Plantae_%s.csv.gz", pclass)

da = read_csv(pa)
da = da %>% mutate(eventDate = as.Date(eventDate))

# Count days since January 1, 2010
da = da %>% mutate(date1 = as.Date('2010-01-01'))
da = da %>% mutate(day=difftime(eventDate, date1, units="days")/1000)

# Convert longitude to radians
da = da %>% mutate(lonrad=pi*decimalLongitude/180)

pdf(sprintf("%s_r.pdf", pclass))

# Plot the locations of the occurrences.
world = map_data("world")
plt = ggplot() + ggtitle(pclass)
plt = plt + geom_map(data=world, map=world, aes(long, lat, map_id=region, fill="grey"))
plt = plt + geom_point(data=da, aes(decimalLongitude, decimalLatitude), alpha=0.2)
plt = plt + theme(legend.position="none")
plt = rasterize(plt, dpi=100)
print(plt)

# A multilevel linear regression with random intercepts. 
m0 = lmer(decimalLatitude ~ (1 | species) + day + sin(lonrad) + cos(lonrad) + sin(lonrad/2) + cos(lonrad/2), da)

# A multilevel linear regression with random intercepts and ramdom slopes.
m1 = lmer(decimalLatitude ~ (1 + day | species) + day + sin(lonrad) + cos(lonrad) + sin(lonrad/2) + cos(lonrad/2), da)

# Estimate the marginal mean latitude at Detroit
dp = head(da, 100)
dp = dp %>% mutate(day=seq(min(da$day), max(da$day), length.out=100))
dp = dp %>% mutate(lonrad=-pi*83/180) # Detroit longitude
dp$yy = predict(m0, dp, re.form=NA)

# Plot the marginal mean latitude
plt = ggplot(dp, aes(x=day, y=yy)) + geom_line() + xlab("Day (x100)") + ylab("Marginal mean latitude")
print(plt)

# Get a day 0 value for each species, based on its longitude
dx = da
dx$day = 0
dx$Latitude0 = predict(m1, dx, re.form=NA)
dx = dx %>% group_by(species) %>% summarize(Latitude0=mean(Latitude0))

# Plot species-level linear trends over this time range
day1 = as.numeric(min(da$day))
day2 = as.numeric(max(da$day))

# Use the predicted random effects to predict species-level
# mean latitudes on the first and last day of the dataset.
rr = ranef(m1)$species
rr = as.data.frame(rr)
rr$species = row.names(rr)
row.names(rr) = NULL
rr = left_join(rr, dx, by="species")
rr = rr %>% mutate(Latitude1=Latitude0+day*(day2-day1))
rr = rr %>% select(species, Latitude0, Latitude1)

# Reshape the predictions for plotting
pr = pivot_longer(rr, cols=Latitude0:Latitude1)
pr = rename(pr, day=name)
pr = rename(pr, latitude=value)
pr = pr %>% mutate(name=ifelse(day=="Latitude0", day1, day2))

# Plot species-level trends in latitude
plt = ggplot(pr, aes(x=day, y=latitude, group=species)) + geom_line(alpha=0.5) 
plt = plt + xlab("Day (x100)") + ylab("Species latitude")
print(plt)

dev.off()
