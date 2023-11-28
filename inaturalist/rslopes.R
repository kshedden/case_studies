library(readr)
library(lme4)
library(dplyr)
library(splines)
library(tidyr)
library(ggplot2)
library(ggrastr)
library(ggfortify)
library(mapdata)
library(RLRsim)

# Use multilevel regression to assess for the presence of species-specific
# temporal trends in the latitudes of occurrences of species within a class
# of plants.  These trends may be consistent with climate change, especially
# if species are generally found to move away from the equator.

# Select a class of plants
pclass = "Pinopsida"
#pclass = "Polypodiopsida"

pa = sprintf("/home/kshedden/data/Teaching/inaturalist/Plantae_%s.csv.gz", pclass)

da = read_csv(pa)
da = da %>% mutate(eventDate = as.Date(eventDate))
da = da %>% filter(eventDate >= as.Date('2010-01-01'))
da = da %>% select(!elevation) %>% drop_na()

# Count days since January 1, 2010
da = da %>% mutate(date1 = as.Date('2010-01-01'))
da = da %>% mutate(day=difftime(eventDate, date1, units="days")/1000)

# Get the mean latitude for each species
dm = da %>% group_by(species) %>% summarize(species_latitude=mean(decimalLatitude))
da = inner_join(da, dm, by=join_by(species))

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

# A multilevel linear regression with random species intercepts.
m1 = lmer(decimalLatitude ~ (1 | species) + day + sin(lonrad) + cos(lonrad) + sin(lonrad/2) + cos(lonrad/2), da)

# A multilevel linear regression with random species intercepts and random day slopes for species.
m2 = lmer(decimalLatitude ~ (1 + day | species) + day + sin(lonrad) + cos(lonrad) + sin(lonrad/2) + cos(lonrad/2), da)

# Conduct a formal test that the day (within species) variance is significantly positive.
# This requires us to fit two alternative models.
m3 = lmer(decimalLatitude ~ (0 + day | species) + day + sin(lonrad) + cos(lonrad) + sin(lonrad/2) + cos(lonrad/2), da)
m4 = lmer(decimalLatitude ~ (0 + day | species) + (1 | species) + day + sin(lonrad) + cos(lonrad) + sin(lonrad/2) + cos(lonrad/2), da)
ss = exactRLRT(m3, m4, m1)

# Estimate the marginal mean latitude at Detroit
dp = head(da, 100)
dp = dp %>% mutate(day=seq(min(da$day), max(da$day), length.out=100))
dp = dp %>% mutate(lonrad=-pi*83/180) # Detroit longitude
dp$yy = predict(m1, dp, re.form=NA)

# Plot the marginal mean latitude
plt = ggplot(dp, aes(x=day, y=yy)) + geom_line() + xlab("Day (x1000)") + ylab("Marginal mean latitude")
print(plt)

# Get a day 0 value for each species, based on its longitude
dx = da
dx$day = 0
dx$Latitude0 = predict(m2, dx, re.form=NA)
dx = dx %>% group_by(species) %>% summarize(Latitude0=mean(Latitude0))

# Plot species-level linear trends over this time range
day1 = as.numeric(min(da$day))
day2 = as.numeric(max(da$day))

# Species random effects (random day slopes and intercepts)
r0 = ranef(m2)$species
r0 = as.data.frame(r0)
r0$species = row.names(r0)
row.names(r0) = NULL
r0 = r0 %>% left_join(group_by(da, species) %>% summarize(species_latitude=mean(species_latitude)))

# Use the predicted random effects to predict species-level
# mean latitudes on the first and last day of the dataset.
rr = left_join(r0, dx, by="species")
rr = rr %>% mutate(Latitude1=Latitude0+day*(day2-day1))
rr = rr %>% select(species, Latitude0, Latitude1)

# Reshape the predictions for plotting
pr = pivot_longer(rr, cols=Latitude0:Latitude1)
pr = rename(pr, day=name, latitude=value)
pr = pr %>% mutate(day=recode(day, Latitude0='Day0', Latitude1='Day1'))
pr = pr %>% mutate(name=ifelse(day=="Day0", day1, day2))

# Plot species-level trends in latitude
plt = ggplot(pr, aes(x=day, y=latitude, group=species)) + geom_line(alpha=0.5)

# Plot the population averaged trend line in red.
ang = -pi*83/180
f = fixef(m2)
a = f[["(Intercept)"]] + f[["sin(lonrad)"]]*sin(ang) + f[["cos(lonrad)"]]*cos(ang) +
    f["sin(lonrad/2)"]*sin(ang/2) + f["cos(lonrad/2)"]*cos(ang/2)
b = f[["day"]]
plt = plt + geom_abline(aes(intercept=a, slope=b, color="red"))

plt = plt + xlab("Day (x1000)") + ylab("Species latitude")
print(plt)

dev.off()
