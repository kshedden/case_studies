# Examine the lifespans of notable people using the BHHT data.
#
# This analysis uses survival analysis methods, allowing us
# to use information from still-living people.

library(dplyr)
library(survival)
library(readr)
library(splines)
library(ggplot2)

pdf("lifespan_survival_R.pdf")

# Change this as needed to point to the directory holding the data file.
pa = "/home/kshedden/mynfs/data/Teaching/bhht"

# Load the dataset.  Use the latin-1 encoding since there is some non-UTF
# data in the file.  Add "nrows=100000" when developing to reduce the run
# time (but use the complete data to get final results).
df = read_csv(file.path(pa, "cross-verified-database.csv.gz"),
              locale=readr::locale(encoding="latin1"))#, n_max=100000)

# Create a lifespan variable (years of life).  It will be missing for people who are currently living.
df = df %>% mutate(lifespan = death - birth)

# Exclude people born before 1500, there is too little data to gain a meaningful
# understanding of the trends in lifespan prior to this year.
dx = df %>% filter(birth >= 1500)
dx = dx %>% select(birth, lifespan, gender, un_region, level1_main_occ)

# There are a small number of people with missing or "Other" gender but it
# is too small of a sample to draw conclusions.
dx = dx %>% filter(gender %in% c("Female", "Male"))

# Drop uninformative occupation codes.
dx = dx %>% filter(!(level1_main_occ %in% c("Missing", "Other")))

# Censor at 2022
censor_year = 2022
dx = dx %>% mutate(clifespan = ifelse(is.na(lifespan), censor_year - birth, lifespan))
dx = dx %>% mutate(died = is.finite(dx$lifespan))

# Now we can drop all rows with missing data
dx = dx %>% select(-lifespan)
dx = dx[complete.cases(dx),]

# A categorical variable indicating the century in which a person was born.
dx = dx %>% mutate(era = floor((birth - 1500) / 100))

# Plot the survival functions for people born in each century
sf = survfit(Surv(clifespan, died) ~ era, dx)
plot(sf, xlab="Age", ylab="Proportion alive", col=1:6, xlim=c(0, 100))
legend(x="topright", c("1500s", "1600s", "1700s", "1800s", "1900s", "2000s"),
       col=1:6, lty=1)

# Fit a proportional hazards regression model
dx = dx %>% mutate(time = sqrt(birth - 1500))
fml = "Surv(clifespan, died) ~ bs(time, 4) + gender + level1_main_occ + un_region"
m0 = coxph(as.formula(fml), dx)

# Plot the estimated baseline cumulative hazard function
bh = basehaz(m0)
plt = ggplot(bh, aes(x=time, y=hazard)) + geom_line() + xlim(0, 100) +
             xlab("Age") + ylab("Cumulative hazard")
print(plt)

# Plot the estimated baseline cumulative hazard function on the log scale
plt = ggplot(bh, aes(x=time, y=log(hazard))) + geom_line() + xlim(0, 100) +
             xlab("Age") + ylab("Log cumulative hazard")
print(plt)

# Plot the estimated baseline hazard function using numerical differentiation
haz = diff(bh$hazard) / diff(bh$time)
haz = data.frame(hazard=haz, time=bh$time[1:(length(bh$time)-1)])
plt = ggplot(haz, aes(x=time, y=hazard)) + geom_line() + xlim(0, 100) +
             xlab("Age") + ylab("Log hazard")
print(plt)

# Fit a sex-stratified proportional hazards regression model
dx = dx %>% mutate(time = sqrt(birth - 1500))
fml = "Surv(clifespan, died) ~ bs(time, 4) + level1_main_occ + un_region + strata(gender)"
m1 = coxph(as.formula(fml), dx)

# Plot the baseline hazard function for each sex
bh = basehaz(m1)
plt = ggplot(bh, aes(x=time, y=log(hazard), color=strata)) + geom_line() + xlim(0, 100) +
             xlab("Age") + ylab("Log cumulative hazard")
print(plt)

dev.off()
