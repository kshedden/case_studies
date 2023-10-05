library(ggplot2)
library(stringr)
library(purrr)

source("prep.R")

# Get the demographics data
mv = births %>% group_by(FIPS) %>% summarize(births_mean=mean(Births), births_var=var(Births))
mv = mv %>% mutate(log_births_mean=log(births_mean), log_births_var=log(births_var))
demog = demog %>% mutate(across(where(anyNA), ~ replace_na(., 0)))
demog[,2:dim(demog)[2]] = sqrt(demog[,2:dim(demog)[2]])
va = colnames(demog)[2:dim(demog)[2]]
demog = demog %>% filter(FIPS %in% unique(mv$FIPS))
demog = demog %>% mutate_at(va, scale, scale=FALSE)

demogx = demog[,2:dim(demog)[2]]
ss = svd(demogx)

u = ss$u
s = ss$d
v = ss$v

#alpha = 1 # distance interpretation for rows (counties)
#alpha = 0 # distance interpretation for columns (demographic categories)
alpha = 0.5
uu = u %*% diag(s^alpha)
vv = v %*% diag(s^(1-alpha))

iqr = function(x){ quantile(x, 0.75) - quantile(x, 0.25) }

# Remove outlier counties to make the plot easier to read
for (j in 1:2) {
    ii = abs(uu[,j] - median(uu[,j])) < 4*iqr(uu[,j])
    uu = uu[ii,]
}

cn = names(demogx)
cx = str_split(cn, "_")

pdf("biplots_r.pdf")

ud = data.frame(u1=uu[,1], u2=uu[,2])

vd = data.frame(v1=vv[,1], v2=vv[,2], race=cx%>%map(1)%>%as_vector, eth=cx%>%map(2)%>%as_vector,
                sex=cx%>%map(3)%>%as_vector)

# To reduce overplotting, plot each sex separately.
for (sex in c("F", "M")) {
    vd1 = vd %>% filter(sex==sex)

    # Plot grey points for the counties
    plt = ggplot(data=ud, aes(x=u1, y=u2)) + geom_point(alpha=0.2)

    # Plot colored points for the demographic variables
    plt = plt + geom_line(data=vd1, aes(x=v1, y=v2, group=interaction(race, eth), color=interaction(race, eth)))

    plt = plt + ggtitle(ifelse(sex=="F", "Female", "Male"))

    print(plt)
}

dev.off()
