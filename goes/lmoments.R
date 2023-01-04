library(dplyr)
library(tidyr)
library(ggplot2)
library(Lmoments)
source("read.R")

pdf("lmoments_R.pdf")

df = get_goes(2017)

dlm = df %>% group_by(Year, Month, Day) %>% group_modify(~ {
    data.frame(Lmoments(.x$Flux1))
})

dlm$L3s = dlm$L3 / dlm$L2
dlm$L4s = dlm$L4 / dlm$L2

v = c("L1", "L2", "L3s", "L4s")

for (j in 1:4) {
    for (k in 1:4) {
        if (j < k) {
            plt = ggplot(aes(x=!!sym(v[j]), y=!!sym(v[k])), data=dlm) + geom_point()
            print(plt)
        }
    }
}

dev.off()
