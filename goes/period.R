library(dplyr)
library(lomb)

pdf("period_r.pdf")

source("read.R")

ii = sample(1:dim(df)[1], 10000, replace=F)

period = 10
w = 1 / period
flux1 = cos(w*2*pi*df$Time) + 0.1*rnorm(dim(df)[1])
m = lsp(flux1[ii], df$Time[ii], type="frequency", from=0.01, to=2)
plot(m)

#m = lsp(df$Flux1, df$Time, type="frequency", from=2*1e-4, to=1e-3)
#plot(m)

dev.off()
