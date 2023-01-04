# Calculate Kendall tau autocorrelations for blocks of consecutive observations
# taken from the much longer X-ray flux time series.

# The block size can be set below via the variable 'bs'.

library(pcaPP)
library(ggplot2)

source("read.R")

pdf("autocor_R.pdf")

df = get_goes(2017)

# Block size, if bs=4000 the overall time of each block
# is around 2 hours and 15 minutes.
bs = 4000

# Make blocks of 'bs' consecutive time points with
# approximately 2-second spacing.
rr = make_blocks(df, 2000, 0)
tix = rr$time
flx = rr$flux

n = dim(flx)[1]
p = dim(flx)[2]

# Consider autocorrelation at these time lags
d = 25
dlags = floor(seq(1, 200, length.out=10))

# Convert lags to time in minutes
dtime = dlags * 2 / 60

# Calculate these quantiles across blocks of the autocorrelations.
pr = c(0.25, 0.5, 0.75)

qa = array(0, c(length(dlags), length(pr)))
for (i in 1:length(dlags)) {

	d = dlags[i]

	# Get the autocorrelation for each block
	r = array(0, p)
	for (j in 1:p) {
		r[j] = cor.fk(flx[1:(n-d), j], flx[(d+1):n, j])
	}
	qa[i,] = quantile(r[is.finite(r)], pr)
}

# Make a data frame for plotting
qp = data.frame()
for (j in 1:dim(qa)[2]) {
	qq = data.frame(autocor=qa[,j])
	qq$prob = pr[j]
	qq$dlag = dtime
	qp = rbind(qp, qq)
}
qp$prob = as.factor(qp$prob)

plt = ggplot(aes(x=dlag, y=autocor, group=prob, color=prob), data=qp) + geom_line()
plt = plt + labs(x="Lag (minutes)", y="Autocorrelation")
print(plt)

dev.off()
