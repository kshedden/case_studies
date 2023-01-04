# Use ridge regression to predict X-ray flux at a given time
# distance into the future, using a block of consecutive values.

library(dplyr)
library(ggplot2)

pdf("autoreg_r.pdf")
source("read.R")

# Use blocks of size m, and use the first q observations to predict
# the final observation.
m = 1000
q = 200

# The time points of the predictor information relative to the time
# being predicted.
tax = seq(-2*m, -2*(m-q))[1:q]
tax = tax / 60

# Regress y on x using ridge regression, with penalty parameter f.
ridge = function(x, y, f) {
	usv = svd(x)
	u = usv$u
	d = usv$d
	v = usv$v
	g = d^2 / (d^2 + f)
	b = v %*% ((t(u) %*% y) * g)
	return(b)
}

# Prepare a year for autoregression modeling
prep_data = function(year) {
	dx = get_goes(year)
	rr = make_blocks(dx, m, 0)
	tix = rr$time
	flx = rr$flux
	flx = log(1e-8 + flx)
	n = dim(flx)[1]
	x = t(flx[1:q,])
	y = flx[n,]

	# Center the data
	n = dim(x)[1]
	y = y - mean(y)
	x = x - outer(array(1, n), colMeans(x))
	
	return(list(y=y, x=x))
}

pt = prep_data(2017)
y = pt$y
x = pt$x

pt = prep_data(2019)
ytest = pt$y
xtest = pt$x

tr = quantile(ytest, 0.5)
ytestbin = 1*(ytest >= tr)

# Consider how the regression coefficients look for various values
# of the penalty parameter.
for (f in c(0, 0.001, 0.1, 1, 10, 100, 1000, 10000)) {
	b = ridge(x, y, f)
	dp = data.frame(b=b, tax=tax)
	plt = ggplot(aes(x=tax, y=b), data=dp) + geom_line()
	plt = plt + ggtitle(sprintf("f=%.4f", f))
	plt = plt + labs(x="Minutes before current time", y="Coefficient")
	print(plt)
}

dev.off()
