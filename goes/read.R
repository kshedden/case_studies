library(readr)
library(dplyr)
library(lubridate)

# Location of the csv file produced by prep.py
qpath = "/home/kshedden/data/Teaching/goes"

get_goes = function(year) {
	df = read_csv(sprintf("%s/goes%4d.csv.gz", qpath, year))
	df["Date"] = make_date(df$Year, df$Month, df$Day)
	return(df)
}

make_blocks = function(df, m, d) {

	# Create a measure of time with units seconds that 
	# increases indefinitely without resetting to zero.
	df[,"yday"] = yday(df$Date)
	df = df %>% filter(Time > 0)
	df = df %>% mutate(Timex = yday * 60 * 60 * 24 + Time)
	df = df %>% arrange(df, Timex)
	ti = df$Timex
	fl = df$Flux1

	# Restructure the series into arrays with each column being
	# a block
	q = floor(length(ti) / m)
	n = q * m
	ti = ti[1:n]
	fl = fl[1:n]
	g = floor(n / m)
	tix = array(ti, c(m, g))
	flx = array(fl, c(m, g))

	# Time difference within block
	td = tix[m,] - tix[1,]

	# Exclude the blocks that contain skips
	ii = abs(td - median(td)) < 1
	tix = tix[, ii]
	flx = flx[, ii]

	if (d > 0) {
		for (j in 1:d) {
			flx = diff(flx, dims=1)
		}
	}

	return(list(time=tix, flux=flx))
}
