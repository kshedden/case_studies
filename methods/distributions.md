# Characterizing distributions

Probability distributions are the central object of study in probability theory
and statistics.  Probability theory provides us with several ways to
represent a probability distribution, such as the [probability density
function](https://en.wikipedia.org/wiki/Probability_density_function)
(pdf), [cumulative distribution
function](https://en.wikipedia.org/wiki/Cumulative_distribution_function)
(cdf), and [moment generating
function](https://en.wikipedia.org/wiki/Moment-generating_function)
(mgf).  It also provides us with a means to summarize probability
distributions using numerical characteristics such as the mean and variance.
The field of statistics (as opposed to the field of probability)
focuses on estimating these quantities from data, using estimators such as the
[empirical
cdf](https://en.wikipedia.org/wiki/Empirical_distribution_function),
the [histogram](https://en.wikipedia.org/wiki/Histogram) (an estimator
of the pdf), and the sample mean.

The two most common characteristics used to summarize univariate
distributions of a quantitative value (i.e. probability distributions
on the real line) are
[moments](https://en.wikipedia.org/wiki/Moment_(mathematics)) and
[quantiles](https://en.wikipedia.org/wiki/Quantile).  Both of these
approaches can provide us with measures of *location* (also known as
*centrality* or *central tendency*), and
measures of
[dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion)
(also known as *scale*).  A moment-based measure of location is
the mean, while a quantile-based measure of location is the median.
A moment-based measure of dispersion is the [standard
deviation](https://en.wikipedia.org/wiki/Standard_deviation),
while quantile-based measures of dispersion include the [inter-quartile
range](https://en.wikipedia.org/wiki/Interquartile_range) (IQR) and the
[MAD](https://en.wikipedia.org/wiki/Median_absolute_deviation) (median
absolute deviation).

"Higher order" characteristics of a distribution such as [skew](https://en.wikipedia.org/wiki/Skewness)
and [kurtosis](https://en.wikipedia.org/wiki/Kurtosis) are less commonly
encountered but remain of interest in some settings.

Below we summarize some less familiar characteristics of probability
distributions, and ways to estimate these characteristics from data.

## Heavy-tailed distributions and tail parameter estimation

The tail parameter of a random variable $X$ describes how rapidly the
tail probability $P(X>x)$ (the complementary CDF) converges to zero as
$x$ grows.  In a [heavy tailed
distribution](https://en.wikipedia.org/wiki/Heavy-tailed_distribution),
these probabilities do not shrink exponentially fast, which means that

$$
\lim_{x\rightarrow \infty} \exp(tx) \cdot P(X>x) = \infty
$$

for all $t > 0$.  To understand this definition, suppose that it does
not hold, so there is a value $t>0$ and a constant $c$ such that

$$
\lim_{x\rightarrow \infty} \exp(tx) \cdot P(X>x) = c
$$

Roughly speaking this means that $P(X>x)$ behaves like $\exp(-tx)$ for
large $x$, or if $c = 0$ it means that $P(X>x)$ is dominated by
$\exp(-tx)$ for large $x$.

If a distribution is not heavy-tailed, then it may have a [power
law](https://en.wikipedia.org/wiki/Power_law) tail, meaning that
$P(X>x) \sim x^{-\alpha}$.  The value of $\alpha$ is called the *tail
index*.  To estimate the tail index based on a sample of data
$\\{X_i\\}$, consider the
[order statistics](https://en.wikipedia.org/wiki/Order_statistic)
$X_{(1)}\le X_{(2)} \le \cdots \le X_{(n)}$.  The "Hill slope estimate" of
$\alpha$ is

$$
k^{-1}\sum_{i=0}^{k-1} \log(X_{(n-i)}) - \log(X_{(n-k)}),
$$

where $k$ is a chosen tuning parameter.  We won't be able to justify
this entirely, but note that this is the average in log space of the
differences between the upper $k$ order statistics, and the $n-k^{\rm
th}$ order statistic.

### Pareto tail plots

Rather than estimating the tail index, we can make a plot that captures
the tail behavior of a distribution.  One of the most common ways
to do this proceeds as follows.  Given a sample $X_1, \ldots, X_n$,
form the order statistics $X_{(1)} \le X_{(2)} \le \cdots$.  Recall
that the $j^{\rm th}$ order statistic $X_{(j)}$ is an estimate of
the $j/n$ quantile of the distribution.  The probability of
observing a value greater than the $p^{\rm th}$ quantile of a
distribution is $1-p$, so the probability of observing a value
greater than $X_{(j)}$ is $1 - j/n$.

A *Pareto tail plot* plots $\log(1 - j/n)$ against $\log(X_{(j)})$.  If
the tail is exactly $P(X>x) = c/x^\alpha$ then $\log P(X>x) = \log(c) - \alpha\log(x)$.
This implies that the Pareto tail plot will have a linear pattern with a
slope of $-\alpha$.

## L-moments

A classical statistical moment is defined to be the expected value of
a random variable raised to a power.  For example, $E[X]$ is the
expected value and $E[X^2]$ is the raw second moment.  In practice we
usually work with the *centered moments*, for example $E[(X-EX)^2]$ is
the centered second moment, which is better known as the variance.

In principle, if you know all the moments of a distribution, then you
know everything about the distribution (there are some technical
conditions for this to be literally true, it is based on the
invertibility of the moment generating function).  But this fact is
not very useful in practice because it is nearly impossible to
estimate high order moments $E[(X-EX)^k]$ for large values of $k$.
The sample estimator of this moment is $n^{-1}\sum_{i=1}^n (X_i -
\bar{X})^k$, and this estimator is consistent and asymptotically
unbiased, but it has huge mean squared error for practically realistic
sample sizes.

As noted earlier, many descriptive statistics are either moments or
quantiles.  If high order moments are hard to estimate, perhaps there
is a quantile-based analogue to these moments that is easier to
estimate?  This idea led to the development of
[L-moments](https://en.wikipedia.org/wiki/L-moment) which are linear
combinations of order statistics (order statistics in turn are
estimates of quantiles).

The definition of an L-moment of arbitrary order is complex, so we
focus here on the first four L-moments.

The first L-moment $\lambda_1$ is the same as the usual mean.

The second L-moment of a distribution represented through the random
variable $X$ is defined as $\lambda_2 = (EX_{2:2} - EX_{1:2}) / 2$.
Here, $X_{j:k}$ is defined to be the random variable obtained by
sampling $k$ independent values from the distribution of $X$, and then
taking the $j^{\rm th}$ largest among them.

The third L-moment is

$$
\lambda_3 = (EX_{3:3} - 2EX_{2:3} + EX_{1:3}) / 3.
$$

Finally, the fourth L-moment is

$$
\lambda_4 = (EX_{4:4} - 3EX_{3:4} + 3EX_{2:4} + EX_{1:4}) / 4.
$$

Often we work with the standardized third and fourth L-moments,
$\lambda_3^s = \lambda_3/\lambda_2$ and $\lambda_4^s =
\lambda_4/\lambda_2$.  Note that these standardized L-moments are
*scale invariant* meaning that their value is not changed by scaling
the data.  All L-moments except for the first L-moment are
*translation invariant*, meaning that their values are not changed by
adding a constant to all values.  Scale and translation invariance are
important becuase they imply that the result does not depend on the
units or origin of the measurement scale.
