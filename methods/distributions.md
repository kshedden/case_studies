# Characterizing distributions

Probability distributions are the central object of interest in probability theory
and statistics.  There are several ways to
represent a probability distribution, such as the [probability density
function](https://en.wikipedia.org/wiki/Probability_density_function)
(pdf), the [cumulative distribution
function](https://en.wikipedia.org/wiki/Cumulative_distribution_function)
(cdf), and the [moment generating
function](https://en.wikipedia.org/wiki/Moment-generating_function)
(mgf).  There are many useful ways to summarize probability
distributions using numerical characteristics such as the mean and variance.
The field of statistics (as opposed to the field of probability)
focuses on estimating these quantities from data, using estimators such as the
[empirical
cdf](https://en.wikipedia.org/wiki/Empirical_distribution_function),
the [histogram](https://en.wikipedia.org/wiki/Histogram) (an estimator
of the pdf), and the sample mean.

At a high level, the two most common characteristics used to summarize univariate
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
encountered but are of great interest in some settings.

Below we summarize some less familiar characteristics of probability
distributions, and discuss how to estimate these characteristics from data.

## Extremes, heavy-tailed distributions, and tail parameter estimation

Many important research questions involve the frequency of "extreme" events, for
example major earthquakes, large movements in financial markets, people with
extremely long lifespans, or extreme rainfall events.  The study of
extremes naturally leads us to focus on the right tail of a probability
distribution.  In some cases the extremes of interest may involve the left
tail, but in that case we can flip the distribution -- therefore by convention,
methods for studying extremes focus on the right tail.

The tail parameter of a random variable $X$ describes how rapidly the
tail probability $P(X > t)$ (the complementary CDF, CCDF, or *survival function*)
converges to zero as $t$ grows.  In many familiar distributions, the tails
are *exponential*  meaning that $P(X > t) = L(t)\cdot \exp(-t/\mu)$, where
$L(t)$ is a [slowly varying function](https://en.wikipedia.org/wiki/Slowly_varying_function).
For our purposes, we can treat $L(t)$ as a constant, so having exponential tails
implies that $P(L > t) \propto \exp(-t/\mu)$,
for some *scale parameter* $\mu > 0$.
In a [heavy tailed distribution](https://en.wikipedia.org/wiki/Heavy-tailed_distribution),
the tail probabilities do not shrink exponentially fast, which means that
for all $k > 0$,

$$
\lim_{t\rightarrow \infty} \exp(k\cdot t) \cdot P(X > t) = \infty.
$$

Many heavy-tailed distributions have [power law](https://en.wikipedia.org/wiki/Power_law) tail, meaning
that $P(X > t) = L(t) / t^\alpha$, where $\alpha$ is the *tail index*
or *shape parameter*. In such distributions, the $k^{\rm th}$
moment only exists and is finite if $k < \alpha$ (the greater the value of $\alpha$,
the more moments exist).

The simplest family of distributions with power law tails is the
[Pareto distribution](https://en.wikipedia.org/wiki/Pareto_distribution).
The sample space of the Pareto distribution is $[1, \infty)$, and the
CCDF is $P(X > t) = 1/t^\alpha$. This distribution has a tail index of
$\alpha$ as defined above.

If $U$ follows a uniform distribution, then $1/U^{1/\alpha}$ is Pareto.
Alternativey, if $Y$ follows a standard exponential
distribution, then $Z = \exp(Y / \alpha)$ follows a Pareto distribution.

### Exceedances

The Pareto and exponential distributions are simple one-parameter families
and in general will not fit many datasets well.  When focusing on extreme
values we usually don't want to become distracted by the structure of the
center of the distribution.  One way to focus on the tail is to convert
the data to *exceedances*.  This means that we find a parameter $T$ and replace
the dataset $\{X_i\}$ with the dataset $\\{X_i-T | X_i \ge T\\}$.  For example,
if $T$ is appropriately selected then the exceedances may follow a
Pareto distribution, even though the full dataset is a poor fit to the
Pareto model.

### Tail plots

Let $X_{(j)}$ denote the $j^{\rm th}$ order statistic either of our data,
or of the exceedances derived from our data.  Recall that the $j^{\rm th}$
order statistic is the $j^{\rm th}$ sorted value in our data, sorted in
increasing order.  The $j^{\rm th}$ order statistic corresponds
to *probability point* $j/(n+1)$.  That is, $P(X \le X_{(j)}) = j/(n+1)$
and $P(X > X_{(j)}) = 1 - j/(n+1)$.

Suppose that the tail of $X$ is a power law with tail index $\alpha$.
Then we have

$$
P(X > X_{(j)}) = c/X_{(j)}^\alpha = 1 - j/(n+1).
$$

Therefore in log-space we have

$$
\log(1 - j/(n+1)) = \log(x) - \alpha\log(X_{(j)}).
$$

This implies that if we plot the probability points $1 - j/(n+1)$
against the order statistics $X_{(j)}$ in log-space, we obtain
a linear relationship with slope $-\alpha$.  This type of
probability plot can be used both as a means for estimating
$\alpha$, and as a diagnostic for whether the tails actually
follow a power law.

Alternatively, if our tails are exponential we have the
relationship

$$
P(X > X_{(j)}) = c\cdot\exp(-X_{(j)}/\mu) = 1 - j/(n+1),
$$

therefore

$$
\log(1 - j/(n+1)) = \log(c) - X_{(j)}/\mu.
$$

In a *semi-log* plot (log transforming the probability points
but not the order statistics), we obtain a linear relationship
with slope $-1/\mu$.

Using such probability tail plots, we can estimate distributional
parameters ($\alpha$ or $\mu$) by fitting a least squares regression
line to the points in the plot.  The number of points used in the
least squares fit is a tuning parameter that must be selected, typically
in the range $20-200$.  These estimators are convenient, intuitive,
and *distributionally robust* (since the depend on the assumed form of
the tail but do not require a complete specification of a probability
model).  However these estimators may not be very efficient.
Alternative estimators will be discussed
below.

### The Hill estimate of the tail parameter

If a distribution has a power-law tail, we can solve for the
upper quantiles, yielding the quantile function

$$
Q(p) = (c/(1-p))^{1/\alpha}.
$$

Since the order statistics estimate quantiles, we have

$$
X_{(j)} \approx (c / (1 - j/(n+1)))^{1/\alpha}.
$$

An estimator known as the *Hill estimator* begins by considering
the ratios of upper order statistics

$$
X_{(n-j)} / X_{(n-k)} \approx ((j+1)/(k+1))^{1/\alpha}.
$$

Thus

$$
\log X_{(n-j)} / X_{(n-k)} \approx -\alpha^{-1} \log((k+1)/(j+1)).
$$

If we hold $k$ fixed and average these ratios for $1 \le j < k$,
we get

$$
\hat{A} \equiv (k-1)^{-1}\sum_{j=1}^{k-1} \log X_{(n-j)} / X_{(n-k)} \approx -\alpha^{-1} \sum_{j=1}^{k-1}\log((k+1)/(j+1)).
$$

This establishes a relationship between the statistic $\hat{A}$ and the quantity of
intetrest $\alpha$.  The constant of proportionality turns out to be nearly
equal to $1$:

$$
\sum_{j=1}^{k-1}\log((k+1)/(j+1)) \approx \int_0^k \log((k+1)/(x+1))dx \rightarrow -1.
$$

Therefore, the *Hill estimate* of the tail parameter is

$$
\hat{\alpha}_{\rm Hill} = 1/\hat{A}.
$$

In the Hill estimate, the value of $k$ is a tuning parameter.  To select
$k$, we usually calculate $\hat{\alpha}_{\rm Hill}$ for various values
of $k$ (typically $k \approx 20-200$) and choose a value that corresponds
to a stable range of values of the estimate.

### The generalized Pareto distribution

As noted above, the one-parameter Pareto distribution may not fit many data
sets well, and further has an awkward sample space of $[1, infty)$.  To
address these issues, the
[generalized Pareto distribution](https://en.wikipedia.org/wiki/Generalized_Pareto_distribution)
was developed, which has sample space $[0, \infty)$ and complementary CDF

$$
P(X > t) = \sigma^{-1}(1 + t/(\sigma\alpha))^{-\alpha}.
$$

Note that the ratio of the CCDF's for the Pareto and generalized Pareto
satisfies

$$
\frac{\sigma^{-1}(1 + t/(\sigma\alpha))^{-\alpha}}{1/t^{\alpha}} \rightarrow (\sigma\alpha)^\alpha \ne 0,
$$

the generalized Pareto and Pareto distributions both have power law tails with the
same tail parameter $\alpha$.

When using the generalized Pareto distribution, we may be able to model the population
of interest directly rather than converting the data to exceedances.

### Maximum likelihood estimation

### Return levels

The *m-observation return level* is a value $x$ that is expected to be
exceeded once out of every $m$ observations.  If $F(x)$ is the cumulative
distribution function (CDF) of a random variable $X$, then the m-observation return level is
the solution to $F(x) = 1 - 1/m$.  In other words, the m-observation return level
is the $1 - 1/m$ quantile of $X$.

Suppose that $X$ is a random variable that may have heavy
tails.  We are interested in m-observation return levels of $X$ but do not
wish to model the distribution of $X$, instead focusing only on the right tail.
Therefore, we consider $P(X-T | X>T)$. The m-observation return level of $X | X>T$
is $T+X^*$, where $F(X^*) = 1 - 1/m$, with $F$ being the CDF of $X-T | X>T$.
Note that $X^*$ is the $1 - 1/m$ quantile of $X-T | X>T$.

If we want the m-observation return level of the marginal distribution of $X$,
let $q = P(X>T)$, and
let $F(\tilde{X}) = 1 - 1/(q\cdot m)$.  Then $T+\tilde{X}$ is the desired
m-observation return level.

## L-moments

A classical statistical moment is defined to be the expected value of
a random variable raised to a power.  For example, $E[X]$ is the
expected value and $E[X^2]$ is the raw second moment.  In practice we
usually work with the *centered moments*, for example $E[(X-EX)^2]$ is
the centered second moment, which is better known as the variance.

In principle, if you know all the moments of a distribution, then you
know everything about the distribution (there are some technical
conditions for this to be literally true, as it is based on the
invertibility of the moment generating function).  But this fact is
not very useful in practice because it is nearly impossible to
estimate high order moments $E[(X-EX)^k]$ for large values of $k$.
The sample estimator of this moment is
$n^{-1}\sum (X_i - \bar{X})^k$,
and this estimator is consistent and asymptotically
unbiased, but if $k > 2$ it has huge mean squared error for practically realistic
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
$\lambda_3^s = \lambda_3/\lambda_2$ and
$\lambda_4^s = \lambda_4/\lambda_2$.
Note that these standardized L-moments are
*scale invariant* meaning that their value is not changed by scaling
the data.  All L-moments except for the first L-moment are
*translation invariant*, meaning that their values are not changed by
adding a constant to all data values.  Scale and translation invariance
(also known as *affine invariance*) are
important becuase they imply that the result does not depend on the
units or origin of the measurement scale.

L-moments are useful descriptive statistics that capture the shape
of distributions.  They are more robust (less sensitive to contamination)
than the classical moments, and one can estimate much higher order
L-moments than is practical with classical moments.

